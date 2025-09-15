clear; clc; close all;

if ~exist('Link', 'class')
    error('请先安装MATLAB Robotics Toolbox（Peter Corke版）和Robotics System Toolbox');
end

% UR5标准DH参数
a = [0, -0.42500, -0.39225, 0, 0, 0];
d = [0.089159, 0, 0, 0.10915, 0.09465, 0.08230];
alpha = [pi/2, 0, 0, pi/2, -pi/2, 0];

% 建立6轴UR5模型
L1 = Link('d', d(1), 'a', a(1), 'alpha', alpha(1), 'standard');
L2 = Link('d', d(2), 'a', a(2), 'alpha', alpha(2), 'standard');
L3 = Link('d', d(3), 'a', a(3), 'alpha', alpha(3), 'standard');
L4 = Link('d', d(4), 'a', a(4), 'alpha', alpha(4), 'standard');
L5 = Link('d', d(5), 'a', a(5), 'alpha', alpha(5), 'standard');
L6 = Link('d', d(6), 'a', a(6), 'alpha', alpha(6), 'standard');
tool_robot = SerialLink([L1, L2, L3, L4, L5, L6], 'name', 'UR5');
tool_robot.display();  % 打印机械臂DH参数，验证模型


%定义B样条轨迹的核心参数：控制点+时间参数
% 这里设计5个笛卡尔空间控制点（x,y,z），确保在UR5工作空间内（半径~850mm，z>0.1m）
% 控制点cpts：3行n列矩阵（3对应x/y/z，n为控制点数，这里n=5）
cpts = [
    0.3, 0.5, 0.6, 0.4, 0.2;   % x轴坐标（单位：m）
    0.1, 0.2, 0.0, -0.2, -0.1;  % y轴坐标（单位：m）
    0.4, 0.5, 0.6, 0.5, 0.4    % z轴坐标（单位：m）
];
n_cpts = size(cpts, 2);  % 控制点数：5

% 时间参数：定义每个控制点对应的时间
v_avg = 0.2;  % 末端平均速度（m/s，比原直线插补慢，突出B样条平滑性）
% 计算控制点之间的距离，估算总运动时间
dist_cpts = 0;
for i = 2:n_cpts
    dist_cpts = dist_cpts + norm(cpts(:,i) - cpts(:,i-1));
end
T_total = dist_cpts / v_avg;  % 总运动时间（s）：根据距离和速度估算
tpts = linspace(0, T_total, n_cpts);  % 控制点对应的时间（均匀分配）

% 输出轨迹的时间向量tvec：按固定步长dt生成，用于动画和逆解
dt = 0.05;  % 时间步长（s）
tvec = linspace(0, T_total, round(T_total/dt) + 1);  % 所有时间节点
N = length(tvec);  % 轨迹总步数（用于循环逆解）


%用bsplinepolytraj生成B样条位置轨迹（笛卡尔空间）
% 函数功能：输入控制点cpts和时间tpts，输出位置/速度/加速度轨迹
% [p_traj, pd_traj, pdd_traj, pp] = bsplinepolytraj(cpts, tpts, tvec)
% - p_traj：3×N矩阵，每个时间节点的末端位置（x/y/z）
% - pd_traj：3×N矩阵，末端速度
% - pdd_traj：3×N矩阵，末端加速度
% - pp：分段多项式对象（可用于后续分析）
[p_traj, pd_traj, pdd_traj, pp] = bsplinepolytraj(cpts, tpts, tvec);

% 可选：打印B样条轨迹信息
fprintf('B样条轨迹总步数：%d\n', N);
fprintf('总运动时间：%.2f s\n', T_total);
fprintf('末端平均速度：%.2f m/s\n', dist_cpts / T_total);

%逆解：将笛卡尔位置轨迹转换为关节角度轨迹
% 保持末端姿态不变：与第一个控制点的姿态一致（简化实现，后续可扩展姿态B样条）
% 步骤：1. 计算第一个控制点的姿态R_start；2. 对每个位置p_traj(:,k)构造齐次矩阵；3. 逆解求关节角
q_prev = zeros(1,6);  % 初始关节角：先求第一个控制点的逆解
% 构造第一个控制点的齐次矩阵T_start（位置cpts(:,1) + 姿态R_start）
T_start = forwardmotion(q_prev);  % 先用初始关节角求正解（后续替换为逆解结果）
R_start = T_start(1:3,1:3);  % 固定姿态：后续所有轨迹点都用这个旋转矩阵

% 初始化关节角度轨迹q_traj（N×6矩阵）
q_traj = zeros(N, 6);
for k = 1:N
    % 第k步的末端位置：p_traj(:,k)
    p_k = p_traj(:,k);
    
    % 构造第k步的齐次变换矩阵T_k（固定姿态R_start + 位置p_k）
    T_k = [
        R_start,       p_k(:);  % 3×3旋转矩阵 + 3×1位置向量
        0, 0, 0,       1        % 齐次矩阵第四行
    ];
    
    % 调用自定义逆解函数inversemotion，得到8组关节角解
    q_all = inversemotion(T_k);
    
    % 选择与前一步关节角最接近的解（避免关节突变，保证运动平滑）
    diffs = vecnorm(q_all - q_prev, 2, 2);  % 计算每组解与前一步的2范数（差异大小）
    [~, idx] = min(diffs);  % 找差异最小的解的索引
    q_pick = q_all(idx, :);  % 选中的关节角
    
    % 保存当前关节角，更新前一步关节角
    q_traj(k, :) = q_pick;
    q_prev = q_pick;
end


%动画可视化：验证B样条轨迹
figure('Position', [100, 100, 1000, 600]);  % 新建图形窗口
tool_robot.plot(q_traj(1, :));  % 初始姿态绘制
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('UR5 B样条曲线轨迹动画');

% 绘制B样条的「控制点」（黑色圆点）和「理论曲线」（灰色虚线），方便对比
plot3(cpts(1,:), cpts(2,:), cpts(3,:), 'ko', 'MarkerSize',5 , 'DisplayName', '控制点');
plot3(p_traj(1,:), p_traj(2,:), p_traj(3,:), 'k--', 'LineWidth', 1, 'DisplayName', '理论B样条');

% 更新机械臂姿态+绘制实际末端轨迹
for k = 1:N
    tool_robot.animate(q_traj(k, :));
    
    %绘制累积的末端实际轨迹
    plot3(p_traj(1, 1:k), p_traj(2, 1:k), p_traj(3, 1:k), 'r-', 'LineWidth', 2, 'DisplayName', '实际轨迹');
    drawnow;    
    % pause(0.01);
end

legend('Location', 'best');
hold off;


