clear; clc; close all;

%检查robotics toolbox工具箱
if ~exist('Link', 'class')
    error('请先安装MATLAB Robotics Toolbox');
end

%UR5 standard_DH parameter
a = [0,-0.42500,-0.39225,0,0,0];
d = [0.089159,0,0,0.10915,0.09465,0.08230];
alpha = [pi/2,0,0,pi/2,-pi/2,0];

% 建立UR5机械臂模型
L1 = Link('d', d(1),  'a', a(1), 'alpha', alpha(1),  'standard');
L2 = Link('d', d(2),  'a', a(2), 'alpha', alpha(2),  'standard');
L3 = Link('d', d(3),  'a', a(3), 'alpha', alpha(3),  'standard');
L4 = Link('d', d(4),  'a', a(4), 'alpha', alpha(4),  'standard');
L5 = Link('d', d(5),  'a', a(5), 'alpha', alpha(5),  'standard');
L6 = Link('d', d(6),  'a', a(6), 'alpha', alpha(6),  'standard');
tool_robot = SerialLink([L1,L2,L3,L4,L5,L6], 'name', 'UR5');
tool_robot.display();

%起始关节角度
q0 = [0, 0, 0, 0, 0, 0]; %0°，-30°，-15°，0°，30°，0°

%最终关节角度
qf = [pi, -pi/2, pi/3, -pi/3, -pi/6, pi/4]; %45°，-60°，-30°，30°，15°，45°

t_all = 2;        %总运动时间2s
dt = 0.01;        %时间步长
t = 0 : dt : t_all;   %201个时间点

joints_dof = 6;
q_t = zeros(length(t), joints_dof);      %每个时间点的每个关节角状态 （行为时间点，列为关节角）
v_t = zeros(length(t), joints_dof);      %关节角速度状态
a_t = zeros(length(t), joints_dof);

for i = 1:joints_dof    %遍历每个关节角
    %每个关节的起始值/目标值
    q0_i = q0(i);
    qf_i = qf(i);

    %三次多项式系数
    a0 = q0_i;
    a1 = 0;
    a2 = 0;
    a3 = 10 * (qf_i - q0_i) / (t_all^3);
    a4 = -15 * (qf_i - q0_i) / (t_all^4);
    a5 = 6 * (qf_i - q0_i) / (t_all^5);

    for j = 1:length(t)
        t_j = t(j);
        q_t(j, i) = a0 + a1 * t_j + a2 * (t_j^2) + a3 * (t_j^3) + a4 * (t_j^4) + a5 * (t_j^5);   %角度公式
        v_t(j, i) = a1 + 2 * a2 * t_j + 3 * a3 * (t_j^2) + 4 * a4 * (t_j^3) + 5 * a5 * (t_j^4);          %速度公式
        a_t(j, i) = 2 * a2 + 6 * a3 * t_j + 12 * a4 * (t_j^2) + 20 * a5 * (t_j^3);
    end
end

% 初始化绘图
figure('Position', [100, 100, 1000, 600]);
tool_robot.plot(q_t(1,:), 'fps', 1/dt);  % 只绘制一次
title('UR5机械臂五次多项式插值轨迹规划仿真', 'FontSize', 14);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
grid on; hold on;

% 循环更新关节角
for j = 1:length(t)
    tool_robot.animate(q_t(j,:));  % 用 animate 只更新姿态
    % plot3(p_traj(1:k,1), p_traj(1:k,2), p_traj(1:k,3), 'r-', 'LineWidth', 2);
    drawnow;
end

figure('Position', [200, 200, 1200, 800]);  % 1200×800窗口，用于显示曲线
for i = 1:joints_dof  % 遍历每个关节

    % 子图1：关节角随时间变化曲线
    subplot(3, joints_dof, i);  % 划分2行6列的子图，当前绘制第i个（第一行）
    plot(t, q_t(:, i) * 180 / pi, 'LineWidth', 1.5);  % 绘制时间t vs 关节角
    xlabel('时间 t (s)');  % x轴：时间
    ylabel(['关节', num2str(i), '角度 (°)']);  % y轴：关节i的角度
    title(['关节', num2str(i), '角度变化']);  % 子图标题
    grid on; hold on;  % 显示网格，保持当前图形
    % 绘制起始角度参考线（红色虚线）和目标角度参考线（绿色虚线）
    plot([0, t_all], [q0(i) * 180/pi, q0(i) * 180/pi], 'r--', 'LineWidth', 1);
    plot([0, t_all], [qf(i) * 180/pi, qf(i) * 180/pi], 'g--', 'LineWidth', 1);
    
    % 子图2：关节角速度随时间变化曲线
    subplot(3, joints_dof, i + joints_dof);  % 绘制第i+6个（第二行）
    plot(t, v_t(:, i) * 180 / pi, 'LineWidth', 1.5, 'Color', [0.2, 0.6, 0.8]);  % 角速度（弧度/s转角度/s）
    xlabel('时间 t (s)');  % x轴：时间
    ylabel(['关节', num2str(i), '角速度 (°/s)']);  % y轴：关节i的角速度
    title(['关节', num2str(i), '角速度变化']);  % 子图标题
    grid on; hold on;
    plot([0, t_all], [0, 0], 'r--', 'LineWidth', 1);  % 绘制零参考线（验证起始和结束时速度为0）

    % 子图3：关节角加速度随时间变化曲线
    subplot(3, joints_dof, i + 2 * joints_dof);  % 绘制第i+12个（第三行）
    plot(t, a_t(:, i) * 180 / pi, 'LineWidth', 1.5, 'Color', [0.2, 0.6, 0.8]);  % 角加速度
    xlabel('时间 t (s)');  % x轴：时间
    ylabel(['关节', num2str(i), '角加速度 (°/s²)']);  % y轴：关节i的角加速度
    title(['关节', num2str(i), '角加速度变化']);  % 子图标题
    grid on; hold on;
    plot([0, t_all], [0, 0], 'r--', 'LineWidth', 1);  % 绘制零参考线（验证起始和结束时速度为0）
end

sgtitle('UR5机械臂五次多项式插值轨迹规划结果', 'FontSize', 16);  % 总标题