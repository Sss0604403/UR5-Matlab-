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
q0 = [0, -pi/4, 0, -pi/2, 0, pi/4]; %0°，-30°，-15°，0°，30°，0°

%最终关节角度
qf = [pi/3, -pi/6, -pi/6, -pi/3, 0, pi/6]; %45°，-60°，-30°，30°，15°，45°

T_start = forwardmotion(q0);
T_end   = forwardmotion(qf);
% T_start = tool_robot.fkine(q0);  % 基座到TCP的正解（4x4齐次矩阵）
% T_end   = tool_robot.fkine(qf);

p_start = T_start(1:3,4);
p_end = T_end(1:3,4);
R_start = T_start(1:3,1:3);

L = norm(p_end - p_start);
v = 0.5;
dt = 0.05;
d = dt * v;
N = round(L / d) + 1;

delta_x = (p_end(1) - p_start(1)) / N;
delta_y = (p_end(2) - p_start(2)) / N;
delta_z = (p_end(3) - p_start(3)) / N;

p_traj = zeros(N,3);
q_traj = zeros(N,6);
q_prev = q0;
for k = 1:N
    p_traj(k,1) = p_start(1) + (k-1) * delta_x;
    p_traj(k,2) = p_start(2) + (k-1) * delta_y;
    p_traj(k,3) = p_start(3) + (k-1) * delta_z;

    T_k = [R_start, p_traj(k,:)';   % 3x3旋转矩阵 + 3x1位置列向量
       0, 0, 0, 1];
    q_all = inversemotion(T_k);
    diffs = vecnorm(q_all - q_prev,2,2);
    [~, idx] = min(diffs);
    q_pick = q_all(idx,:);
    
    q_traj(k,:) = q_pick;
    q_prev = q_pick;
end

figure;
tool_robot.plot(q_traj(1,:));
hold on; grid on;
plot3(p_traj(:,1), p_traj(:,2), p_traj(:,3), 'r-', 'LineWidth', 2);

for k = 1:N
    tool_robot.animate(q_traj(k,:));
    plot3(p_traj(1:k,1), p_traj(1:k,2), p_traj(1:k,3), 'r-', 'LineWidth', 2);
    drawnow;
end
