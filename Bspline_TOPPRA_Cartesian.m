clear; clc; close all;

if ~exist('Link', 'class')
    error('请先安装MATLAB Robotics Toolbox');
end

% 定义机器人参数
a = [0, -0.42500, -0.39225, 0, 0, 0];
d = [0.089159, 0, 0, 0.10915, 0.09465, 0.08230];
alpha = [pi/2, 0, 0, pi/2, -pi/2, 0];

L(1) = Link('d', d(1), 'a', a(1), 'alpha', alpha(1), 'standard');
L(2) = Link('d', d(2), 'a', a(2), 'alpha', alpha(2), 'standard');
L(3) = Link('d', d(3), 'a', a(3), 'alpha', alpha(3), 'standard');
L(4) = Link('d', d(4), 'a', a(4), 'alpha', alpha(4), 'standard');
L(5) = Link('d', d(5), 'a', a(5), 'alpha', alpha(5), 'standard');
L(6) = Link('d', d(6), 'a', a(6), 'alpha', alpha(6), 'standard');

L(1).qlim = [-2*pi, 2*pi]; 
L(2).qlim = [-2*pi, 2*pi];
L(3).qlim = [-pi, pi];     
L(4).qlim = [-2*pi, 2*pi];
L(5).qlim = [-2*pi, 2*pi]; 
L(6).qlim = [-2*pi, 2*pi];
ur5 = SerialLink(L, 'name', 'UR5');

% 数据准备: 转换途径点 
fprintf('\n===== 1. 数据准备: 转换途径点 =====\n');

q_via_joint = [
    linspace(0, pi/2, 7);          
    linspace(-pi/3, -pi/6, 7);     
    linspace(pi/3, pi/4, 7);       
    linspace(-pi/2, -pi/3, 7);     
    linspace(pi/2, 2*pi/3, 7);     
    linspace(0, pi/3, 7)           
];
num_via_pts = size(q_via_joint, 2);

pos_via = zeros(3, num_via_pts);
quat_via = zeros(4, num_via_pts); 

fprintf('正在计算正运动学...\n');
for i = 1:num_via_pts
    T = ur5.fkine(q_via_joint(:, i)');
    pos_via(:, i) = T.t; 
    Q = UnitQuaternion(T); 
    quat_via(:, i) = double(Q)'; 
end
fprintf('笛卡尔途径点准备完毕。\n');

% 笛卡尔空间 B样条规划 
fprintf('\n===== 2. 笛卡尔空间 B样条规划 =====\n');

p = 5; 
num_ctrl_pts = 12;
lambda = 1e-4; 
epsilon = 1e-8;

dists = sqrt(sum(diff(pos_via, 1, 2).^2, 1));
u_bar = [0, cumsum(dists) / sum(dists)];

U = BuildKnotVector(p, num_ctrl_pts, u_bar, num_via_pts);
[C_mat, A_mat] = BuildSmoothMatrix(num_ctrl_pts, p, U);

% 规划位置
v0_pos = [0;0;0]; 
v6_pos = [0;0;0];
a0_pos = [0;0;0]; 
a6_pos = [0;0;0];
[B_pos, W_pos, Q_target_pos] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, pos_via, v0_pos, v6_pos, a0_pos, a6_pos, U);

M_pos = B_pos' * W_pos * B_pos + lambda * (C_mat' * A_mat * C_mat) + epsilon * eye(num_ctrl_pts);
F_pos = B_pos' * W_pos * Q_target_pos;
P_pos_ctrl = M_pos \ F_pos;
P_pos_ctrl_T = P_pos_ctrl'; 

% 规划姿态 
v0_quat = [0;0;0;0]; 
v6_quat = [0;0;0;0];
a0_quat = [0;0;0;0]; 
a6_quat = [0;0;0;0];
[B_quat, W_quat, Q_target_quat] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, quat_via, v0_quat, v6_quat, a0_quat, a6_quat, U);

M_quat = B_quat' * W_quat * B_quat + lambda * (C_mat' * A_mat * C_mat) + epsilon * eye(num_ctrl_pts);
F_quat = B_quat' * W_quat * Q_target_quat;
P_quat_ctrl = M_quat \ F_quat;
P_quat_ctrl_T = P_quat_ctrl'; 

fprintf('位置与姿态拟合完成。\n');

% TOPP-RA (集成关节速度约束) 
fprintf('\n===== 3. TOPP-RA (精确版可达性分析) =====\n');

% 物理约束设置 
v_max_lin = 1.5;    % 最大线速度 (m/s)
a_max_lin = 3.0;    % 最大线加速度 (m/s^2)
w_max = 2.5;        % 最大角速度 (rad/s)
alpha_max = 10.0;   % 最大角加速度 (rad/s^2)

limit_q_vel = 3.0;  % 关节最大速度限制 (rad/s)

% 离散化设置
N_grid = 2000; 
ds = 1 / (N_grid - 1);  
s_grid = linspace(0, 1, N_grid);

% 预计算几何导数
fprintf('预计算曲线几何导数...\n');
P_p = zeros(3, N_grid); 
P_pp = zeros(3, N_grid);
Q_p = zeros(4, N_grid); 
Q_pp = zeros(4, N_grid);

for i = 1:N_grid
    CK_pos = CurveDerivs1(s_grid(i), U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p(:, i) = CK_pos(2, :)';
    P_pp(:, i) = CK_pos(3, :)';
    
    CK_quat = CurveDerivs1(s_grid(i), U, p, P_quat_ctrl_T, 2, num_ctrl_pts - 1);
    Q_p(:, i) = CK_quat(2, :)';
    Q_pp(:, i) = CK_quat(3, :)';
end

% 计算 MVC (包含笛卡尔与关节约束)
fprintf('计算 MVC (笛卡尔 + 关节 Jacobian 映射)...\n');
x_limit = zeros(1, N_grid); 

% IK 热启动变量
last_q_ik = q_via_joint(:,1)'; 

for i = 1:N_grid
    % A. 笛卡尔约束
    
    % 线速度: ||P'||^2 * x <= v_max^2
    norm_Pp_sq = sum(P_p(:, i).^2); 
    if norm_Pp_sq > 1e-12
        x_v = (v_max_lin^2) / norm_Pp_sq; 
    else
        x_v = 1e10; 
    end

    % 角速度: ||2*Q'||^2 * x <= w_max^2
    norm_Qp_sq = 4 * sum(Q_p(:, i).^2);
    if norm_Qp_sq > 1e-12
        x_w = (w_max^2) / norm_Qp_sq; 
    else
        x_w = 1e10;
    end

    % 加速度(离心项)约束
    norm_Ppp = norm(P_pp(:, i));
    if norm_Ppp > 1e-12
        x_a_lin = a_max_lin / norm_Ppp;
    else
        x_a_lin = 1e10; 
    end
    
    norm_Qpp = 2 * norm(Q_pp(:, i));
    if norm_Qpp > 1e-12
        x_a_rot = alpha_max / norm_Qpp;
    else
        x_a_rot = 1e10; 
    end
    
    x_cart_limit = min([x_v, x_w, x_a_lin, x_a_rot]);
    
    % B. 关节速度约束 
    % 计算当前笛卡尔位姿
    P_curr = CurvePoint(num_ctrl_pts-1, p, U, P_pos_ctrl_T, s_grid(i));
    Q_curr = CurvePoint(num_ctrl_pts-1, p, U, P_quat_ctrl_T, s_grid(i));
    Q_curr = Q_curr / norm(Q_curr); 
    
    T_curr = eye(4);
    T_curr(1:3,1:3) = quat2rotm(Q_curr');
    T_curr(1:3,4) = P_curr;
    
    % 快速逆解 
    try
        q_sol = ur5.ikine(T_curr, 'q0', last_q_ik, 'mask', [1 1 1 1 1 1], 'tol', 1e-4, 'maxiter', 20);
        last_q_ik = q_sol;
    catch
        q_sol = last_q_ik; 
    end
    
    % 3. 计算雅可比
    J = ur5.jacob0(q_sol);
    
    % 4. 计算几何 Twist (v_geom, w_geom)
    % v_geom = dP/ds = P_p
    % w_geom = 2 * dQ/ds * Q_inv
    v_vec = P_p(:, i);
    
    % 计算角速度向量 (严格公式: w = 2 * Q_dot * Q_conj)
    % Q = [w, x, y, z], MATLAB Quaternion format [w, v]
    % 注意: Q_p 和 Q_curr 都是 [w; x; y; z] 列向量
    Q_val = Q_curr;     
    Q_dot = Q_p(:, i);  
    
    % Q_conj = [w; -x; -y; -z]
    % 四元数乘积的虚部 = w1*v2 + w2*v1 + cross(v1, v2)
    % w_vec = 2 * Im(Q_dot * Q_conj)
    w_vec = 2 * (Q_dot(1)*(-Q_val(2:4)) + Q_val(1)*Q_dot(2:4) + cross(Q_dot(2:4), -Q_val(2:4)));
    
    twist_geom = [v_vec; w_vec];
    
    % 5. 映射到关节几何速度: q_dot_geom = J \ twist
    q_dot_geom = pinv(J) * twist_geom;
    
    % 6. 计算关节限制
    % |q_dot_geom_j| * sqrt(x) <= limit_q_vel
    % x <= (limit / |q_dot_geom_j|)^2
    x_joint_limit = 1e10;
    for j = 1:6
        geom_val = abs(q_dot_geom(j));
        if geom_val > 1e-6
            this_limit = (limit_q_vel / geom_val)^2;
            if this_limit < x_joint_limit
                x_joint_limit = this_limit;
            end
        end
    end
    
    % C. 综合限制 
    x_limit(i) = min(x_cart_limit, x_joint_limit);
end

x_limit(1) = 0; 
x_limit(end) = 0;

% Backward & Forward Pass 
beta = x_limit; 
fprintf('执行反向可达性分析...\n');
for i = N_grid-1 : -1 : 1
    x_next_max = beta(i+1);
    [acc_min, ~] = GetSDDotBounds_Exact(P_p(:,i), P_pp(:,i), Q_p(:,i), Q_pp(:,i), x_next_max, a_max_lin, alpha_max);
    x_curr_max = max(0, x_next_max - 2 * ds * acc_min);
    if x_curr_max < beta(i)
        beta(i) = x_curr_max; 
    end
end

fprintf('执行正向最优积分...\n');
x_opt = zeros(1, N_grid);
for i = 1 : N_grid-1
    x_curr = x_opt(i);
    [~, acc_max] = GetSDDotBounds_Exact(P_p(:,i), P_pp(:,i), Q_p(:,i), Q_pp(:,i), x_curr, a_max_lin, alpha_max);
    x_next_potential = max(0, x_curr + 2 * ds * acc_max);
    x_opt(i+1) = min(x_next_potential, beta(i+1));
end

x_opt = max(0, x_opt); 
sd_optimal = sqrt(x_opt);
x_opt_grad = gradient(x_opt) / ds; 
s_ddot_optimal = 0.5 * x_opt_grad; 

% 恢复时间戳
t_grid = zeros(1, N_grid);
for i = 1:N_grid-1
    v_avg = (sd_optimal(i) + sd_optimal(i+1)) / 2;
    if v_avg < 1e-4
        norm_P = norm(P_p(:,i)) + 1e-6; 
        eff_a = a_max_lin / norm_P; 
        dt_step = sqrt(2 * ds / eff_a);
    else
        dt_step = ds / v_avg;
    end
    t_grid(i+1) = t_grid(i) + max(1e-8, dt_step);
end
t_grid = real(t_grid);
T_total = t_grid(end);
fprintf('TOPP-RA 规划完成，总时间: %.4f 秒\n', T_total);

% 生成轨迹与逆解 
fprintf('\n===== 4. 生成轨迹与逆解 =====\n');

dt_sample = 0.02; 
t_traj = 0 : dt_sample : T_total;
num_steps = length(t_traj);

q_traj_ik = zeros(num_steps, 6);
pos_traj_cart = zeros(num_steps, 3); 
vel_traj_cart = zeros(num_steps, 3);
acc_traj_cart = zeros(num_steps, 3);

q_guess = q_via_joint(:, 1)'; 

for k = 1:num_steps
    t_now = t_traj(k);
    
    % 1. 时间映射
    if t_now > t_grid(end)
        t_now = t_grid(end); 
    end
    s_now = interp1(t_grid, s_grid, t_now, 'pchip');
    s_now = max(0, min(1, s_now)); 
    
    % 2. 状态映射
    s_dot_now = interp1(s_grid, sd_optimal, s_now, 'linear');
    if s_now >= 1 || s_now <= 0
        s_ddot_now = 0;
    else
        s_ddot_now = interp1(s_grid, s_ddot_optimal, s_now, 'linear');
    end
    
    % 3. B样条计算
    CK = CurveDerivs1(s_now, U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_now   = CK(1, :)'; 
    qs_now  = CK(2, :)'; 
    qss_now = CK(3, :)'; 
    
    pos_traj_cart(k, :) = P_now';
    vel_traj_cart(k, :) = (qs_now * s_dot_now)';
    acc_traj_cart(k, :) = (qs_now * s_ddot_now + qss_now * s_dot_now^2)';
    
    % 4. 逆解
    Q_raw = CurvePoint(num_ctrl_pts-1, p, U, P_quat_ctrl_T, s_now);
    Q_norm = Q_raw / norm(Q_raw); 
    Quat_obj = UnitQuaternion(Q_norm'); 
    T_des = SE3(P_now) * Quat_obj.SE3;
    
    try
        q_sol = ur5.ikine(T_des, 'q0', q_guess, 'mask', [1 1 1 1 1 1], 'tol', 1e-4);
    catch
        q_sol = q_guess;
    end
    
    % 关节连续性与限位
    if isempty(q_sol)
        q_sol = q_guess; 
    end
    % 防止关节跳变 (选离上一步最近的解)
    diff_q = q_sol - q_guess;
    for j=1:6
        if diff_q(j) > pi
            q_sol(j) = q_sol(j) - 2*pi; 
        elseif diff_q(j) < -pi
            q_sol(j) = q_sol(j) + 2*pi; 
        end
    end
    
    q_traj_ik(k, :) = q_sol;
    q_guess = q_sol;
end

% 验证最大关节速度
q_vel_check = zeros(num_steps, 6);
for k = 2:num_steps
    q_vel_check(k, :) = (q_traj_ik(k, :) - q_traj_ik(k-1, :)) / dt_sample;
end
max_real_q_vel = max(max(abs(q_vel_check)));
fprintf('最终生成的最大关节速度: %.3f rad/s (限制: %.1f)\n', max_real_q_vel, limit_q_vel);
if max_real_q_vel > limit_q_vel * 1.1
    warning('虽然加上了MVC约束，但由于数值积分误差，仍有轻微超速。');
end

% 绘图与分析 
fprintf('正在绘制分析图...\n');

vel_mag = sqrt(sum(vel_traj_cart.^2, 2));
acc_mag = sqrt(sum(acc_traj_cart.^2, 2));

figure('Position', [100, 100, 1000, 800], 'Color', 'w');
sgtitle('TOPP-RA 规划结果分析');

subplot(3, 1, 1);
plot(t_traj, pos_traj_cart, 'LineWidth', 1.5); 
grid on;
legend('X','Y','Z');
title('末端位置');

subplot(3, 1, 2);
plot(t_traj, vel_mag, 'k-', 'LineWidth', 2);
hold on;
yline(v_max_lin, 'r--', 'LineWidth', 2);
grid on; 
ylabel('速度 (m/s)');
title(sprintf('线速度 (Limit: %.1f)', v_max_lin));

subplot(3, 1, 3);
plot(t_traj, acc_mag, 'b-', 'LineWidth', 2);
hold on;
yline(a_max_lin, 'r--', 'LineWidth', 2);
grid on;
ylabel('加速度 (m/s^2)'); 
title(sprintf('线加速度 (Limit: %.1f)', a_max_lin));

% 动画
figure('Position', [50, 50, 1200, 600], 'Color', 'w');
subplot(1, 2, 1);
hold on; 
grid on; 
axis equal;
view(3);
xlabel('X'); 
ylabel('Y'); 
zlabel('Z'); 
title('笛卡尔轨迹');
plot3(pos_via(1,:), pos_via(2,:), pos_via(3,:), 'ro', 'MarkerFaceColor', 'r');
plot3(pos_traj_cart(:,1), pos_traj_cart(:,2), pos_traj_cart(:,3), 'b-', 'LineWidth', 2);

subplot(1, 2, 2);
ur5.plot(q_traj_ik(1,:));
title('UR5 动画演示');
hold on;
plot3(pos_traj_cart(:,1), pos_traj_cart(:,2), pos_traj_cart(:,3), 'b-', 'LineWidth', 2);

step_size_anim = 2; 
for k = 1:step_size_anim:num_steps
    ur5.animate(q_traj_ik(k,:));
    drawnow;
end

% 辅助函数 
function [alpha_min, alpha_max] = GetSDDotBounds_Exact(Pp, Ppp, Qp, Qpp, x, a_max, wa_max)
    alpha_min = -1e10; 
    alpha_max = 1e10;
    
    % 线加速度
    A_vec = Pp; 
    B_vec = Ppp * x;
    R = a_max;
    a_c = dot(A_vec, A_vec); 
    b_c = 2 * dot(A_vec, B_vec); 
    c_c = dot(B_vec, B_vec) - R^2;
    [min1, max1] = SolveQuadraticInequality(a_c, b_c, c_c);

    alpha_min = max(alpha_min, min1); 
    alpha_max = min(alpha_max, max1);
    
    % 角加速度
    A_vec_w = 2 * Qp; B_vec_w = 2 * Qpp * x; R_w = wa_max;
    a_cw = dot(A_vec_w, A_vec_w); b_cw = 2 * dot(A_vec_w, B_vec_w); c_cw = dot(B_vec_w, B_vec_w) - R_w^2;
    [min2, max2] = SolveQuadraticInequality(a_cw, b_cw, c_cw);
    alpha_min = max(alpha_min, min2); alpha_max = min(alpha_max, max2);
    
    if alpha_min > alpha_max
        alpha_max = alpha_min; 
    end
end

function [xmin, xmax] = SolveQuadraticInequality(a, b, c)
    if abs(a) < 1e-10
        if abs(b) < 1e-10
            if c <= 1e-10 % 稍微放宽一点容差
                xmin = -1e10; 
                xmax = 1e10; % 恒成立
            else
                xmin = 1e10; 
                xmax = -1e10; % 无解
            end
        elseif b > 0
            % x <= -c/b
            xmin = -1e10;  % 负无穷
            xmax = -c/b;
        else
            % x >= -c/b
            xmin = -c/b; 
            xmax = 1e10;   % 正无穷
        end
    else
        % 标准二次不等式
        delta = b^2 - 4*a*c;
        
        % 如果 delta 是微小的负数，强制视为 0
        if delta < 0 && delta > -1e-8
            delta = 0;
        end
        
        if delta < 0
            % 无实根
            if a > 0
                % 开口向上，最小值 > 0 -> 无解
                xmin = 1e10; 
                xmax = -1e10;
            else
                % 开口向下，最大值 < 0 -> 恒成立
                xmin = -1e10; 
                xmax = 1e10;
            end
        else
            sqrt_delta = sqrt(delta);
            r1 = (-b - sqrt_delta) / (2*a);
            r2 = (-b + sqrt_delta) / (2*a);
            
            if a > 0
                % 开口向上，两根之间 [r1, r2] (假设 r1 < r2)
                xmin = min(r1, r2);
                xmax = max(r1, r2);
            else
                % 开口向下，两根之外 (-inf, r1] U [r2, inf)
                % 物理约束下通常 a = ||P'||^2 > 0，这里仅处理理论情况
                xmin = -1e10; xmax = 1e10; 
            end
        end
    end
end

function U = BuildKnotVector(p, num_ctrl_pts, u_bar, num_via_pts)
    U = zeros(1, num_ctrl_pts + p + 1);
    U(1:p+1) = 0; 
    U(end-p:end) = 1;
    num_internal = num_ctrl_pts - p - 1;
    if num_internal > 0
        step = (num_via_pts - 1) / (num_internal + 1);
        for i = 1:num_internal
            idx = 1 + i * step;
            idx_floor = floor(idx); 
            idx_ceil = ceil(idx);
            if idx_floor == idx_ceil
                idx_floor = min(idx_floor, num_via_pts);
                U(p+1+i) = u_bar(idx_floor);
            else
                w = idx - idx_floor;
                U(p+1+i) = (1 - w) * u_bar(idx_floor) + w * u_bar(idx_ceil);
            end
        end
    end
end

function [B, W, Q_target] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, q_via, v0, v6, a0, a6, U)
    dim = size(q_via, 1); 
    num_via = size(q_via, 2);
    num_targets = num_via + 4;
    Q_target = zeros(num_targets, dim);
    Q_target(1:num_via, :) = q_via';
    Q_target(num_via+1, :) = v0'; Q_target(num_via+2, :) = a0';
    Q_target(num_via+3, :) = v6'; Q_target(num_via+4, :) = a6';
    
    W = eye(num_targets);
    W(1,1) = 1e6; W(num_via, num_via) = 1e6;
    W(num_via+1:end, num_via+1:end) = 1e6 * eye(4);
    
    B = zeros(num_targets, num_ctrl_pts);
    for k = 1:num_via
        u_k = u_bar(k);
        span = FindSpan(num_ctrl_pts-1, p, u_k, U);
        N = BasisFuns(span, u_k, p, U);
        for j=0:p
            B(k, span-p+j+1) = N(j+1);
        end
    end
    
    % 起点一阶导数约束 (速度)
    span0 = FindSpan(num_ctrl_pts-1, p, u_bar(1), U);
    ders0 = DersBasisFuns(span0, u_bar(1), p, 2, U);
    for j=0:p
        B(num_via+1, span0-p+j+1) = ders0(2, j+1); 
        B(num_via+2, span0-p+j+1) = ders0(3, j+1); % 起点二阶导
    end
    
    % 终点一阶导数约束 (速度)
    span6 = FindSpan(num_ctrl_pts-1, p, u_bar(end), U);
    ders6 = DersBasisFuns(span6, u_bar(end), p, 2, U);
    for j=0:p
        B(num_via+3, span6-p+j+1) = ders6(2, j+1); 
        B(num_via+4, span6-p+j+1) = ders6(3, j+1); % 终点二阶导
    end
end

function [C, A] = BuildSmoothMatrix(num_ctrl_pts, p, U)
    % 构建平滑矩阵 (用于最小化加速度变化)
    n = num_ctrl_pts;
    C = zeros(n - 2, n);
    for k = 1:(n - 2)
        u_ij1 = U(k + p + 1) - U(k + 1);
        u_ij2 = U(k + p + 2) - U(k + 2);
        if abs(u_ij1) > 1e-10
            C(k, k) = 6 / u_ij1; 
        end

        if abs(u_ij1) > 1e-10 && abs(u_ij2) > 1e-10
            C(k, k+1) = -6 * (1/u_ij1 + 1/u_ij2); 
        end

        if abs(u_ij2) > 1e-10
            C(k, k+2) = 6 / u_ij2; 
        end
    end
    A = zeros(n-2, n-2);
    u_term = U(n+1) - U(p+1);
    for k = 1:(n-2)
        A(k,k) = 2 * (U(k+p+2) - U(k+p+1)) * u_term;
        if k < (n-2)
            A(k, k+1) = (U(k+p+2) - U(k+p+1)) * u_term;
            A(k+1, k) = A(k, k+1);
        end
    end
end

function span = FindSpan(n, p, u, U)
    % 寻找节点区间索引
    if u >= U(n+2)
        span = n; 
        return; 
    end

    if u <= U(p+1)
        span = p;
        return; 
    end

    low = p; 
    high = n + 1; 
    mid = floor((low + high) / 2);
    while u < U(mid+1) || u >= U(mid+2)
        if u < U(mid+1)
            high = mid; 
        else
            low = mid; 
        end
        mid = floor((low + high) / 2);
    end
    span = mid;
end

function N = BasisFuns(i, u, p, U)
    % 计算基函数值
    N = zeros(p+1, 1); 
    N(1) = 1.0;
    left = zeros(p+1, 1); 
    right = zeros(p+1, 1);
    for j = 1:p
        left(j+1) = u - U(i+1-j+1); 
        right(j+1) = U(i+j+1) - u;
        saved = 0.0;
        for r = 0:(j-1)
            temp = N(r+1) / (right(r+2) + left(j-r+1));
            N(r+1) = saved + right(r+2) * temp;
            saved = left(j-r+1) * temp;
        end
        N(j+1) = saved;
    end
end

function ders = DersBasisFuns(i, u, p, n, U)
    % 计算基函数导数
    ders = zeros(n+1, p+1);
    ndu = zeros(p+1, p+1); 
    ndu(1,1) = 1.0;
    left = zeros(p+1, 1); 
    right = zeros(p+1, 1);
    for j = 1:p
        left(j+1) = u - U(i+1-j+1); right(j+1) = U(i+j+1) - u;
        saved = 0.0;
        for r = 0:(j-1)
            ndu(j+1, r+1) = right(r+2) + left(j-r+1);
            temp = ndu(r+1, j) / ndu(j+1, r+1);
            ndu(r+1, j+1) = saved + right(r+2) * temp;
            saved = left(j-r+1) * temp;
        end
        ndu(j+1, j+1) = saved;
    end
    for j = 0:p
        ders(1, j+1) = ndu(j+1, p+1); 
    end
    a = zeros(2, p+1);
    for r = 0:p
        s1 = 0; s2 = 1; a(1, 1) = 1.0;
        for k = 1:n
            d = 0.0; 
            rk = r - k; 
            pk = p - k;
            if r >= k
                a(s2+1, 1) = a(s1+1, 1) / ndu(pk+2, rk+1); 
                d = a(s2+1, 1) * ndu(rk+1, pk+1); 
            end

            if rk >= -1
                j1 = 1; 
            else
                j1 = -rk; 
            end

            if r-1 <= pk
                j2 = k - 1; 
            else
                j2 = p - r; 
            end

            for j = j1:j2
                a(s2+1, j+1) = (a(s1+1, j+1) - a(s1+1, j)) / ndu(pk+2, rk+j+1);
                d = d + a(s2+1, j+1) * ndu(rk+j+1, pk+1);
            end
            
            if r <= pk
                a(s2+1, k+1) = -a(s1+1, k) / ndu(pk+2, r+1); 
                d = d + a(s2+1, k+1) * ndu(r+1, pk+1); 
            end
            ders(k+1, r+1) = d; 
            j = s1; 
            s1 = s2; 
            s2 = j;
        end
    end
    r = p;
    for k = 1:n
        for j = 0:p
            ders(k+1, j+1) = ders(k+1, j+1) * r; 
        end
        r = r * (p - k);
    end
end

function C = CurvePoint(n, p, U, P, u)
    % 计算曲线上的点
    span = FindSpan(n, p, u, U);
    N = BasisFuns(span, u, p, U);
    C = zeros(size(P, 1), 1);
    for i = 0:p
        C = C + N(i+1) * P(:, span-p+i+1);
    end
end

function CK = CurveDerivs1(u, U, p, P, d, n)
    % 计算曲线的一阶、二阶导数
    dim = size(P, 1); du = min(d, p);
    CK = zeros(d+1, dim);
    span = FindSpan(n, p, u, U);
    nders = DersBasisFuns(span, u, p, du, U);
    for k = 0:du
        for j = 0:p
            CK(k+1, :) = CK(k+1, :) + nders(k+1, j+1) * P(:, span-p+j+1)';
        end
    end
end
