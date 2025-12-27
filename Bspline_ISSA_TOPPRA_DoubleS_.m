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

for i=1:6
    L(i).qlim = [-2*pi, 2*pi]; 
end
L(3).qlim = [-pi, pi];     
ur5 = SerialLink(L, 'name', 'UR5');

% B样条拟合
fprintf('\n===== 1. 数据准备: B样条拟合 =====\n');
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
for i = 1:num_via_pts
    T = ur5.fkine(q_via_joint(:, i)');
    pos_via(:, i) = T.t; 
    Q = UnitQuaternion(T); 
    quat_via(:, i) = double(Q)'; 
end

% B样条参数
p = 5; 
num_ctrl_pts = 12;
lambda = 1e-4; 
epsilon = 1e-8;

dists = sqrt(sum(diff(pos_via, 1, 2).^2, 1));
u_bar = [0, cumsum(dists) / sum(dists)];
U = BuildKnotVector(p, num_ctrl_pts, u_bar, num_via_pts);
[C_mat, A_mat] = BuildSmoothMatrix(num_ctrl_pts, p, U);

% 位置拟合
[B_pos, W_pos, Q_target_pos] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, pos_via,...
                                                         [0;0;0], [0;0;0], [0;0;0], [0;0;0], U);
M_pos = B_pos' * W_pos * B_pos + lambda * (C_mat' * A_mat * C_mat) + epsilon * eye(num_ctrl_pts);
F_pos = B_pos' * W_pos * Q_target_pos;
P_pos_ctrl_T = (M_pos \ F_pos)';

% 姿态拟合
[B_quat, W_quat, Q_target_quat] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, quat_via,...
                                                            [0;0;0;0], [0;0;0;0], [0;0;0;0], [0;0;0;0], U);
M_quat = B_quat' * W_quat * B_quat + lambda * (C_mat' * A_mat * C_mat) + epsilon * eye(num_ctrl_pts);
F_quat = B_quat' * W_quat * Q_target_quat;
P_quat_ctrl_T = (M_quat \ F_quat)';

fprintf('B样条拟合完成。\n');

% 参数设置 
v_max_lin = 1.5;    % m/s
a_max_lin = 3.0;    % m/s^2
w_max = 2.5;       
alpha_max = 10.0;   
j_max = 50.0;       
limit_q_vel = 3.0;  % 关节最大速度限制 (rad/s)

% 方法一: TOPP-RA 
fprintf('\n===== 2. 方法一: TOPP-RA 计算 =====\n');
tic;
N_grid = 2000; % 稍微降低以平衡速度         
ds = 1 / (N_grid - 1);  
s_grid = linspace(0, 1, N_grid); % 路径等分的临时网格

P_p = zeros(3, N_grid); 
P_pp = zeros(3, N_grid);
Q_p = zeros(4, N_grid); 
Q_pp = zeros(4, N_grid);

% 纯路径几何参数
for i = 1:N_grid
    CK_pos = CurveDerivs1(s_grid(i), U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p(:, i) = CK_pos(2, :)'; 
    P_pp(:, i) = CK_pos(3, :)';
    
    CK_quat = CurveDerivs1(s_grid(i), U, p, P_quat_ctrl_T, 2, num_ctrl_pts - 1);
    Q_p(:, i) = CK_quat(2, :)'; 
    Q_pp(:, i) = CK_quat(3, :)';
end

x_limit = zeros(1, N_grid); 
last_q_ik = q_via_joint(:,1)'; 

for i = 1:N_grid
    % 笛卡尔约束
    % v = dP/dt = (dP/ds) × (ds/dt) = P' × ṡ ≤ v_max
    % 线速度约束: ||P'||^2 * x <= v_max^2
    norm_Pp_sq = sum(P_p(:, i).^2); 
    if norm_Pp_sq > 1e-12
        x_v = (v_max_lin^2) / norm_Pp_sq; 
    else
        x_v = 1e10; 
    end

    % q̇ = dq/dt = (dq/ds) × ṡ = q' × ṡ
    % ω = 2 × q̇ ⊗ q⁻¹ ≤ ω_max
    % 角速度约束: ||2*Q'||^2 * x <= w_max^2        
    norm_Qp_sq = 4 * sum(Q_p(:, i).^2);
    if norm_Qp_sq > 1e-12
        x_w = (w_max^2) / norm_Qp_sq; 
    else
        x_w = 1e10; 
    end

    % 线加速度引起的限速 (离心力约束)
    % 物理原理: 即使切向不加速，离心加速度 ||P'' * x|| 也不能超过 a_max
    % a = dv/dt = P' × s̈ + P'' × ṡ² = P' × s̈ + P'' × x ≤ a_max
    % 即: x <= a_max / ||P''||
    norm_Ppp = norm(P_pp(:, i));
    if norm_Ppp > 1e-12
        x_a_lin = a_max_lin / norm_Ppp;
    else
        x_a_lin = 1e10; 
    end
    
    % 角加速度引起的限速 (旋转离心力约束)
    % α = dω/dt ≈ 2 × (q' × s̈ + q'' × x) ≤ α_max    
    norm_Qpp = 2 * norm(Q_pp(:, i));
    if norm_Qpp > 1e-12
        x_a_rot = alpha_max / norm_Qpp; 
    else
        x_a_rot = 1e10;
    end
    
    x_cart_limit = min([x_v, x_w, x_a_lin, x_a_rot]);
    
    % 关节速度约束 
    P_curr = CurvePoint(num_ctrl_pts-1, p, U, P_pos_ctrl_T, s_grid(i));
    Q_curr = CurvePoint(num_ctrl_pts-1, p, U, P_quat_ctrl_T, s_grid(i)); 
    Q_curr = Q_curr / norm(Q_curr); 

    T_curr = eye(4); 
    T_curr(1:3,1:3) = quat2rotm(Q_curr'); 
    T_curr(1:3,4) = P_curr;
    
    try
        q_sol = ur5.ikine(T_curr, 'q0', last_q_ik, 'mask', [1 1 1 1 1 1], 'tol', 1e-4, 'maxiter', 15);
        last_q_ik = q_sol;
    catch
        q_sol = last_q_ik; 
    end
    J = ur5.jacob0(q_sol);

    v_vec = P_p(:, i);
    Q_val = Q_curr; 
    Q_dot = Q_p(:, i);  
    w_vec = 2 * (Q_dot(1)*(-Q_val(2:4)) + Q_val(1)*Q_dot(2:4) + cross(Q_dot(2:4), -Q_val(2:4)));
    twist_geom = [v_vec; w_vec];  % 两个速度约束
    q_dot_geom = pinv(J) * twist_geom;
    
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
    
    x_limit(i) = min(x_cart_limit, x_joint_limit);
end
x_limit(1) = 0; 
x_limit(end) = 0;

% Backward Pass 求解可达集
beta = x_limit; 
for i = N_grid-1 : -1 : 1
    x_next_max = beta(i+1);
    [acc_min, ~] = GetSDDotBounds_Exact(P_p(:,i), P_pp(:,i), Q_p(:,i), Q_pp(:,i), x_next_max, a_max_lin, alpha_max);
    x_curr_max = max(0, x_next_max - 2 * ds * acc_min);

    if x_curr_max < beta(i)
        beta(i) = x_curr_max;
    end
end

% Forward Pass 求解可控集
x_opt = zeros(1, N_grid);
for i = 1 : N_grid-1
    x_curr = x_opt(i);
    [~, acc_max] = GetSDDotBounds_Exact(P_p(:,i), P_pp(:,i), Q_p(:,i), Q_pp(:,i), x_curr, a_max_lin, alpha_max);
    x_next_potential = max(0, x_curr + 2 * ds * acc_max);
    x_opt(i+1) = min(x_next_potential, beta(i+1));
end

x_opt = max(0, x_opt); 
sd_optimal = sqrt(x_opt);  % ds/dt
% x= ds/dt, 两边对s进行求导
% dx/ds = 2sdot * dsdot/ds, dsdot/ds = dsdot/dt * dt/ds = sddot/sdot
% 故dx/ds = 2sddot
x_opt_grad = gradient(x_opt) / ds; 
s_ddot_optimal = 0.5 * x_opt_grad;  % dds/t 

% 恢复时间信息
t_grid_toppra = zeros(1, N_grid);
for i = 1:N_grid-1
    v_avg = (sd_optimal(i) + sd_optimal(i+1)) / 2;
    if v_avg < 1e-4
        norm_P = norm(P_p(:,i)) + 1e-6;
        eff_a = a_max_lin / norm_P; 
        dt_step = sqrt(2 * ds / eff_a);
    else
        dt_step = ds / v_avg;
    end
    t_grid_toppra(i+1) = t_grid_toppra(i) + max(1e-8, dt_step);
end
T_toppra = t_grid_toppra(end);

fprintf('TOPP-RA 规划完成，总时间: %.4f 秒\n', T_toppra);
toc;

% 方法二: 分段 S-Curve 
fprintf('\n===== 3. 方法二: 分段 S-Curve ====\n');
tic;

safety_factor = 0.6; 
a_design_limit = safety_factor * a_max_lin; 

N_scan = 2000; % 提高扫描密度
u_scan = linspace(0, 1, N_scan);
s_scan = zeros(1, N_scan);     
v_limit_raw = zeros(1, N_scan); 
kappa_scan = zeros(1, N_scan); 

current_arc = 0;
last_q_sc = q_via_joint(:,1)'; % S-Curve 扫描用的热启动 IK

for i = 1:N_scan  % 纯几何路径参数
    CK = CurveDerivs1(u_scan(i), U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1); 
    P_p_val = CK(2,:); 
    P_pp_val = CK(3,:);
    
    CK_q = CurveDerivs1(u_scan(i), U, p, P_quat_ctrl_T, 2, num_ctrl_pts - 1);
    Q_p_val = CK_q(2,:);
    
    if i > 1
        dist_step = norm(P_p_val) * (u_scan(i) - u_scan(i-1));
        current_arc = current_arc + dist_step;
    end
    s_scan(i) = current_arc; % 将参数u转化为弧长s, u -> s
    
    % 1. 曲率/笛卡尔限速
    norm_p = norm(P_p_val); % || P' ||
    norm_cross = norm(cross(P_p_val, P_pp_val)); % || P' x P'' ||
    if norm_p > 1e-6
        kappa = norm_cross / (norm_p^3); % 曲率公式
    else
        kappa = 0;
    end
    kappa_scan(i) = kappa;
    
    if kappa > 1e-3
        v_cart_lim = sqrt(a_design_limit / kappa); % 基于向心加速度的约束
    else
        v_cart_lim = v_max_lin;
    end
    v_cart_lim = min(v_cart_lim, v_max_lin);
    
    % 2. 关节速度限制检查 
    if norm_p > 1e-6
        % 计算位姿
        P_curr_sc = CurvePoint(num_ctrl_pts-1, p, U, P_pos_ctrl_T, u_scan(i));
        Q_curr_sc = CurvePoint(num_ctrl_pts-1, p, U, P_quat_ctrl_T, u_scan(i));
        Q_curr_sc = Q_curr_sc / norm(Q_curr_sc);
        T_curr_sc = eye(4); 
        T_curr_sc(1:3,1:3) = quat2rotm(Q_curr_sc'); 
        T_curr_sc(1:3,4) = P_curr_sc;
        
        try
            q_sol_sc = ur5.ikine(T_curr_sc, 'q0', last_q_sc, 'mask', [1 1 1 1 1 1], 'tol', 1e-3, 'maxiter', 10);
            last_q_sc = q_sol_sc;
        catch
            q_sol_sc = last_q_sc;
        end
        
        J_sc = ur5.jacob0(q_sol_sc);
        
        % 计算几何速度 (ds = 1 的情况下)
        % 注意：这里的 P_p_val 是 dP/du。我们需要 dP/ds = (dP/du) / (ds/du) = P_p_val / norm_p
        % 或者简单理解：机器人以 1m/s 速度沿轨迹运动时，关节速度是多少？
        v_unit = P_p_val' / norm_p; % 单位速度向量
        
        % 计算对应的角速度单位向量
        Q_val = Q_curr_sc; Q_dot = Q_p_val';
        w_raw = 2 * (Q_dot(1)*(-Q_val(2:4)) + Q_val(1)*Q_dot(2:4) + cross(Q_dot(2:4), -Q_val(2:4)));
        % w_unit = w_raw / (ds/du) = w_raw / norm_p
        w_unit = w_raw / norm_p;
        
        twist_unit = [v_unit; w_unit];
        q_dot_unit = pinv(J_sc) * twist_unit; % 单位线速度下的关节速度
        
        max_q_dot_unit = max(abs(q_dot_unit));
        if max_q_dot_unit > 1e-6
            % limit_q_vel >= max_q_dot_unit * v_actual
            v_joint_lim = limit_q_vel / max_q_dot_unit;
        else
            v_joint_lim = 1e10;
        end
    else
        v_joint_lim = 1e10;
    end
    
    % 取二者较小值
    v_limit_raw(i) = min(v_cart_lim, v_joint_lim);
end
total_arc_length = s_scan(end);

% 平滑与关键点提取
window_size = round(N_scan * 0.05); 
v_limit_smooth = movmin(v_limit_raw, window_size);
[pks, locs] = findpeaks(-v_limit_smooth, 'MinPeakProminence', 0.1); % 寻找速度局部最小值
bottleneck_idxs = locs; 
bottleneck_vs = -pks;

key_indices = [1, bottleneck_idxs, N_scan]; % 分段点索引
key_s = s_scan(key_indices);                % 分段点所对应的弧长s
key_v_limit = [0, bottleneck_vs, 0];        % 分段点速度

% 合并短路段
min_seg_dist = 0.15; 
is_merged = true;
while is_merged
    is_merged = false; % 是否因路径太短而合并标志

    if length(key_s) <= 2
        break; 
    end
    new_indices = []; 
    new_s = []; 
    new_v = []; 
    skip_next = false; % 防止重复合并而跳过标志

    for k = 1:length(key_s)-1
        if skip_next
            skip_next = false; 
            continue; 
        end
        dist = key_s(k+1) - key_s(k);
        if k == length(key_s)-1  % 最后一段
            if dist < min_seg_dist
                is_merged = true; 
            else  % 正常拼接
                new_indices = [new_indices, key_indices(k)]; 
                new_s = [new_s, key_s(k)]; 
                new_v = [new_v, key_v_limit(k)]; 
            end
            % 去掉倒数第二个点
            new_indices = [new_indices, key_indices(end)]; 
            new_s = [new_s, key_s(end)]; 
            new_v = [new_v, key_v_limit(end)];
            continue;
        end
        if dist < min_seg_dist  % 剩余段
            v_safe = min(key_v_limit(k), key_v_limit(k+1));
            new_indices = [new_indices, key_indices(k)]; 
            new_s = [new_s, key_s(k)]; 
            new_v = [new_v, v_safe];
            skip_next = true; 
            is_merged = true;
        else
            new_indices = [new_indices, key_indices(k)]; 
            new_s = [new_s, key_s(k)]; 
            new_v = [new_v, key_v_limit(k)];
        end
    end
    if is_merged % 更新去除短段落的值
        key_indices = new_indices; 
        key_s = new_s; 
        key_v_limit = new_v; 
    end
end
num_segments = length(key_indices) - 1; % 段落个数

% 可行性前瞻
key_v_feasible = key_v_limit; 
a_planning = 0.5 * a_design_limit;
% 反向扫描
for k = length(key_s)-1 : -1 : 1
    ds = key_s(k+1) - key_s(k);
    key_v_feasible(k) = min(key_v_feasible(k), sqrt(key_v_feasible(k+1)^2 + 2 * a_planning * ds));
end
% 正向扫描
for k = 1 : length(key_s)-1
    ds = key_s(k+1) - key_s(k);
    key_v_feasible(k+1) = min(key_v_feasible(k+1), sqrt(key_v_feasible(k)^2 + 2 * a_planning * ds));
end

% 生成分段 S 曲线
t_total_segments = []; 
s_total_segments = []; 
v_total_segments = []; 
a_total_segments = [];
current_t_offset = 0; 
current_s_offset = 0;

% 对每一段进行双s速度规划
for i = 1:num_segments
    s_start = key_s(i); 
    s_end = key_s(i+1); 
    dist_seg = s_end - s_start;
    v_start = key_v_feasible(i); 
    v_end = key_v_feasible(i+1);
    idx_start = key_indices(i); 
    idx_end = key_indices(i+1);
    seg_geo_max = max(v_limit_smooth(idx_start:idx_end));
    v_cruise_final = max([min(v_max_lin, seg_geo_max), v_start, v_end]);
    
    [q, v, a, ~, t_seg, T_seg, ~] = double_s_curve_complete(0, dist_seg, v_start, v_end, v_cruise_final, a_design_limit, j_max);
    
    start_idx = 1; 
    if i > 1
        start_idx = 2; 
    end
    t_total_segments = [t_total_segments, t_seg(start_idx:end) + current_t_offset];
    s_total_segments = [s_total_segments, q(start_idx:end) + current_s_offset];
    v_total_segments = [v_total_segments, v(start_idx:end)];
    a_total_segments = [a_total_segments, a(start_idx:end)];
    current_t_offset = current_t_offset + T_seg; 
    current_s_offset = current_s_offset + dist_seg;
end

% 安全检查 
u_mapped_check = interp1(s_scan, u_scan, s_total_segments, 'linear', 'extrap');
kappa_check = interp1(s_scan, kappa_scan, s_total_segments, 'nearest', 'extrap');
a_total_check = sqrt(a_total_segments.^2 + (v_total_segments.^2 .* kappa_check).^2); % 总加速度
max_a_total = max(a_total_check);

time_scaling_ratio = 1.0;
if max_a_total > a_max_lin
    time_scaling_ratio = sqrt(max_a_total / a_max_lin) * 1.01; 
    fprintf('  [修正] 加速度超限 (%.2f)，执行缩放 %.3f\n', max_a_total, time_scaling_ratio);
end
t_total_segments = t_total_segments * time_scaling_ratio;
v_total_segments = v_total_segments / time_scaling_ratio;
a_total_segments = a_total_segments / (time_scaling_ratio^2);
T_scurve = t_total_segments(end);

% 导出数据
[s_scan_unique, unique_idx] = unique(s_scan); 
u_scan_unique = u_scan(unique_idx);
u_scurve_mapped = interp1(s_scan_unique, u_scan_unique, s_total_segments, 'pchip');
u_scurve_mapped = max(0, min(1, u_scurve_mapped)); 
t_grid_scurve = t_total_segments; 
s_vel_scurve = v_total_segments; 
s_acc_scurve = a_total_segments;
fprintf('分段 S-Curve 规划完成，总时间: %.4f 秒\n', T_scurve);
toc;

% 方法三: ISSA 
fprintf('\n===== 4. 方法三: ISSA  ====\n');
tic;

num_segments_issa = 15; 
lb_issa = 0.05 * ones(1, num_segments_issa); 
ub_issa = 1.5 * ones(1, num_segments_issa);
dim_issa = num_segments_issa;

% 预计算评估点
N_eval = 200; 
u_eval = linspace(0, 1, N_eval);
P_p_eval = zeros(3, N_eval); 
P_pp_eval = zeros(3, N_eval);
Q_p_eval = zeros(4, N_eval); 
Q_curr_eval = zeros(4, N_eval);

% 优化：预存储伪逆矩阵，避免在循环中重复计算
pinv_J_all = zeros(6, 6, N_eval); % 存储每个点的伪逆矩阵

last_q_issa = q_via_joint(:,1)';
fprintf('预计算 ISSA 评估点 (IK + Jacobian Pinv)...\n');

for i=1:N_eval
    % 1. 几何导数
    CK = CurveDerivs1(u_eval(i), U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p_eval(:,i) = CK(2,:)'; 
    P_pp_eval(:,i) = CK(3,:)';

    CK_q = CurveDerivs1(u_eval(i), U, p, P_quat_ctrl_T, 2, num_ctrl_pts - 1);
    Q_p_eval(:,i) = CK_q(2,:)';
    
    % 2. 当前位姿
    P_curr = CurvePoint(num_ctrl_pts-1, p, U, P_pos_ctrl_T, u_eval(i));
    Q_curr = CurvePoint(num_ctrl_pts-1, p, U, P_quat_ctrl_T, u_eval(i)); 
    Q_curr = Q_curr / norm(Q_curr);
    Q_curr_eval(:,i) = Q_curr;
    
    % 3. 逆解 IK
    T_curr = eye(4); 
    T_curr(1:3,1:3) = quat2rotm(Q_curr'); 
    T_curr(1:3,4) = P_curr;
    try
        q_sol = ur5.ikine(T_curr, 'q0', last_q_issa, 'mask', [1 1 1 1 1 1], 'tol', 1e-3, 'maxiter', 10);
        last_q_issa = q_sol;
    catch
        q_sol = last_q_issa; 
    end
    
    % 4. 关键优化：在这里就把 Jacobian 的逆算好存起来
    J = ur5.jacob0(q_sol);
    pinv_J_all(:, :, i) = pinv(J); 
end

% 定义代价函数 
fobj = @(dt_vec) TrajectoryCost_Optimized(dt_vec, u_eval, P_p_eval, P_pp_eval, Q_p_eval, Q_curr_eval, pinv_J_all,...
                                          v_max_lin, a_max_lin, limit_q_vel);

% 运行 ISSA
popsize = 30; 
max_iter = 1000; 
fprintf('运行 ISSA (种群: %d, 迭代: %d)...\n', popsize, max_iter);
[best_fitness, best_dt, ~] = ISSA(popsize, max_iter, lb_issa, ub_issa, dim_issa, fobj);

T_issa = sum(best_dt);
fprintf('ISSA 优化完成，总时间: %.4f 秒\n', T_issa);
toc;

% 重构数据用于绘图
time_nodes_opt = [0, cumsum(best_dt)];
u_nodes_opt = linspace(0, 1, num_segments_issa + 1);
dt_plot = 0.01; 
t_plot = 0:dt_plot:max([T_toppra, T_scurve, T_issa]);
vel_mag_issa = zeros(size(t_plot)); 
acc_mag_issa = zeros(size(t_plot));
for k = 1:length(t_plot)
    t_now = t_plot(k);
    if t_now > T_issa
        continue; 
    end
    u_now = interp1(time_nodes_opt, u_nodes_opt, t_now, 'pchip');
    
    delta_t = 1e-4;
    t_minus = max(0, t_now - delta_t); 
    t_plus = min(T_issa, t_now + delta_t);
    u_minus = interp1(time_nodes_opt, u_nodes_opt, t_minus, 'pchip');
    u_plus  = interp1(time_nodes_opt, u_nodes_opt, t_plus, 'pchip');
    u_dot = (u_plus - u_minus) / (t_plus - t_minus);
    
    u_dot_minus = (interp1(time_nodes_opt, u_nodes_opt, t_now, 'pchip') -...
                   interp1(time_nodes_opt, u_nodes_opt, t_now-delta_t, 'pchip'))/delta_t;
    u_dot_plus  = (interp1(time_nodes_opt, u_nodes_opt, t_now+delta_t, 'pchip') -...
                   interp1(time_nodes_opt, u_nodes_opt, t_now, 'pchip'))/delta_t;
    u_ddot = (u_dot_plus - u_dot_minus) / delta_t;
    
    CK = CurveDerivs1(u_now, U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p_vec = CK(2,:)'; 
    P_pp_vec = CK(3,:)';
    vel_mag_issa(k) = norm(P_p_vec * u_dot);
    acc_mag_issa(k) = norm(P_p_vec * u_ddot + P_pp_vec * u_dot^2);
end

% 5. 综合对比绘图
fprintf('===== 5. 综合对比绘图 =====');
% 为了简洁，这里重用上面的 vel_mag_issa 等变量
% 需要补全 TOPP-RA 和 S-Curve 的 t_plot 对应数据
vel_mag_toppra = zeros(size(t_plot)); 
acc_mag_toppra = zeros(size(t_plot));
for k = 1:length(t_plot)
    t_now = min(t_plot(k), T_toppra);
    u_now = interp1(t_grid_toppra, s_grid, t_now, 'pchip');
    s_dot_now = interp1(t_grid_toppra, sd_optimal, t_now, 'pchip');
    idx = find(t_grid_toppra <= t_now, 1, 'last');
    if isempty(idx) || idx >= length(t_grid_toppra)
        idx = length(t_grid_toppra)-1;
    end
    dt_loc = t_grid_toppra(idx+1) - t_grid_toppra(idx);
    s_ddot_now = (sd_optimal(idx+1) - sd_optimal(idx)) / dt_loc;
    CK = CurveDerivs1(u_now, U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);

    vel_mag_toppra(k) = norm(CK(2,:)' * s_dot_now);
    acc_mag_toppra(k) = norm(CK(2,:)' * s_ddot_now + CK(3,:)' * s_dot_now^2);
end

vel_mag_scurve = zeros(size(t_plot)); 
acc_mag_scurve = zeros(size(t_plot));
for k = 1:length(t_plot)
    t_now = min(t_plot(k), T_scurve);
    
    s_dot_now = interp1(t_grid_scurve, s_vel_scurve, t_now, 'linear');
    s_ddot_now = interp1(t_grid_scurve, s_acc_scurve, t_now, 'linear'); 
    
    u_now = interp1(t_grid_scurve, u_scurve_mapped, t_now, 'linear');
    CK = CurveDerivs1(u_now, U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_prime = CK(2,:)'; 
    P_double_prime = CK(3,:)';
    norm_Pp = norm(P_prime);
    
    if norm_Pp > 1e-6
        cross_val = cross(P_prime, P_double_prime);
        kappa_now = norm(cross_val) / (norm_Pp^3);
        a_normal_now = kappa_now * s_dot_now^2;
        acc_mag_scurve(k) = sqrt(s_ddot_now^2 + a_normal_now^2);
    else
        acc_mag_scurve(k) = 0;
    end
    vel_mag_scurve(k) = s_dot_now;
end

% 计算加速度分量 (切向/法向)
t = t_plot(:); 
dt_g = gradient(t);

v_top = vel_mag_toppra(:); 
a_top_tot = acc_mag_toppra(:);
at_top = gradient(v_top) ./ (dt_g + 1e-9); 
an_top = sqrt(max(0, a_top_tot.^2 - at_top.^2));

v_sc = vel_mag_scurve(:); 
a_sc_tot = acc_mag_scurve(:);
at_sc = gradient(v_sc) ./ (dt_g + 1e-9); 
an_sc = sqrt(max(0, a_sc_tot.^2 - at_sc.^2));

v_issa = vel_mag_issa(:); 
a_issa_tot = acc_mag_issa(:);
at_issa = gradient(v_issa) ./ (dt_g + 1e-9); 
an_issa = sqrt(max(0, a_issa_tot.^2 - at_issa.^2));

% 绘图
figure('Position', [50, 50, 1200, 1000], 'Color', 'w');
sgtitle('轨迹规划算法对比', 'FontSize', 16);

subplot(3, 2, 1);
plot(t, v_top, 'b-', 'LineWidth', 1.5); 
hold on;
plot(t, v_sc, 'r--', 'LineWidth', 1.5);
plot(t, v_issa, 'g-.', 'LineWidth', 2);
yline(v_max_lin, 'k:', 'Limit');
ylabel('Speed (m/s)'); 
title('1. 线速度 v(t)'); 
legend('TOPP-RA','S-Curve','ISSA'); 
grid on;

subplot(3, 2, 2);
b = barh([1,2,3], [T_toppra, T_scurve, T_issa]);
b.FaceColor = 'flat'; 
b.CData(1,:) = [0 0 1]; 
b.CData(2,:) = [1 0 0]; 
b.CData(3,:) = [0 1 0];
yticklabels({'TOPP-RA', 'S-Curve', 'ISSA'}); 
title('2. 总耗时对比');
text(T_toppra/2, 1, sprintf('%.3fs',T_toppra),'Color','w');
text(T_scurve/2, 2, sprintf('%.3fs',T_scurve),'Color','w');
text(T_issa/2, 3, sprintf('%.3fs',T_issa),'Color','k'); 

subplot(3, 2, 3);
plot(t, at_top, 'b-', 'LineWidth', 1.5); 
hold on;
plot(t, at_sc, 'r--', 'LineWidth', 1.5);
plot(t, at_issa, 'g-.', 'LineWidth', 1.5);
yline(a_max_lin, 'k:', 'Limit'); 
yline(-a_max_lin, 'k:', 'Limit');
ylabel('Acc_t (m/s^2)'); 
title('3. 切向加速度'); 
grid on; 
ylim([-a_max_lin*1.5, a_max_lin*1.5]);

subplot(3, 2, 4);
plot(t, an_top, 'b-', 'LineWidth', 1.5); 
hold on;
plot(t, an_sc, 'r--', 'LineWidth', 1.5);
plot(t, an_issa, 'g-.', 'LineWidth', 1.5);
yline(a_max_lin, 'k:', 'Limit');
ylabel('Acc_n (m/s^2)'); 
title('4. 法向加速度'); 
grid on; 
ylim([0, a_max_lin*1.5]);

subplot(3, 2, [5,6]);
plot(t, a_top_tot, 'b-', 'LineWidth', 1.5); 
hold on;
plot(t, a_sc_tot, 'r--', 'LineWidth', 1.5);
plot(t, a_issa_tot, 'g-.', 'LineWidth', 1.5);
yline(a_max_lin, 'k:', 'Limit', 'LineWidth', 2);
ylabel('Acc_total (m/s^2)'); 
title('5. 总加速度'); 
legend('TOPP-RA', 'S-Curve', 'ISSA', 'Location', 'Best'); 
grid on; 
ylim([0, a_max_lin*1.5]);

fprintf('绘图完成。');

% 代价函数
function cost = TrajectoryCost_Optimized(dt_vec, u_eval, P_p_eval, P_pp_eval, Q_p_eval, Q_curr_eval, pinv_J_all,...
                                         v_max, a_max, q_vel_limit)
    num_seg = length(dt_vec);
    u_nodes = linspace(0, 1, num_seg + 1);
    t_nodes = [0, cumsum(dt_vec)];
    T_total = t_nodes(end);
    
    % 快速插值
    t_eval = interp1(u_nodes, t_nodes, u_eval, 'pchip');
    dt_grad = gradient(t_eval); 
    du_grad = gradient(u_eval); 
    
    % 计算 u 对 t 的导数
    t_prime_u = dt_grad ./ (du_grad + 1e-12); 
    u_dot = 1.0 ./ (t_prime_u + 1e-12);       
    u_dot_grad_u = gradient(u_dot) ./ (du_grad + 1e-12); 
    u_ddot = u_dot_grad_u .* u_dot;
    
    % 1. 笛卡尔违规 (矩阵化计算，极快)
    v_vec_sq = sum(P_p_eval.^2, 1) .* (u_dot.^2);
    v_norm = sqrt(v_vec_sq);
    a_vec = P_p_eval .* u_ddot + P_pp_eval .* (u_dot.^2);
    a_norm = sqrt(sum(a_vec.^2, 1));
    
    v_violation = max(0, v_norm - v_max);
    a_violation = max(0, a_norm - a_max);
    
    % 2. 关节速度违规 (查表法)
    % 这里无法完全矩阵化，因为 pinv_J 是三维数组，但比调用 robot 对象快得多
    q_vel_violation = zeros(size(u_eval));
    
    % 预计算角速度系数 w_geom_base
    % Q_dot_geom = Q_p_eval; 
    % w = 2 * dQ/du * Q_inv
    % 这一步其实也可以提到主循环外预计算，但这里算也很快
    for i = 1:length(u_eval)
        Q_val = Q_curr_eval(:, i);
        Q_dot_geom = Q_p_eval(:, i);
        
        % 几何角速度向量 (相对于 u)
        w_geom_base = 2 * (Q_dot_geom(1)*(-Q_val(2:4)) + Q_val(1)*Q_dot_geom(2:4) + cross(Q_dot_geom(2:4), -Q_val(2:4)));
        
        % 构建 Twist 向量 [v; w]
        twist = [P_p_eval(:, i); w_geom_base] * u_dot(i);
        
        % 直接矩阵乘法
        q_dot = pinv_J_all(:, :, i) * twist;
        
        max_q_dot = max(abs(q_dot));
        if max_q_dot > q_vel_limit
            q_vel_violation(i) = max_q_dot - q_vel_limit;
        end
    end
    
    penalty_v = sum(v_violation.^2);
    penalty_a = sum(a_violation.^2);
    penalty_q = sum(q_vel_violation.^2);
    
    w1 = 1.0; 
    w2 = 1e5; 
    w3 = 1e5;
    cost = w1 * T_total + w2 * (penalty_v + penalty_a) + w3 * penalty_q;
    
    if isnan(cost)
        cost = 1e10; 
    end
end

% ISSA 优化算法 
function [best_fitness, best_position, convergence_curve] = ISSA(popsize, max_iter, lb, ub, dim, fobj)
    p_percent = 0.2; 
    p_num = round(popsize * p_percent);
    sd_percent = 0.1; 
    sd_num = round(popsize * sd_percent);
    ST = 0.8;

    Y = zeros(popsize, dim); 
    k = 4;
    y_seed = 2 * rand(1, dim) - 1; 
    Y(1, :) = y_seed;
    for i = 2:popsize
        val = cos(k * acos(Y(i-1, :))); 
        val(val > 1) = 1; 
        val(val < -1) = -1; 
        Y(i, :) = val; 
    end

    positions = zeros(popsize, dim);
    for d = 1:dim
        positions(:, d) = lb(d) + (1 + Y(:, d)) .* (ub(d) - lb(d)) / 2; 
    end

    fitness = zeros(popsize, 1);
    for i = 1:popsize
        fitness(i) = fobj(positions(i, :));
    end
    [best_fitness, best_idx] = min(fitness);
    best_position = positions(best_idx, :);
    convergence_curve = zeros(max_iter, 1);

    for t = 1:max_iter
        [~, sort_idx] = sort(fitness);
        best_pos_t = positions(sort_idx(1), :); worst_pos_t = positions(sort_idx(end), :); 
        for i = 1:p_num
            idx = sort_idx(i); 
            if rand() < ST
                positions(idx, :) = positions(idx, :) .* exp(-i / (rand() * max_iter));
            else
                positions(idx, :) = positions(idx, :) + randn() * ones(1, dim); 
            end
        end

        for i = (p_num + 1):popsize
            idx = sort_idx(i); 
            if i > (popsize + p_num) / 2
                positions(idx, :) = randn() * exp((worst_pos_t - positions(idx, :)) / (i^2));
            else
                A = ones(1, dim);
                A_plus = 1 ./ (A' * (A * A')^(-1)); 
                positions(idx, :) = best_pos_t + abs(positions(idx, :) - best_pos_t) .* A_plus' .* randn(); 
            end
        end

        rand_indices = randperm(popsize, sd_num); 
        for i = 1:sd_num
            idx = rand_indices(i); 
            if fitness(idx) > best_fitness
                positions(idx, :) = best_position + randn() * abs(positions(idx, :) - best_position);
            else
                denominator = fitness(idx) - fitness(sort_idx(end)) + 1e-8; 
                K = 2 * rand() - 1; 
                positions(idx, :) = positions(idx, :) + K * (abs(positions(idx, :) - worst_pos_t) ./ denominator);
            end
        end

        if rand < 0.2
            start_index = floor(popsize / 2) + 1; 
            for i = start_index : popsize
                t_rand_vec = trnd(t, [1, dim]); 
                positions(i, :) = positions(i, :) + positions(i, :) .* t_rand_vec; 
            end
        else
            alpha = 0.5; 
            candidate_position = best_position + alpha * (rand(1, dim) - 0.5);
            candidate_position = max(candidate_position, lb); 
            candidate_position = min(candidate_position, ub);
            candidate_fitness = fobj(candidate_position);
            if candidate_fitness < best_fitness
                best_fitness = candidate_fitness;
                best_position = candidate_position; 
            end
        end

        for i = 1:popsize
            positions(i, :) = max(positions(i, :), lb);
            positions(i, :) = min(positions(i, :), ub); 
            fitness(i) = fobj(positions(i, :));
        end

        [current_best_fitness, current_best_idx] = min(fitness);
        if current_best_fitness < best_fitness
            best_fitness = current_best_fitness; 
            best_position = positions(current_best_idx, :); 
        end
        convergence_curve(t) = best_fitness;

        if mod(t, 20) == 0 || t == 1
            fprintf('  \nISSA Iter %d: \n Best Cost = %.4f', t, best_fitness); 
        end
    end
end

% 双S曲线函数 
function [q_all, q_d_all, q_dd_all, q_ddd_all, t_segment, T, params] = double_s_curve_complete(q0, q1, v0, v1, vmax, amax, jmax)
    sigma = sign(q1 - q0); 
    if sigma == 0
        sigma = 1; 
    end

    dist = abs(q1 - q0); 
    v0 = abs(v0); 
    v1 = abs(v1); 
    vmax = abs(vmax); 
    amax = abs(amax); 
    jmax = abs(jmax);
    isValid = false; 
    v_target = vmax;

    for iter = 1:50 
        [Ta, Tv, Td, Tj1, Tj2, vlim, alim_a, alim_d] = try_calculate_times(dist, v0, v1, v_target, amax, jmax);
        if Tv >= -1e-6
            Tv = max(0, Tv); 
            isValid = true; 
            break; 
        else
            v_target = v_target * 0.9; 
        end
    end

    if ~isValid
        Ta=0.1; 
        Tv=0; 
        Td=0.1; 
        Tj1=0; Tj2=0; 
        vlim=v0; 
        alim_a=0; 
        alim_d=0;
    end

    T = Ta + Tv + Td;
    dt = 0.01; 
    t_segment = 0:dt:T; 
    if isempty(t_segment)
        t_segment=[0,T]; 
    end

    n = length(t_segment);
    q = zeros(1, n);
    v = zeros(1, n);
    a = zeros(1, n); 
    j = zeros(1, n); 
    jmin = -jmax;
    for i = 1:n
        t = t_segment(i);
        if t <= Tj1
            q(i) = v0*t + jmax*t^3/6;
            v(i) = v0 + jmax*t^2/2;
            a(i) = jmax*t;
        elseif t <= Ta - Tj1
            t1 = t;
            q(i) = v0*t1 + alim_a/6*(3*t1^2 - 3*Tj1*t1 + Tj1^2);
            v(i) = v0 + alim_a*(t - Tj1/2);
            a(i) = alim_a;
        elseif t <= Ta
            t2 = Ta - t; 
            q(i) = (v0 + vlim)*Ta/2 - vlim*t2 - jmin*t2^3/6;
            v(i) = vlim + jmin*t2^2/2;
            a(i) = -jmin*t2;
        elseif t <= Ta + Tv
            t3 = t - Ta;
            q(i) = (v0 + vlim)*Ta/2 + vlim*t3; 
            v(i) = vlim;
        elseif t <= Ta + Tv + Tj2
            t4 = t - Ta - Tv; 
            q(i) = (v0 + vlim)*Ta/2 + vlim*Tv + vlim*t4 - jmax*t4^3/6;
            v(i) = vlim - jmax*t4^2/2; 
            a(i) = -jmax*t4;
        elseif t <= T - Tj2
            t5 = t - Ta - Tv;
            delta_s = (v0+vlim)*Ta/2 + vlim*Tv; 
            q(i) = delta_s + vlim*t5 + alim_d/6*(3*t5^2 - 3*Tj2*t5 + Tj2^2);
            v(i) = vlim + alim_d*(t5 - Tj2/2); 
            a(i) = alim_d;
        else
            t6 = T - t;
            q(i) = dist - v1*t6 - jmax*t6^3/6; 
            v(i) = v1 + jmax*t6^2/2; 
            a(i) = -jmax*t6;
        end
    end

    if sigma == 1
        q_all = q0 + q;
        q_d_all = v;
        q_dd_all = a;
        q_ddd_all = j;
    else
        q_all = q0-q;
        q_d_all = -v;
        q_dd_all = -a;
        q_ddd_all = -j;
    end
    params.T = T;
end

function [Ta, Tv, Td, Tj1, Tj2, vlim, alim_a, alim_d] = try_calculate_times(dist, v0, v1, v_target, amax, jmax)
    if (v_target - v0) * jmax < amax^2
        Tj1 = sqrt(abs(v_target - v0) / jmax); 
        Ta = 2 * Tj1;
        alim_a = jmax * Tj1 * sign(v_target - v0);
    else
        Tj1 = amax / jmax; 
        Ta = Tj1 + abs(v_target - v0) / amax; 
        alim_a = amax * sign(v_target - v0); 
    end

    if abs(v_target - v1) * jmax < amax^2
        Tj2 = sqrt(abs(v_target - v1) / jmax);
        Td = 2 * Tj2; 
        alim_d = -jmax * Tj2 * sign(v_target - v1);
    else
        Tj2 = amax / jmax; 
        Td = Tj2 + abs(v_target - v1) / amax;
        alim_d = -amax * sign(v_target - v1); 
    end
    
    S_acc = (v0 + v_target) * Ta / 2; 
    S_dec = (v_target + v1) * Td / 2; 
    S_remaining = dist - S_acc - S_dec; 
    Tv = S_remaining / v_target; 
    vlim = v_target;
end

% 求解满足约束的最大最小加速度
function [alpha_min, alpha_max] = GetSDDotBounds_Exact(Pp, Ppp, Qp, Qpp, x, a_max, wa_max)
    % a = dv/dt = P' × s̈ + P'' × ṡ² = P' × s̈ + P'' × x ≤ a_max 两边平方
    alpha_min = -1e10; 
    alpha_max = 1e10;
    [min1, max1] = SolveQuad(dot(Pp, Pp), 2*dot(Pp, Ppp*x), dot(Ppp*x, Ppp*x) - a_max^2);
    alpha_min = max(alpha_min, min1); 
    alpha_max = min(alpha_max, max1);

    % α = dω/dt ≈ 2 × (q' × s̈ + q'' × x) ≤ α_max 两边平方
    A_w = 2*Qp; 
    B_w = 2*Qpp*x;
    [min2, max2] = SolveQuad(dot(A_w, A_w), 2*dot(A_w, B_w), dot(B_w, B_w) - wa_max^2);
    alpha_min = max(alpha_min, min2); 
    alpha_max = min(alpha_max, max2);

    if alpha_min > alpha_max
        alpha_max = alpha_min;
    end
end

% 求解一元二次方程
function [xmin, xmax] = SolveQuad(a, b, c)
    if abs(a) < 1e-10      % a为0
        if abs(b) < 1e-10  % b为0
            if c <= 1e-10  % c为0
                xmin = -1e10;
                xmax = 1e10; 
            else
                xmin = 1e10;
                xmax = -1e10;
            end
        elseif b > 0
            xmin = -1e10; 
            xmax = -c/b;
        else
            xmin = -c/b;
            xmax = 1e10; 
        end
    else  % 正常情况
        delta = b^2 - 4*a*c; 
        if delta < 0 && delta > -1e-8
            delta = 0;
        end

        if delta < 0  % 无解
            if a > 0
                xmin = 1e10; 
                xmax = -1e10; 
            else
                xmin = -1e10; 
                xmax = 1e10; 
            end
        else
            d = sqrt(delta); 
            r1 = (-b - d)/(2*a); 
            r2 = (-b + d)/(2*a);
            if a > 0
                xmin = min(r1,r2); 
                xmax = max(r1,r2); 
            else
                xmin = -1e10; 
                xmax = 1e10; 
            end
        end
    end
end

function U = BuildKnotVector(p, n, u_bar, m)
    U = zeros(1, n + p + 1); 
    U(end-p:end) = 1; 
    step = (m - 1) / (n - p);
    for i = 1:(n - p - 1)
        idx = 1 + i * step;
        w = idx - floor(idx);
        U(p+1+i) = (1 - w) * u_bar(floor(idx)) + w * u_bar(ceil(idx));
    end
end

function [B, W, Q_target] = BuildConstraints_Cartesian(p, n, u_bar, q_via, v0, v6, a0, a6, U)
    dim = size(q_via, 1); 
    m = size(q_via, 2); 
    num_t = m + 4;
    Q_target = [q_via'; v0'; a0'; v6'; a6']; 
    W = diag([1e6, ones(1, m-2), 1e6, 1e6*ones(1,4)]);
    B = zeros(num_t, n);

    for k = 1:m
        s = FindSpan(n-1, p, u_bar(k), U); 
        N = BasisFuns(s, u_bar(k), p, U); 
        B(k, s-p+1:s+1) = N'; 
    end

    s0 = FindSpan(n-1, p, u_bar(1), U);
    d0 = DersBasisFuns(s0, u_bar(1), p, 2, U); 
    B(m+1, s0-p+1:s0+1) = d0(2, :); 
    B(m+2, s0-p+1:s0+1) = d0(3, :);

    s6 = FindSpan(n-1, p, u_bar(end), U); 
    d6 = DersBasisFuns(s6, u_bar(end), p, 2, U);
    B(m+3, s6-p+1:s6+1) = d6(2, :); 
    B(m+4, s6-p+1:s6+1) = d6(3, :);
end

function [C, A] = BuildSmoothMatrix(n, p, U)
    C = zeros(n-2, n);
    A = zeros(n-2);
    u_t = U(n+1) - U(p+1);

    for k = 1:(n-2)
        d1 = U(k+p+1)-U(k+1);
        d2 = U(k+p+2)-U(k+2);
        if d1>1e-10
            C(k,k) = 6/d1;
            C(k,k+1) = -6/d1;
        end
        if d2>1e-10
            C(k,k+1) = C(k,k+1) - 6/d2;
            C(k,k+2) = 6/d2; 
        end
        
        A(k,k) = 2*(U(k+p+2) - U(k+p+1))*u_t;
        if k<n-2
            A(k,k+1) = (U(k+p+2)-U(k+p+1))*u_t;
            A(k+1,k) = A(k,k+1); 
        end
    end
end

function span = FindSpan(n, p, u, U)
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
    N = zeros(p+1, 1); 
    N(1) = 1.0; 
    L = zeros(p+1, 1); 
    R = zeros(p+1, 1);

    for j = 1:p
        L(j+1) = u - U(i+1-j+1);
        R(j+1) = U(i+j+1) - u;
        saved = 0.0;

        for r = 0:(j-1)
            t = N(r+1) / (R(r+2) + L(j-r+1));
            N(r+1) = saved + R(r+2) * t; 
            saved = L(j-r+1) * t;
        end
        N(j+1) = saved;
    end
end

function ders = DersBasisFuns(i, u, p, n, U)
    ders = zeros(n+1, p+1);
    ndu = zeros(p+1, p+1);
    ndu(1,1) = 1.0;
    L = zeros(p+1, 1);
    R = zeros(p+1, 1);
    for j = 1:p
        L(j+1) = u - U(i+1-j+1);
        R(j+1) = U(i+j+1) - u; 
        saved = 0.0;
        for r = 0:(j-1)
            ndu(j+1, r+1) = R(r+2) + L(j-r+1); 
            t = ndu(r+1, j) / ndu(j+1, r+1); 
            ndu(r+1, j+1) = saved + R(r+2) * t; 
            saved = L(j-r+1) * t;
        end
        ndu(j+1, j+1) = saved;
    end

    for j=0:p
        ders(1, j+1) = ndu(j+1, p+1); 
    end

    a = zeros(2, p+1);
    for r = 0:p
        s1=0; 
        s2=1; 
        a(1,1) = 1.0;
        for k = 1:n
            d = 0.0; 
            rk = r-k; 
            pk = p-k;
            if r >= k
                a(s2+1,1)=a(s1+1,1)/ndu(pk+2,rk+1);
                d=a(s2+1,1)*ndu(rk+1,pk+1); 
            end

            j1 = max(1, -rk); 
            j2 = min(k-1, p-r);
            for j = j1:j2
                a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1); 
                d = d + a(s2 + 1,j + 1)*ndu(rk+j+1,pk+1);
            end
            if r <= pk
                a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
                d = d+a(s2+1,k+1)*ndu(r+1,pk+1); 
            end
            ders(k+1,r+1)=d; 
            j=s1; 
            s1=s2; 
            s2=j;
        end
    end
    r=p; 
    for k=1:n
        for j=0:p
            ders(k+1,j+1)=ders(k+1,j+1)*r;
        end
        r=r*(p-k); 
    end
end

function CK = CurveDerivs1(u, U, p, P, d, n)
    du = min(d, p);
    CK = zeros(d+1, size(P, 1)); 
    s = FindSpan(n, p, u, U); 
    nders = DersBasisFuns(s, u, p, du, U);
    for k = 0:du
        CK(k+1, :) = nders(k+1, :) * P(:, s-p+1:s+1)';
    end
end

function C = CurvePoint(n, p, U, P, u)
    span = FindSpan(n, p, u, U);
    N = BasisFuns(span, u, p, U); 
    C = zeros(size(P, 1), 1);

    for i = 0:p
        C = C + N(i+1) * P(:, span-p+i+1); 
    end
end
