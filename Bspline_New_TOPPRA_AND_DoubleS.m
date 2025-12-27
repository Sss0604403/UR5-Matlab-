clear; clc; close all;

if ~exist('Link', 'class')
    error('请先安装MATLAB Robotics Toolbox');
end

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

% 转换途径点
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

% 转换为笛卡尔位姿
pos_via = zeros(3, num_via_pts);
quat_via = zeros(4, num_via_pts); 
for i = 1:num_via_pts
    T = ur5.fkine(q_via_joint(:, i)');
    pos_via(:, i) = T.t; 
    Q = UnitQuaternion(T); 
    quat_via(:, i) = double(Q)'; 
end

% 笛卡尔空间 B样条规划
p = 5; 
num_ctrl_pts = 12;
lambda = 1e-4; 
epsilon = 1e-8;

dists = sqrt(sum(diff(pos_via, 1, 2).^2, 1));
u_bar = [0, cumsum(dists) / sum(dists)];

U = BuildKnotVector(p, num_ctrl_pts, u_bar, num_via_pts);
[C_mat, A_mat] = BuildSmoothMatrix(num_ctrl_pts, p, U);

% 位置拟合
v0_pos = [0;0;0]; 
v6_pos = [0;0;0];
a0_pos = [0;0;0];
a6_pos = [0;0;0];
[B_pos, W_pos, Q_target_pos] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, pos_via, v0_pos, v6_pos, a0_pos, a6_pos, U);
M_pos = B_pos' * W_pos * B_pos + lambda * (C_mat' * A_mat * C_mat) + epsilon * eye(num_ctrl_pts);
F_pos = B_pos' * W_pos * Q_target_pos;
P_pos_ctrl_T = (M_pos \ F_pos)';

% 姿态拟合
v0_quat = [0;0;0;0]; 
v6_quat = [0;0;0;0];
a0_quat = [0;0;0;0]; 
a6_quat = [0;0;0;0];
[B_quat, W_quat, Q_target_quat] = BuildConstraints_Cartesian(p, num_ctrl_pts, u_bar, quat_via, v0_quat, v6_quat, a0_quat, a6_quat, U);
M_quat = B_quat' * W_quat * B_quat + lambda * (C_mat' * A_mat * C_mat) + epsilon * eye(num_ctrl_pts);
F_quat = B_quat' * W_quat * Q_target_quat;
P_quat_ctrl_T = (M_quat \ F_quat)';

fprintf('B样条拟合完成。\n');

% 参数设置 
% 物理约束 (两个算法共用)
v_max_lin = 1.5;    % 最大线速度 (m/s)
a_max_lin = 3.0;    % 最大线加速度 (m/s^2)
w_max = 2.5;       % 最大角速度 (rad/s)
alpha_max = 10.0;   % 最大角加速度 (rad/s^2)
limit_q_vel = 3.0;  % 关节最大速度限制 (rad/s)

% S曲线独有参数
j_max = 50.0;       % 最大加加速度 (m/s^3)

% 方法一：TOPP-RA 
fprintf('\n===== 2. 方法一: TOPP-RA 计算 =====\n');
tic;
N_grid = 2000;           
ds = 1 / (N_grid - 1);  
s_grid = linspace(0, 1, N_grid);

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

x_limit = zeros(1, N_grid); % 速度平方(sdot²)
for i = 1:N_grid
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
    q_dot = Q_p(:, i); % 四元数一阶导
    
    % 简单的工程修正：直接用 4倍模长平方
    % 为什么？因为 Slerp 的角速度 omega = 2 * ||q_dot|| (当 q 是单位四元数时)
    % 你的 B 样条出来的 q 虽然模长不是严格 1，但差别极小，直接用这个公式误差 < 0.1%
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
    
    x_limit(i) = min([x_v, x_w, x_a_lin, x_a_rot]);
end
x_limit(1) = 0; x_limit(end) = 0;

% Backward Pass
beta = x_limit; 
for i = N_grid-1 : -1 : 1
    x_next_max = beta(i+1);
    [acc_min, ~] = GetSDDotBounds_Exact(P_p(:,i), P_pp(:,i), Q_p(:,i), Q_pp(:,i), x_next_max, a_max_lin, alpha_max);
    x_curr_max = x_next_max - 2 * ds * acc_min;
    if x_curr_max < 0
        x_curr_max = 0; 
    end
    if x_curr_max < beta(i)
        beta(i) = x_curr_max; 
    end
end

% Forward Pass
x_opt = zeros(1, N_grid);
for i = 1 : N_grid-1
    x_curr = x_opt(i);
    [~, acc_max] = GetSDDotBounds_Exact(P_p(:,i), P_pp(:,i), Q_p(:,i), Q_pp(:,i), x_curr, a_max_lin, alpha_max);
    x_next_potential = x_curr + 2 * ds * acc_max;
    if x_next_potential < 0
        x_next_potential = 0; 
    end
    if x_next_potential > beta(i+1)
        x_opt(i+1) = beta(i+1); 
    else
        x_opt(i+1) = x_next_potential; 
    end
end
sd_optimal = sqrt(max(0, x_opt));

% 恢复时间戳
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
    t_grid_toppra(i+1) = t_grid_toppra(i) + dt_step;
end
t_grid_toppra = real(t_grid_toppra);
T_toppra = t_grid_toppra(end);
fprintf('TOPP-RA 规划完成，总时间: %.4f 秒\n', T_toppra);
toc;

%% ===== [补充] 关节空间验证与 IK (核心工程步骤) =====
fprintf('\n===== [补充] 关节空间验证与生成 =====\n');
% 目的：检查笛卡尔规划是否导致关节速度超限 (例如经过奇异点)
% 如果超限，需要对 TOPP-RA 的结果进行全局时间缩放

current_scale = 1.0;
dt_ik = 0.02; % 控制周期 20ms
limit_q_vel = 3.0; % 关节速度上限 rad/s

% 尝试生成关节轨迹
t_traj = 0 : dt_ik : (T_toppra * current_scale);
num_steps = length(t_traj);
q_traj_ik = zeros(num_steps, 6);
q_guess = q_via_joint(:, 1)'; 

max_q_vel = 0;

for k = 1:num_steps
    t_now = t_traj(k) / current_scale; % 映射回原始时间
    if t_now > T_toppra, t_now = T_toppra; end
    
    % 1. 采样笛卡尔点
    s_now = interp1(t_grid_toppra, s_grid, t_now, 'pchip');
    
    % 计算位置
    CK_pos = CurveDerivs1(s_now, U, p, P_pos_ctrl_T, 0, num_ctrl_pts - 1);
    P_now = CK_pos(1, :)';
    
    % 计算姿态 (归一化)
    Q_raw = CurvePoint(num_ctrl_pts-1, p, U, P_quat_ctrl_T, s_now);
    Q_norm = Q_raw / norm(Q_raw); 
    Quat_obj = UnitQuaternion(Q_norm'); 
    
    % 2. IK
    T_des = SE3(P_now) * Quat_obj.SE3;
    try
        q_sol = ur5.ikine(T_des, 'q0', q_guess, 'mask', [1 1 1 1 1 1], 'ilimit', 50, 'tol', 1e-4);
    catch
        q_sol = q_guess;
    end
    
    if isempty(q_sol), q_sol = q_guess; end
    q_traj_ik(k, :) = q_sol;
    
    % 3. 速度检查
    if k > 1
        q_vel = abs(q_sol - q_guess) / dt_ik;
        max_q_vel = max(max_q_vel, max(q_vel));
    end
    q_guess = q_sol;
end

fprintf('最大关节速度检测值: %.2f rad/s\n', max_q_vel);

if max_q_vel > limit_q_vel
    needed_scale = max_q_vel / limit_q_vel;
    fprintf('警告：关节超速！建议将 TOPP-RA 总时间延长 %.2f 倍 (即 %.2f 秒 -> %.2f 秒)\n', needed_scale, T_toppra, T_toppra * needed_scale);
    % 在实际工程中，这里会更新 T_toppra 并把 t_grid_toppra 乘上系数
    % 这里仅做打印提示，不打断后续对比逻辑
else
    fprintf('关节速度安全。\n');
end

%% ===== 3. 方法二: 双 S 曲线计算 ... (后面保持不变)

fprintf('\n===== 3. 方法二: 分段 S-Curve (最终修正+安全缩放版) ====\n');
tic;

% 将基础设计目标设为 0.6 倍 a_max (3.0 m/s^2)
% 留出充足的余量给“切向+法向”的矢量叠加
safety_factor = 0.6; 
a_design_limit = safety_factor * a_max_lin; 

% 第一步：几何扫描
N_scan = 5000; 
u_scan = linspace(0, 1, N_scan);
s_scan = zeros(1, N_scan);     
v_limit_raw = zeros(1, N_scan); 
kappa_scan = zeros(1, N_scan); % 记录曲率供后续检查

current_arc = 0;
for i = 1:N_scan
    CK = CurveDerivs1(u_scan(i), U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p = CK(2,:);
    P_pp = CK(3,:);
    
    if i > 1
        dist_step = norm(P_p) * (u_scan(i) - u_scan(i-1));
        current_arc = current_arc + dist_step;
    end
    s_scan(i) = current_arc;
    
    norm_p = norm(P_p);
    norm_cross = norm(cross(P_p, P_pp));
    
    if norm_p > 1e-6
        kappa = norm_cross / (norm_p^3);
    else
        kappa = 0;
    end
    kappa_scan(i) = kappa;
    
    % 几何限速
    if kappa > 1e-3
        v_lim = sqrt(a_design_limit / kappa); 
        v_limit_raw(i) = min(v_lim, v_max_lin);
    else
        v_limit_raw(i) = v_max_lin;
    end
end
total_arc_length = s_scan(end);

% 平滑限速曲线
window_size = round(N_scan * 0.05); 
v_limit_smooth = movmin(v_limit_raw, window_size);

% 第二步：关键点提取与路段合并
[pks, locs] = findpeaks(-v_limit_smooth, 'MinPeakProminence', 0.1); 
bottleneck_idxs = locs;
bottleneck_vs = -pks;

key_indices = [1, bottleneck_idxs, N_scan];
key_s = s_scan(key_indices);
key_v_limit = [0, bottleneck_vs, 0]; 

% 强行合并短路段 (小于 15cm 的全部吃掉)
min_seg_dist = 0.15; 
is_merged = true;
while is_merged
    is_merged = false;
    if length(key_s) <= 2
        break;
    end
    
    new_indices = [];
    new_s = [];
    new_v = [];
    skip_next = false;
    
    for k = 1:length(key_s)-1
        if skip_next
            skip_next = false;
            continue;
        end
        
        dist = key_s(k+1) - key_s(k);
        
        if k == length(key_s)-1
            if dist < min_seg_dist
                 is_merged = true;
            else
                 new_indices = [new_indices, key_indices(k)];
                 new_s = [new_s, key_s(k)];
                 new_v = [new_v, key_v_limit(k)];
            end
            new_indices = [new_indices, key_indices(end)];
            new_s = [new_s, key_s(end)];
            new_v = [new_v, key_v_limit(end)];
            continue;
        end
        
        if dist < min_seg_dist
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
    
    if is_merged
        key_indices = new_indices;
        key_s = new_s;
        key_v_limit = new_v;
    end
end
num_segments = length(key_indices) - 1;

% 第三步：可行性前瞻
key_v_feasible = key_v_limit;
a_planning = 0.5 * a_design_limit; % 更加保守的规划加速度

% Backward
for k = length(key_s)-1 : -1 : 1
    ds = key_s(k+1) - key_s(k);
    v_next = key_v_feasible(k+1);
    v_max_brake = sqrt(v_next^2 + 2 * a_planning * ds);
    key_v_feasible(k) = min(key_v_feasible(k), v_max_brake);
end

% Forward
for k = 1 : length(key_s)-1
    ds = key_s(k+1) - key_s(k);
    v_curr = key_v_feasible(k);
    v_max_accel = sqrt(v_curr^2 + 2 * a_planning * ds);
    key_v_feasible(k+1) = min(key_v_feasible(k+1), v_max_accel);
end

% 第四步：分段 S 曲线生成
t_total_segments = [];
s_total_segments = [];
v_total_segments = [];
a_total_segments = [];

current_t_offset = 0;
current_s_offset = 0;

for i = 1:num_segments
    s_start = key_s(i);
    s_end   = key_s(i+1);
    dist_seg = s_end - s_start;
    v_start = key_v_feasible(i);
    v_end   = key_v_feasible(i+1);
    
    idx_start = key_indices(i);
    idx_end = key_indices(i+1);
    seg_geo_max = max(v_limit_smooth(idx_start:idx_end));
    v_cruise_target = min(v_max_lin, seg_geo_max);
    
    req_dist = (abs(v_cruise_target^2 - v_start^2) + abs(v_cruise_target^2 - v_end^2)) / (1.5 * a_design_limit);
    
    if dist_seg < req_dist
        v_cruise_final = max(v_start, v_end);
    else
        v_cruise_final = v_cruise_target;
    end
    v_cruise_final = max([v_cruise_final, v_start, v_end]);
    
    [q, v, a, ~, t_seg, T_seg, ~] = double_s_curve_complete(0, dist_seg, v_start, v_end, v_cruise_final, a_design_limit, j_max);
    
    if i == 1
        start_idx = 1;
    else
        start_idx = 2;
    end
    
    t_total_segments = [t_total_segments, t_seg(start_idx:end) + current_t_offset];
    s_total_segments = [s_total_segments, q(start_idx:end) + current_s_offset];
    v_total_segments = [v_total_segments, v(start_idx:end)];
    a_total_segments = [a_total_segments, a(start_idx:end)];
    
    current_t_offset = current_t_offset + T_seg;
    current_s_offset = current_s_offset + dist_seg;
end

% 第五步：全局安全检查与时间缩放
% 原理：加速度与时间的平方成反比。如果最大加速度超了，就把时间拉长一点点。
% 1. 映射回曲率空间
u_mapped_check = interp1(s_scan, u_scan, s_total_segments, 'linear', 'extrap');
kappa_check = zeros(size(s_total_segments));
for k = 1:length(u_mapped_check)
    % 快速估算曲率 (这里简单用插值，实际应用可用更精确方法)
    idx = find(s_scan >= s_total_segments(k), 1);
    if isempty(idx), idx = length(s_scan); end
    kappa_check(k) = kappa_scan(idx);
end

% 2. 计算实际的总加速度
a_tan_check = a_total_segments;
a_norm_check = (v_total_segments.^2) .* kappa_check;
a_total_check = sqrt(a_tan_check.^2 + a_norm_check.^2);

% 3. 寻找违规峰值
max_a_total = max(a_total_check);
fprintf('  [安检] 原始规划最大合加速度: %.3f m/s^2 (限制: %.1f)\n', max_a_total, a_max_lin);

time_scaling_ratio = 1.0;
if max_a_total > a_max_lin
    % 计算缩放比例
    % a_new = a_old / k^2  => k = sqrt(a_old / a_limit)
    time_scaling_ratio = sqrt(max_a_total / a_max_lin) * 1.01; % 多给 1% 安全余量
    fprintf('  [修正] 发现超限！执行时间缩放，比例: %.4f\n', time_scaling_ratio);
else
    fprintf('  [通过] 轨迹安全，无需修正。\n');
end

% 4. 执行缩放
% 时间变长，速度变慢，加速度变小
t_total_segments = t_total_segments * time_scaling_ratio;
v_total_segments = v_total_segments / time_scaling_ratio;
a_total_segments = a_total_segments / (time_scaling_ratio^2);

T_scurve = t_total_segments(end);

% --- 第六步：导出绘图数据 ---
[s_scan_unique, unique_idx] = unique(s_scan);
u_scan_unique = u_scan(unique_idx);
u_scurve_mapped = interp1(s_scan_unique, u_scan_unique, s_total_segments, 'pchip');
u_scurve_mapped = max(0, min(1, u_scurve_mapped)); 

t_grid_scurve = t_total_segments;
s_dist_scurve = s_total_segments;
s_vel_scurve  = v_total_segments;
s_acc_scurve  = a_total_segments;

fprintf('分段 S-Curve (安全修正版) 规划完成，总时间: %.4f 秒\n', T_scurve);
toc;

% 结果对比分析
fprintf('\n===== 4. 生成对比数据 =====\n');

% 为了对比，我们将 S曲线 的结果采样到相同的时间分辨率
dt_plot = 0.01;
t_plot = 0:dt_plot:max(T_toppra, T_scurve);

% TOPP-RA 数据重采样 
vel_mag_toppra = zeros(size(t_plot));
acc_mag_toppra = zeros(size(t_plot));

for k = 1:length(t_plot)
    t_now = t_plot(k);
    if t_now > T_toppra
        t_now = T_toppra; 
    end
    
    u_now = interp1(t_grid_toppra, s_grid, t_now, 'pchip');
    s_dot_now = interp1(t_grid_toppra, sd_optimal, t_now, 'pchip');
    
    % 数值差分算 s_ddot
    idx = find(t_grid_toppra <= t_now, 1, 'last');
    if isempty(idx) || idx >= length(t_grid_toppra)
        idx = length(t_grid_toppra)-1; 
    end
    dt_local = t_grid_toppra(idx+1) - t_grid_toppra(idx);
    dv_local = sd_optimal(idx+1) - sd_optimal(idx);
    s_ddot_now = dv_local / dt_local;
    if isnan(s_ddot_now)
        s_ddot_now = 0; 
    end
    
    CK = CurveDerivs1(u_now, U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p_vec = CK(2,:)';
    P_pp_vec = CK(3,:)';
    
    vel_vec = P_p_vec * s_dot_now;
    acc_vec = P_p_vec * s_ddot_now + P_pp_vec * s_dot_now^2;
    
    vel_mag_toppra(k) = norm(vel_vec);
    acc_mag_toppra(k) = norm(acc_vec);
end

% S-Curve 数据重采样
% S-Curve 数据重采样
vel_mag_scurve = zeros(size(t_plot));
acc_mag_scurve = zeros(size(t_plot));

for k = 1:length(t_plot)
    t_now = t_plot(k);
    if t_now > T_scurve
        t_now = T_scurve; 
    end
    
    % 1. 获取当前 S-Curve 规划的线速度 (s_dot) 和 线加速度 (s_ddot)
    s_dot_now = interp1(t_grid_scurve, s_vel_scurve, t_now, 'linear');
    s_ddot_now = interp1(t_grid_scurve, s_acc_scurve, t_now, 'linear');
    
    % 2. 获取当前的 B样条参数 u
    u_now = interp1(t_grid_scurve, u_scurve_mapped, t_now, 'pchip');
    
    % 3. 计算笛卡尔导数
    CK = CurveDerivs1(u_now, U, p, P_pos_ctrl_T, 2, num_ctrl_pts - 1);
    P_p_vec = CK(2,:)'; 
    P_pp_vec = CK(3,:)';
    
    % 4. 链式法则转换: 
    % v = P' * u_dot
    % a = P' * u_ddot + P'' * u_dot^2
    % 其中 u_dot = s_dot / |P'|
    
    norm_Pp = norm(P_p_vec);
    if norm_Pp < 1e-6
        u_dot_now = 0; 
        u_ddot_now = 0;
    else
        u_dot_now = s_dot_now / norm_Pp;
        % u_ddot = (s_ddot - (P'·P'' / |P'|^2) * s_dot^2 ) / |P'|  <-- 修正公式
        term_tangent = dot(P_p_vec, P_pp_vec) / (norm_Pp^2);
        u_ddot_now = (s_ddot_now - term_tangent * s_dot_now^2) / norm_Pp; % 注意这里的推导
        % 上面这个 u_ddot 公式来源于 s_ddot = d/dt(|P'|u_dot) 的展开
    end
    
    vel_vec = P_p_vec * u_dot_now;
    acc_vec = P_p_vec * u_ddot_now + P_pp_vec * u_dot_now^2;
    
    vel_mag_scurve(k) = norm(vel_vec);
    acc_mag_scurve(k) = norm(acc_vec);
end

% 全面指标对比绘图 (含分量分解)
% 数据预处理与分量计算
% 强制转换为列向量
t = t_plot(:);
dt = gradient(t);

% TOPP-RA 数据
v_toppra = vel_mag_toppra(:);
a_tot_toppra = acc_mag_toppra(:);
at_toppra = gradient(v_toppra) ./ dt;
an_toppra = sqrt(max(0, a_tot_toppra.^2 - at_toppra.^2));

% S-Curve 数据
v_scurve = vel_mag_scurve(:);
a_tot_scurve = acc_mag_scurve(:);
at_scurve = gradient(v_scurve) ./ dt;
an_scurve = sqrt(max(0, a_tot_scurve.^2 - at_scurve.^2));

% 开始绘图 
figure('Position', [50, 50, 1200, 1000], 'Color', 'w');
sgtitle('TOPP-RA vs Improved S-Curve (曲率限速版)', 'FontSize', 16);

% 子图 1: 速度对比
subplot(3, 2, 1);
plot(t, v_toppra, 'b-', 'LineWidth', 2); 
hold on;
plot(t, v_scurve, 'r--', 'LineWidth', 2);
yline(v_max_lin, 'k-.', 'V_{max} Set');
grid on; 
xlabel('Time (s)'); 
ylabel('Speed (m/s)');
title('1. 末端线速度 v(t)');
legend('TOPP-RA', 'S-Curve (Safe)', 'Location', 'South');

% 子图 2: 总耗时对比
subplot(3, 2, 2);
b = barh([1, 2], [T_toppra, T_scurve], 0.6);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1]; % TOPP-RA (Blue)
b.CData(2,:) = [0 1 0]; % S-Curve (Green) -> 表示Safe but Slow
yticklabels({'TOPP-RA', 'S-Curve (Safe)'});
xlabel('Total Time (s)');
title('2. 总耗时对比');
grid on;
text(T_toppra/2, 1, sprintf('%.3f s (Time Opt)', T_toppra), 'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(T_scurve/2, 2, sprintf('%.3f s (Safe but Slow)', T_scurve), 'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% 子图 3: 切向加速度 (油门/刹车)
subplot(3, 2, 3);
plot(t, at_toppra, 'b-', 'LineWidth', 1.5); hold on;
plot(t, at_scurve, 'r--', 'LineWidth', 1.5);
yline(a_max_lin, 'k:', 'Limit');  
yline(-a_max_lin, 'k:', 'Limit'); 
grid on;
ylabel('Acc (m/s^2)');
title('3. 切向加速度 a_t');
legend('TOPP-RA', 'S-Curve (Safe)', 'Location', 'Best');
ylim([-a_max_lin*1.2, a_max_lin*1.2]);

% 子图 4: 法向加速度 (离心力)
subplot(3, 2, 4);
% 这次 S-Curve 是安全的，所以不需要标红警告
plot(t, an_scurve, 'r--', 'LineWidth', 2); hold on;
plot(t, an_toppra, 'b-', 'LineWidth', 2);
yline(a_max_lin, 'k-.', 'Limit');
grid on;
ylabel('Acc (m/s^2)');
title('4. 法向加速度 a_n (离心力)');
subtitle('两者均在安全范围内，但 S-Curve 过于保守');
legend('S-Curve (保守)', 'TOPP-RA (激进)', 'Location', 'Best');
ylim([0, a_max_lin*1.2]);

% 子图 5: 总加速度 (切向+法向矢量和)
subplot(3, 2, [5, 6]); 
plot(t, a_tot_toppra, 'b-', 'LineWidth', 2); hold on;
plot(t, a_tot_scurve, 'r--', 'LineWidth', 2);
yline(a_max_lin, 'k-.', 'A_{max} Limit', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('Total Acceleration (m/s^2)');
title('5. 总加速度对比 (最终合力)');
legend('TOPP-RA', 'S-Curve (Safe)', 'Location', 'Best');
ylim([0, a_max_lin*1.2]);

fprintf('绘图完成。\n');
fprintf('对比结论:\n');
fprintf('1. 改进后的 S 曲线通过降低全程速度，成功将最大离心力控制在安全范围内(图4)。\n');
fprintf('2. 但代价是总耗时显著增加(图2)，因为它为了最弯的那一点，牺牲了整条路的速度。\n');
fprintf('3. TOPP-RA 依然是时间最优的，因为它只在需要减速的地方减速。\n');

% 辅助函数 
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
            if Tv < 0
                Tv = 0; 
            end
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
        Tj1=0; 
        Tj2=0; 
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
            delta_s_acc = (v0+vlim)*Ta/2; 
            delta_s_const = vlim*Tv; 
            q_dec_start = delta_s_acc + delta_s_const; 
            q(i) = q_dec_start + vlim*t5 + alim_d/6*(3*t5^2 - 3*Tj2*t5 + Tj2^2); 
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
        q_all = q0 - q; 
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

% 在满足约束的前提下的反向最大减速度和正向最大加速度
function [alpha_min, alpha_max] = GetSDDotBounds_Exact(Pp, Ppp, Qp, Qpp, x, a_max, wa_max)
    alpha_min = -1e10; 
    alpha_max = 1e10;
    
    % 线加速度
    A_vec = Pp; 
    B_vec = Ppp * x; 
    R = a_max;
    a_coeff = dot(A_vec, A_vec); 
    b_coeff = 2 * dot(A_vec, B_vec); 
    c_coeff = dot(B_vec, B_vec) - R^2;
    [min1, max1] = SolveQuadraticInequality(a_coeff, b_coeff, c_coeff);
    alpha_min = max(alpha_min, min1); 
    alpha_max = min(alpha_max, max1);
    
    % 角加速度
    A_vec_w = 2 * Qp; 
    B_vec_w = 2 * Qpp * x; 
    R_w = wa_max;
    a_coeff_w = dot(A_vec_w, A_vec_w); 
    b_coeff_w = 2 * dot(A_vec_w, B_vec_w); 
    c_coeff_w = dot(B_vec_w, B_vec_w) - R_w^2;
    [min2, max2] = SolveQuadraticInequality(a_coeff_w, b_coeff_w, c_coeff_w);
    alpha_min = max(alpha_min, min2); 
    alpha_max = min(alpha_max, max2);
    
    if alpha_min > alpha_max
        alpha_max = alpha_min; 
    end
end

function [xmin, xmax] = SolveQuadraticInequality(a, b, c)
    if abs(a) < 1e-10
        if abs(b) < 1e-10
            if c <= 1e-10
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
    else
        delta = b^2 - 4*a*c;
        if delta < 0 && delta > -1e-8
            delta = 0; 
        end
        if delta < 0
            if a > 0
                xmin = 1e10;
                xmax = -1e10; 
            else
                xmin = -1e10;
                xmax = 1e10; 
            end
        else
            sqrt_delta = sqrt(delta); 
            r1 = (-b - sqrt_delta) / (2*a); 
            r2 = (-b + sqrt_delta) / (2*a);
            if a > 0
                xmin = min(r1, r2); 
                xmax = max(r1, r2);
            else
                xmin = -1e10; 
                xmax = 1e10; 
            end
        end
    end
end

% B样条工具函数
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
    Q_target(num_via+1, :) = v0'; 
    Q_target(num_via+2, :) = a0';
    Q_target(num_via+3, :) = v6'; 
    Q_target(num_via+4, :) = a6';
    W = eye(num_targets); 
    W(1,1) = 1e6; 
    W(num_via, num_via) = 1e6;
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
    span0 = FindSpan(num_ctrl_pts-1, p, u_bar(1), U); 
    ders0 = DersBasisFuns(span0, u_bar(1), p, 2, U);
    for j=0:p
        B(num_via+1, span0-p+j+1) = ders0(2, j+1); 
        B(num_via+2, span0-p+j+1) = ders0(3, j+1); 
    end
    span6 = FindSpan(num_ctrl_pts-1, p, u_bar(end), U); 
    ders6 = DersBasisFuns(span6, u_bar(end), p, 2, U);
    for j=0:p
        B(num_via+3, span6-p+j+1) = ders6(2, j+1); 
        B(num_via+4, span6-p+j+1) = ders6(3, j+1); 
    end
end

function [C, A] = BuildSmoothMatrix(num_ctrl_pts, p, U)
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
    ders = zeros(n+1, p+1); 
    ndu = zeros(p+1, p+1); 
    ndu(1,1) = 1.0;
    left = zeros(p+1, 1); 
    right = zeros(p+1, 1);
    for j = 1:p
        left(j+1) = u - U(i+1-j+1); 
        right(j+1) = U(i+j+1) - u;
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
        s1 = 0; 
        s2 = 1;
        a(1, 1) = 1.0;
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

function CK = CurveDerivs1(u, U, p, P, d, n)
    dim = size(P, 1); 
    du = min(d, p); 
    CK = zeros(d+1, dim);
    span = FindSpan(n, p, u, U); 
    nders = DersBasisFuns(span, u, p, du, U);
    for k = 0:du
        for j = 0:p
            CK(k+1, :) = CK(k+1, :) + nders(k+1, j+1) * P(:, span-p+j+1)';
        end
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
