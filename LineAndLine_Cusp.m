clear; clc; close all;

% 途径点位置 (8个点) [x, y, z]
Q_pos = [0.0 0.0 0.0;
         0.0 0.4 0.0;
         0.6 0.4 0.0;
         0.6 0.0 0.0;
         0.6 0.0 0.4;
         0.0 0.0 0.4;
         0.0 0.4 0.4;
         0.0 0.4 0.0];

% 途径点姿态 (8个点) [Roll, Pitch, Yaw]
Q_eul = [0, 0, 0;
         0, 0, pi/4;
         0, 0, pi/2;
         pi/4, 0, pi/2;
         pi/2, 0, pi/2;
         pi/2, 0, 0;
         0, pi/4, 0;
         0, 0, 0];

% 将欧拉角转换为四元数
num_pts = size(Q_pos, 1);
Q_quat = zeros(num_pts, 4);
for i = 1:num_pts
    Q_quat(i,:) = my_eul2quat(Q_eul(i,:));
end

delta = 0.02;   % 切角半径
v_max_global = 1.0;  % 全局最大速度
a_max = 20.0; 
j_max = 100.0; 

% 几何路径生成 (位置+姿态)
segments = struct('type', {}, 'length', {}, 'max_kappa', {}, 'v_limit', {}, 'p_start', {}, 'p_end', {}, 'vec', {}, 'CtrlPts', {}, ...
                  'quat_start', {}, 'quat_end', {});
cnt = 1;

% 计算切点
q_in_pos = zeros(num_pts, 3);
q_out_pos = zeros(num_pts, 3);
q_in_quat = zeros(num_pts, 4);
q_out_quat = zeros(num_pts, 4);

for i = 2:num_pts-1
    vec_prev = Q_pos(i,:) - Q_pos(i-1,:);
    len_prev = norm(vec_prev);
    vec_next = Q_pos(i+1,:) - Q_pos(i,:);
    len_next = norm(vec_next);
    dir_prev = vec_prev / len_prev;  % 切入圆弧的方向
    dir_next = vec_next / len_next;  % 切出圆弧的方向
    
    d = min([delta, len_prev/2, len_next/2]); 
    
    q_in_pos(i,:) = Q_pos(i,:) - d * dir_prev;  % 切入圆弧的起点(位置)
    q_out_pos(i,:) = Q_pos(i,:) + d * dir_next; % 切出圆弧的终点(位置)
    
    ratio_in = (len_prev - d) / len_prev;
    q_in_quat(i,:) = my_slerp(Q_quat(i-1,:), Q_quat(i,:), ratio_in);  % 切入圆弧的起点(姿态)
    
    ratio_out = d / len_next;
    q_out_quat(i,:) = my_slerp(Q_quat(i,:), Q_quat(i+1,:), ratio_out);% 切出圆弧的终点(姿态)
end

% 构建路径链表
current_pos = Q_pos(1,:);
current_quat = Q_quat(1,:);

for i = 1:num_pts-1
    seg = struct('type', [], 'length', [], 'max_kappa', [], 'v_limit', [], 'p_start', [], 'p_end', [], 'vec', [], 'CtrlPts', [], ...
                 'quat_start', [], 'quat_end', []);
    
    % 直线段
    if i == 1
        target_pos = q_in_pos(2,:); 
        target_quat = q_in_quat(2,:);
    elseif i == num_pts-1
        target_pos = Q_pos(end,:); 
        target_quat = Q_quat(end,:);
    else
        target_pos = q_in_pos(i+1,:); 
        target_quat = q_in_quat(i+1,:);
    end
    
    seg.type = 'Line';
    seg.p_start = current_pos; 
    seg.p_end = target_pos;
    seg.vec = target_pos - current_pos;
    seg.length = norm(seg.vec);
    seg.max_kappa = 0; 
    seg.v_limit = v_max_global;
    seg.quat_start = current_quat; 
    seg.quat_end = target_quat;
    segments(cnt) = seg; 
    cnt = cnt + 1;
    
    current_pos = target_pos; 
    current_quat = target_quat;
    
    % 贝塞尔曲线段 
    if i < num_pts-1
        seg = struct('type', [], 'length', [], 'max_kappa', [], 'v_limit', [], 'p_start', [], 'p_end', [], 'vec', [], 'CtrlPts', [], ...
                     'quat_start', [], 'quat_end', []);
                 
        P0 = q_in_pos(i+1,:); 
        P5 = q_out_pos(i+1,:);
        dist = norm(P5 - P0);
        dir_in = (Q_pos(i+1,:) - Q_pos(i,:)) / norm(Q_pos(i+1,:) - Q_pos(i,:));       % 圆弧切入方向
        dir_out = (Q_pos(i+2,:) - Q_pos(i+1,:)) / norm(Q_pos(i+2,:) - Q_pos(i+1,:));  % 圆弧切出方向
        
        P1 = P0 + 0.2 * dist * dir_in; 
        P2 = P1; 
        P4 = P5 - 0.2 * dist * dir_out; 
        P3 = P4;
        CtrlPts = [P0; P1; P2; P3; P4; P5];
        
        [L_arc, kappa_max] = analyze_bezier(CtrlPts);
        
        seg.type = 'Bezier';
        seg.CtrlPts = CtrlPts;
        seg.length = L_arc;
        seg.max_kappa = kappa_max;
        
        if kappa_max > 1e-6
            v_curvature = sqrt(a_max / kappa_max);
        else
            v_curvature = v_max_global;
        end
        seg.v_limit = min(v_max_global, v_curvature);
        seg.quat_start = current_quat; 
        seg.quat_end = q_out_quat(i+1,:);
        
        segments(cnt) = seg; 
        cnt = cnt + 1;
        current_pos = P5; 
        current_quat = seg.quat_end;
    end
end

% 分段速度合成 
num_segs = length(segments);
segment_v_constr = zeros(1, num_segs);

for k = 1:num_segs  % 速度限制
    if strcmp(segments(k).type, 'Bezier')
        segment_v_constr(k) = segments(k).v_limit;
    else
        segment_v_constr(k) = v_max_global; 
    end
end

trajectory_data = [];
t_current = 0;

for k = 1:num_segs
    seg = segments(k);
    
    if strcmp(seg.type, 'Bezier')
        % 弯道
        v_run = segment_v_constr(k);
        T_seg = seg.length / v_run;
        dt = 0.01; 
        t_local = 0:dt:T_seg;
        if isempty(t_local)
            t_local = [0, T_seg]; 
        end
        
        s_local = v_run * t_local; 
        u_local = s_local / seg.length; 
        
        pos = zeros(length(u_local), 3);
        quat = zeros(length(u_local), 4);
        acc_mag = zeros(1, length(u_local)); % 存储总加速度
        
        for m = 1:length(u_local)
            % 获取贝塞尔的一阶导 d1 (切向) 和二阶导 d2 (弯曲程度)
            [p, d1, d2] = eval_bezier(seg.CtrlPts, u_local(m));
            pos(m,:) = p;
            
            % 计算瞬时曲率 k = |d1 x d2| / |d1|^3
            c_cross = cross(d1, d2);
            num = sqrt(sum(c_cross.^2));
            den = (sqrt(sum(d1.^2)))^3 + 1e-9;
            kappa = num / den;
            
            % 计算总加速度
            % 弯道是匀速运动，所以切向加速度=0
            % 只有法向加速度(离心力) a = v^2 * kappa
            acc_mag(m) = v_run^2 * kappa;
            
            quat(m,:) = my_slerp(seg.quat_start, seg.quat_end, u_local(m));
        end
        
        chunk.t = t_current + t_local;
        chunk.pos = pos;
        chunk.quat = quat;
        chunk.v_mag = repmat(v_run, size(t_local));
        chunk.a_mag = acc_mag; % 保存加速度
        
        trajectory_data = [trajectory_data, chunk];
        t_current = t_current + T_seg;
        
    else
        % 直线
        if k == 1
            v_start = 0; 
        else
            v_start = segments(k-1).v_limit; 
        end
        if k == num_segs
            v_end = 0; 
        else
            v_end = segments(k+1).v_limit; 
        end
        
        q0 = 0; 
        q1 = seg.length;
        % 注意这里要把 sdd (加速度) 也接出来
        [s, sd, sdd, ~, t_local, T_seg] = double_s_curve_complete(q0, q1, v_start, v_end, v_max_global, a_max, j_max);
        
        pos = zeros(length(t_local), 3);
        quat = zeros(length(t_local), 4);
        dir_vec = seg.vec / seg.length;
        
        for m = 1:length(t_local)
            pos(m,:) = seg.p_start + dir_vec * s(m);
            if seg.length > 1e-6
                u = s(m) / seg.length; 
            else
                u = 1; 
            end
            u = max(0, min(1, u));
            quat(m,:) = my_slerp(seg.quat_start, seg.quat_end, u);
        end
        
        chunk.t = t_current + t_local;
        chunk.pos = pos;
        chunk.quat = quat;
        chunk.v_mag = sd;
        chunk.a_mag = abs(sdd); % 直线加速度就是 Double S 算出来的切向加速度
        
        trajectory_data = [trajectory_data, chunk];
        t_current = t_current + T_seg;
    end
end

% 绘图展示 (3个子图: 3D路径, 速度, 加速度)
fprintf('总耗时: %.4f s\n', t_current);

% 提取合并数据
T_all = [trajectory_data.t];
P_all = vertcat(trajectory_data.pos);
Q_all = vertcat(trajectory_data.quat);
V_all = [trajectory_data.v_mag];
A_all = [trajectory_data.a_mag];

figure('Color','w','Position',[100 50 1200 800]);

% 子图1: 3D 轨迹与姿态
subplot(2, 2, [1, 3]); % 占据左半边
plot3(Q_pos(:,1), Q_pos(:,2), Q_pos(:,3), 'ko--', 'LineWidth', 1, 'MarkerFaceColor', 'k'); 
hold on;
plot3(P_all(:,1), P_all(:,2), P_all(:,3), 'b-', 'LineWidth', 2);

% 绘制坐标系
step = 25; 
scale = 0.08; 
for i = 1:step:size(P_all, 1)
    p = P_all(i,:); 
    q = Q_all(i,:);
    R = my_quat2rotm(q);
    quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), scale, 'r', 'LineWidth', 1.5);
    quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), scale, 'g', 'LineWidth', 1.5);
    quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), scale, 'b', 'LineWidth', 1.5);
end
grid on; axis equal; view(3);
title('3D Trajectory with Orientation');
xlabel('X'); ylabel('Y'); zlabel('Z');

% 子图2: 速度曲线
subplot(2, 2, 2);
plot(T_all, V_all, 'k-', 'LineWidth', 1.5);
yline(v_max_global, 'r--', 'V max');
grid on;
title('End-Effector Speed Profile');
xlabel('Time (s)'); 
ylabel('Speed (m/s)');
ylim([0, v_max_global * 1.2]);

% 子图3: 加速度曲线
subplot(2, 2, 4);
plot(T_all, A_all, 'b-', 'LineWidth', 1.5);
yline(a_max, 'r--', 'A max limit');
grid on;
title('End-Effector Acceleration Profile (Total)');
xlabel('Time (s)'); 
ylabel('Accel (m/s^2)');
ylim([0, max(A_all) * 1.2 + 0.1]);


% 辅助函数 
% 欧拉角转为四元数
function q = my_eul2quat(eul)
    phi = eul(1)/2; theta = eul(2)/2; psi = eul(3)/2;
    c1 = cos(phi); s1 = sin(phi); c2 = cos(theta); s2 = sin(theta); c3 = cos(psi); s3 = sin(psi);

    w = c1*c2*c3 + s1*s2*s3; x = s1*c2*c3 - c1*s2*s3; y = c1*s2*c3 + s1*c2*s3; z = c1*c2*s3 - s1*s2*c3;
    q = [w, x, y, z]; q = q / norm(q);
end

% 姿态的球面线性插值
function q_interp = my_slerp(q1, q2, t)
    dot_val = sum(q1 .* q2);
    if dot_val < 0
        q1 = -q1; 
        dot_val = -dot_val; 
    end
    if dot_val > 0.9995 
        q_interp = (1-t)*q1 + t*q2; 
        q_interp = q_interp / norm(q_interp); 
        return; 
    end
    theta_0 = acos(dot_val); 
    theta = theta_0 * t;
    s0 = cos(theta) - dot_val * sin(theta) / sin(theta_0); 
    s1 = sin(theta) / sin(theta_0);
    q_interp = s0 * q1 + s1 * q2; 
    q_interp = q_interp / norm(q_interp);
end

% 四元数转化为旋转矩阵
function R = my_quat2rotm(q)
    w = q(1); 
    x = q(2); 
    y = q(3); 
    z = q(4);
    xx = x*x; 
    yy = y*y; 
    zz = z*z; 
    xy = x*y; 
    xz = x*z; 
    yz = y*z; 
    wx = w*x; 
    wy = w*y; 
    wz = w*z;
    R = [1 - 2*(yy+zz), 2*(xy-wz), 2*(xz+wy); 2*(xy+wz), 1 - 2*(xx+zz), 2*(yz-wx); 2*(xz-wy), 2*(yz+wx), 1 - 2*(xx+yy)];
end

% 计算曲线曲率最大处
function [len, k_max] = analyze_bezier(CP)
    u = linspace(0, 1, 100); 
    pts = zeros(100, 3); 
    d1 = zeros(100, 3); 
    d2 = zeros(100, 3);
    for i=1:100
        [pts(i,:), d1(i,:), d2(i,:)] = eval_bezier(CP, u(i)); 
    end
    segs = diff(pts, 1, 1); 
    len = sum(sqrt(sum(segs.^2, 2)));
    k_vals = zeros(1, 100);
    for i=1:100
        nom = norm(cross(d1(i,:), d2(i,:))); 
        den = norm(d1(i,:))^3 + 1e-9; 
        k_vals(i) = nom / den;
    end
    k_max = max(k_vals);
end

% 计算贝塞尔曲线
function [p, d1, d2] = eval_bezier(P, u)
    B = [(1-u)^5, 5*u*(1-u)^4, 10*u^2*(1-u)^3, 10*u^3*(1-u)^2, 5*u^4*(1-u), u^5];
    B_d = [-5*(1-u)^4, 5*(1-u)^4 - 20*u*(1-u)^3, 20*u*(1-u)^3 - 30*u^2*(1-u)^2, 30*u^2*(1-u)^2 - 20*u^3*(1-u), 20*u^3*(1-u) - 5*u^4, 5*u^4];
    p = B * P; 
    d1 = B_d * P;  % 一阶导
    B3 = [(1-u)^3, 3*u*(1-u)^2, 3*u^2*(1-u), u^3];
    Q_dd = zeros(4, 3);  % 二阶差分控制点
    for i=1:4
        Q_dd(i,:) = P(i+2,:) - 2*P(i+1,:) + P(i,:); 
    end
    d2 = 20 * B3 * Q_dd; % 二阶导
end

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
            % t_in_dec = t - (Ta + Tv); 
            q(i) = q_dec_start + vlim*t5 + alim_d/6*(3*t5^2 - 3*Tj2*t5 + Tj2^2); 
            v(i) = vlim + alim_d*(t5 - Tj2/2); 
            a(i) = alim_d;
        else
            t6 = T - t; 
            q(i) = dist - v1*t6 - jmax*t6^3/6; 
            v(i) = v1 + jmax*t6^2/2; 
            a(i) = -jmax*t6; 
        end
        q(i) = q(i) + 0;
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
        Td = 2 * Tj2; alim_d = -jmax * Tj2 * sign(v_target - v1); 
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
