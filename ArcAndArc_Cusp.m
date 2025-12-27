clear; clc; close all;

% 圆弧1: 水平进场
Arc1.Center = [-1, 0, 0];
Arc1.Radius = 1.0;
Arc1.Normal = [0, 0, 1]; 
Arc1.StartAngle = -pi/2; 
Arc1.EndAngle = 0;       
Arc1.Dir = 1; % 1=逆时针

% 圆弧2: 垂直向上出场
Arc2.Center = [0, 1, 0];
Arc2.Radius = 1.0;
Arc2.Normal = [0, 0, 1]; 
Arc2.StartAngle = -pi/2; 
Arc2.EndAngle = 0; 
Arc2.Dir = 1; 

% 姿态定义
Quat_Start = my_eul2quat([0, 0, 0]);
Quat_End   = my_eul2quat([0, pi/4, pi/2]);

% 物理限制参数
delta = 0.08;        % 过渡区半长
v_max_global = 0.8;  % 全局最大速度
a_max = 5.0;        
j_max = 20.0;       

% 几何运算，生成平滑路径
fprintf('正在计算 G2 连续过渡曲线...\n');

% 计算切入点回退角度 (P_in) 和 切出点多余角度 (P_out) 
d_theta1 = delta / Arc1.Radius;
angle_in = Arc1.EndAngle - Arc1.Dir * d_theta1;

P_in = Arc1.Center + Arc1.Radius * [cos(angle_in), sin(angle_in), 0]; % 起始切入点位置向量
t_in = [-sin(angle_in), cos(angle_in), 0] * Arc1.Dir;                 % 起始切线方向（求导）
n_in_vec = (Arc1.Center - P_in) / norm(Arc1.Center - P_in);           % 起始曲率方向

d_theta2 = delta / Arc2.Radius;
angle_out = Arc2.StartAngle + Arc2.Dir * d_theta2;

P_out = Arc2.Center + Arc2.Radius * [cos(angle_out), sin(angle_out), 0];  % 终点切出点位置向量
t_out = [-sin(angle_out), cos(angle_out), 0] * Arc2.Dir;                  % 终点切线方向
n_out_vec = (Arc2.Center - P_out) / norm(Arc2.Center - P_out);            % 终点曲率方向

% 求解5次贝塞尔控制点 
CtrlPts = QuinticBezierSolver_Core(P_in, P_out, t_in, t_out, Arc1.Radius, n_in_vec, Arc2.Radius, n_out_vec);

% 构建三段式路径链表 
template_seg = struct('type', [], ...
                      'length', [], ...
                      'v_limit', [], ...
                      'Center', [], ...      % 圆弧专用字段
                      'Radius', [], ...      % 圆弧专用字段
                      'AngleRange', [], ...  % 圆弧专用字段
                      'CtrlPts', []);        % 贝塞尔专用字段

% 预分配 3 个结构体
segments = repmat(template_seg, 1, 3);

% 第一段圆弧
segments(1).type = 'Arc';
segments(1).Center = Arc1.Center;
segments(1).Radius = Arc1.Radius;
segments(1).AngleRange = [Arc1.StartAngle, angle_in]; 
segments(1).length = abs(segments(1).AngleRange(2) - segments(1).AngleRange(1)) * Arc1.Radius;
segments(1).v_limit = v_max_global;

% 第二段贝塞尔
segments(2).type = 'Bezier';
segments(2).CtrlPts = CtrlPts;
[len_bez, k_max_bez] = analyze_bezier(CtrlPts);
segments(2).length = len_bez;

if k_max_bez > 1e-6
    v_curve = sqrt(a_max / k_max_bez); % 拐弯处恒定速度
else
    v_curve = v_max_global;
end
segments(2).v_limit = min(v_max_global, v_curve);
% segments(2).v_limit = v_max_global; 
% 第三段圆弧
segments(3).type = 'Arc';
segments(3).Center = Arc2.Center;
segments(3).Radius = Arc2.Radius;
segments(3).AngleRange = [angle_out, Arc2.EndAngle]; 
segments(3).length = abs(segments(3).AngleRange(2) - segments(3).AngleRange(1)) * Arc2.Radius;
segments(3).v_limit = v_max_global;

% Double S 速度规划
fprintf('正在进行速度与姿态规划...\n');

trajectory_data = [];
t_current = 0;
total_len = sum([segments.length]);
dist_accumulated = 0;

for k = 1:length(segments)
    seg = segments(k);
    
    % 确定起止速度 
    if k == 1
        v_start = 0; 
    else
        v_start = min(segments(k-1).v_limit, seg.v_limit); 
    end

    if k == length(segments)
        v_end = 0; 
    else
        v_end = min(segments(k+1).v_limit, seg.v_limit); 
    end
    
    v_cruise = seg.v_limit; % 每段的最大速度
    
    [s, v, a, ~, t_local, T_seg] = double_s_curve_complete(0, seg.length, v_start, v_end, v_cruise, a_max, j_max);
    
    num_steps = length(t_local);
    pos = zeros(num_steps, 3);
    quat = zeros(num_steps, 4);
    acc_total = zeros(1, num_steps); 
    
    for m = 1:num_steps
        u_progress = s(m) / seg.length; % 参数u
        
        if strcmp(seg.type, 'Arc')
            % 圆弧参数方程
            ang = seg.AngleRange(1) + u_progress * (seg.AngleRange(2) - seg.AngleRange(1));
            pos(m,:) = seg.Center + seg.Radius * [cos(ang), sin(ang), 0];
            kappa = 1 / seg.Radius;
            
        elseif strcmp(seg.type, 'Bezier')
            % 贝塞尔求值
            [p, d1, d2] = eval_bezier(seg.CtrlPts, u_progress);
            pos(m,:) = p;
            c_cross = cross(d1, d2);
            kappa = norm(c_cross) / (norm(d1)^3 + 1e-9);
        end
        
        % 总加速度
        a_tan = a(m);
        a_norm = v(m)^2 * kappa;
        acc_total(m) = sqrt(a_tan^2 + a_norm^2);
        
        % 姿态插值
        global_progress = (dist_accumulated + s(m)) / total_len;
        global_progress = max(0, min(1, global_progress));
        quat(m,:) = my_slerp(Quat_Start, Quat_End, global_progress);
    end
    
    chunk.t = t_current + t_local;
    chunk.pos = pos;
    chunk.quat = quat;
    chunk.v = v;
    chunk.a = acc_total;
    
    trajectory_data = [trajectory_data, chunk];
    t_current = t_current + T_seg;
    dist_accumulated = dist_accumulated + seg.length;
end

% 结果可视化 
figure('Color','w','Position',[50 50 1400 800]);

% 左图：3D 轨迹 
subplot(2, 2, [1, 3]);
hold on; 
grid on; 
axis equal; 
view(3);
P_all = vertcat(trajectory_data.pos);
Q_all = vertcat(trajectory_data.quat);
plot3(P_all(:,1), P_all(:,2), P_all(:,3), 'b-', 'LineWidth', 2);
% plot3([Arc1.Center(1)+1, 0, Arc2.Center(1)], [0, 0, 1], [0,0,0], 'k--', 'LineWidth', 0.5);
scatter3(0,0,0, 50, 'r', 'filled'); 
text(0.1,0,0.1,'圆弧交点处');

step = 25; 
scale = 0.2;
for i=1:step:size(P_all,1)
    R = my_quat2rotm(Q_all(i,:)); 
    p = P_all(i,:);
    quiver3(p(1),p(2),p(3),R(1,1),R(2,1),R(3,1), scale, 'r');
    quiver3(p(1),p(2),p(3),R(1,2),R(2,2),R(3,2), scale, 'g');
    quiver3(p(1),p(2),p(3),R(1,3),R(2,3),R(3,3), scale, 'b');
end
title('Arc-Bezier-Arc Smooth Trajectory'); 
xlabel('X'); ylabel('Y'); zlabel('Z');

% 数据提取 
t_all = [trajectory_data.t];
v_all = [trajectory_data.v];
a_all = [trajectory_data.a];

% 右上图：速度曲线
subplot(2, 2, 2);
plot(t_all, v_all, 'b-', 'LineWidth', 1.5);
yline(v_max_global, 'b--', 'V Limit', 'LabelHorizontalAlignment','left');
grid on;
title('Velocity Profile');
ylabel('Speed (m/s)');
xlim([0, max(t_all)]);
ylim([0, v_max_global * 1.2]); % 留点余量

% 右下图：加速度曲线
subplot(2, 2, 4);
plot(t_all, a_all, 'r-', 'LineWidth', 1.5);
yline(a_max, 'r--', 'A Limit', 'LabelHorizontalAlignment','left');
grid on;
title('Total Acceleration Profile');
xlabel('Time (s)');
ylabel('Accel (m/s^2)');
xlim([0, max(t_all)]);
ylim([0, max(max(a_all), a_max) * 1.2]); 

% 辅助函数
% alpha法计算贝塞尔曲线控制点
function P_all = QuinticBezierSolver_Core(P0, P5, t0, t5, R0, n0, R5, n5)
    % 先用标准的 Alpha 法计算最佳切长系数 alpha
    vec_t_sum = t0 + t5; 
    vec_q_diff = P5 - P0; 
    val_a = 256 - 49 * sum(vec_t_sum.^2); 
    val_b = 420 * dot(vec_q_diff, vec_t_sum); 
    val_c = -900 * sum(vec_q_diff.^2);
    
    % 解一元二次方程
    delta_val = val_b^2 - 4 * val_a * val_c;
    if delta_val < 0
        delta_val = 0; 
    end
    alpha = (-val_b + sqrt(delta_val)) / (2 * val_a); 
    
    % 计算基础控制点 P1, P4 (切线上的点)
    P1 = P0 + (alpha/5)*t0; 
    P4 = P5 - (alpha/5)*t5; 
    
    % 计算 P2, P3 
    % 公式：P2 = (2*P1 - P0) + [法向偏移量]
    % 偏移量 = (alpha^2 / 20R) * n
    offset_start = (alpha^2 / (20 * R0)) * n0;
    P2 = 2*P1 - P0 + offset_start;
    
    offset_end   = (alpha^2 / (20 * R5)) * n5;
    P3 = 2*P4 - P5 + offset_end;
    
    P_all = [P0; P1; P2; P3; P4; P5];
end

% 双S曲线算法
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

% 寻找贝塞尔曲线的最大曲率
function [len, k_max] = analyze_bezier(CP)
    u = linspace(0, 1, 50); 
    pts = zeros(50, 3); 
    d1 = zeros(50, 3); 
    d2 = zeros(50, 3);
    for i=1:50
        [pts(i,:), d1(i,:), d2(i,:)] = eval_bezier(CP, u(i)); 
    end
    segs = diff(pts, 1, 1); 
    len = sum(sqrt(sum(segs.^2, 2)));

    k_vals = zeros(1, 50);
    for i = 1:50
        nom = norm(cross(d1(i,:), d2(i,:))); 
        den = norm(d1(i,:))^3 + 1e-9;
        k_vals(i) = nom / den; % 曲率公式
    end
    k_max = max(k_vals);
end

% 计算贝塞尔曲线及其导数
function [p, d1, d2] = eval_bezier(P, u)
    B = [(1-u)^5, 5*u*(1-u)^4, 10*u^2*(1-u)^3, 10*u^3*(1-u)^2, 5*u^4*(1-u), u^5];
    B_d = [-5*(1-u)^4, 5*(1-u)^4 - 20*u*(1-u)^3, 20*u*(1-u)^3 - 30*u^2*(1-u)^2, 30*u^2*(1-u)^2 ...
           - 20*u^3*(1-u), 20*u^3*(1-u) - 5*u^4, 5*u^4];
    p = B * P;    % 贝塞尔曲线
    d1 = B_d * P; % 一阶导
    B3 = [(1-u)^3, 3*u*(1-u)^2, 3*u^2*(1-u), u^3];
    Q_dd = zeros(4, 3); % 二阶导
    for i=1:4
        Q_dd(i,:) = P(i+2,:) - 2*P(i+1,:) + P(i,:); 
    end
    d2 = 20 * B3 * Q_dd; 
end

% 欧拉角转化为四元数
function q = my_eul2quat(eul)
    phi = eul(1)/2; 
    theta = eul(2)/2; 
    psi = eul(3)/2;
    q = [cos(phi)*cos(theta)*cos(psi)+sin(phi)*sin(theta)*sin(psi), ...
         sin(phi)*cos(theta)*cos(psi)-cos(phi)*sin(theta)*sin(psi), ...
         cos(phi)*sin(theta)*cos(psi)+sin(phi)*cos(theta)*sin(psi), ...
         cos(phi)*cos(theta)*sin(psi)-sin(phi)*sin(theta)*cos(psi)];
end

% 姿态球面线性插值
function q_interp = my_slerp(q1, q2, t)
    dot_val = sum(q1 .* q2);
    if dot_val < 0
        q1 = -q1; 
        dot_val = -dot_val; 
    end
    if dot_val > 0.9995
        q_interp = (1-t)*q1 + t*q2; 
        q_interp = q_interp/norm(q_interp); 
        return; 
    end
    theta_0 = acos(dot_val); 
    theta = theta_0 * t;
    q_interp = (sin(theta_0 - theta)*q1 + sin(theta)*q2) / sin(theta_0);
end

% 四元数转化为旋转矩阵
function R = my_quat2rotm(q)
    w=q(1); 
    x=q(2); 
    y=q(3); 
    z=q(4);
    R=[1-2*(y^2+z^2), 2*(x*y-z*w), 2*(x*z+y*w);
       2*(x*y+z*w), 1-2*(x^2+z^2), 2*(y*z-x*w);
       2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x^2+y^2)];
end
