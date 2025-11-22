clear; clc; close all;

p = 5;
q_via = [3, -2, -5, 0, 6, 12, 8;
        -1,  0,  2, 4,-9, 7,  3;
         0,  0,  0,-2,-1, 3,  0];
% q_via = [3, -2, -5, 0, 6, 12, 8;
%         -1,  0,  2, 4,-9, 7,  3;
%          0,  0,  0,0,0, 0,  0];

v0 = [-30, 10, 0]; 
v6 = [-20, 0, 0];
a0 = [-200, 10, 0]; 
a6 = [0, 300, 0]; 

% 弦长参数化
dists = sqrt(sum(diff(q_via,1,2).^2,1));
u_bar = [0, cumsum(dists) / sum(dists)];

% 关键参数
num_ctrl_pts = 16;  % 控制点个数
lambda = 5e-6;      % 平滑权重
num_via_pts = size(q_via, 2);
weights = ones(1, num_ctrl_pts);  % 初始权重全为1
weights(5) = 2.0;  % 增加某个控制点的权重

% 构建节点矢量 - 使用均匀分布
U = zeros(1, num_ctrl_pts + p + 1);
U(1:p+1) = 0;
U(end-p:end) = 1;
num_internal_knots = num_ctrl_pts - p - 1;
if num_internal_knots > 0
    internal_knots = linspace(0, 1, num_internal_knots + 2);
    internal_knots = internal_knots(2:end-1);
    U(p+2:num_ctrl_pts) = internal_knots;
end
fprintf('节点矢量 U = [');
fprintf('%.4f ', U);
fprintf(']');
fprintf('途径点参数 u_bar = [');
fprintf('%.4f ', u_bar);
fprintf(']\n');

% 构建约束
num_targets = 11;
Q_target = zeros(num_targets, 3);
Q_target(1:7, :) = q_via';
Q_target(8, :) = v0;
Q_target(9, :) = a0;
Q_target(10, :) = v6;
Q_target(11, :) = a6;

% 权重矩阵
W_diag = ones(num_targets, 1);
weight_boundary = 1e6;
boundary_indices = [1, 7, 8, 9, 10, 11];
W_diag(boundary_indices) = weight_boundary;
W = diag(W_diag);

% 构建B矩阵
B = zeros(num_targets, num_ctrl_pts);
% 位置约束（途径点）
for k = 1:7
    u_k = u_bar(k);
    span_k = FindSpan(num_ctrl_pts - 1, p, u_k, U);
    N_k = BasisFuns(span_k, u_k, p, U);
    for j = 0:p
        B(k, span_k - p + j + 1) = N_k(j + 1);
    end
    fprintf('途径点 %d: u=%.4f, span=%d\n', k-1, u_k, span_k);
end

% 起点速度和加速度约束
span0 = FindSpan(num_ctrl_pts - 1, p, u_bar(1), U);
ders0 = DersBasisFuns(span0, u_bar(1), p, 2, U);
for j = 0:p
    B(8, span0 - p + j + 1) = ders0(2, j + 1);  % 一阶导数
    B(9, span0 - p + j + 1) = ders0(3, j + 1);  % 二阶导数
end
fprintf('起点: u=%.4f, span=%d', u_bar(1), span0);

% 终点速度和加速度约束
span6 = FindSpan(num_ctrl_pts - 1, p, u_bar(end), U);
ders6 = DersBasisFuns(span6, u_bar(end), p, 2, U);
for j = 0:p
    B(10, span6 - p + j + 1) = ders6(2, j + 1);  % 一阶导数
    B(11, span6 - p + j + 1) = ders6(3, j + 1);  % 二阶导数
end
fprintf('终点: u=%.4f, span=%d\n', u_bar(end), span6);

% 构建平滑矩阵
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

% 求解
M = B' * W * B + lambda * (C' * A * C);
F = B' * W * Q_target;
P_calc = M \ F;
P_calc_T = P_calc';

% 验证边界条件
fprintf('===== 边界条件验证 =====');
CK_start = CurveDerivs1(u_bar(1), U, p, P_calc_T, 2, num_ctrl_pts - 1);
CK_end = CurveDerivs1(u_bar(end), U, p, P_calc_T, 2, num_ctrl_pts - 1);

fprintf('起点:\n');
fprintf('  位置: [%.3f, %.3f, %.3f], 目标: [%.3f, %.3f, %.3f], 误差: %.6e', ...
    CK_start(1,:), q_via(:,1)', norm(CK_start(1,:) - q_via(:,1)'));
fprintf('  速度: [%.3f, %.3f, %.3f], 目标: [%.3f, %.3f, %.3f], 误差: %.6e', ...
    CK_start(2,:), v0, norm(CK_start(2,:) - v0));
fprintf('  加速度: [%.3f, %.3f, %.3f], 目标: [%.3f, %.3f, %.3f], 误差: %.6e', ...
    CK_start(3,:), a0, norm(CK_start(3,:) - a0));

fprintf('\n终点:');
fprintf('  位置: [%.3f, %.3f, %.3f], 目标: [%.3f, %.3f, %.3f], 误差: %.6e', ...
    CK_end(1,:), q_via(:,end)', norm(CK_end(1,:) - q_via(:,end)'));
fprintf('  速度: [%.3f, %.3f, %.3f], 目标: [%.3f, %.3f, %.3f], 误差: %.6e', ...
    CK_end(2,:), v6, norm(CK_end(2,:) - v6));
fprintf('  加速度: [%.3f, %.3f, %.3f], 目标: [%.3f, %.3f, %.3f], 误差: %.6e\n', ...
    CK_end(3,:), a6, norm(CK_end(3,:) - a6));

% 绘制曲线
num_plot_pts = 500;
u_plot = linspace(0, 1, num_plot_pts);
curve_pts = zeros(3, num_plot_pts);
for k = 1:num_plot_pts
    curve_pts(:, k) = CurvePoint(num_ctrl_pts - 1, p, U, P_calc_T, u_plot(k)); % 无理
end

% num_plot_pts = 500;
% u_plot = linspace(0, 1, num_plot_pts);
curve_pts_uniform = zeros(3, num_plot_pts);
for k = 1:num_plot_pts
    curve_pts_uniform(:, k) = NURBSCurvePoint(num_ctrl_pts - 1, p, U, P_calc_T, weights, u_plot(k)); % 有理
end

figure('Position', [100, 100, 800, 600]);
subplot(1,2,1);
hold on;
plot3(curve_pts(1,:), curve_pts(2,:), curve_pts(3,:), 'b-', 'LineWidth', 2.5);
plot3(q_via(1,:), q_via(2,:), q_via(3,:), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
plot3(P_calc_T(1,:), P_calc_T(2,:), P_calc_T(3,:), 'k--s', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
for i = 1:num_via_pts
    text(q_via(1,i), q_via(2,i), q_via(3,i), sprintf('  q%d', i-1), 'FontSize', 10);
end
hold off;
title(sprintf('B样条曲线 (n=%d, \\lambda=%.0e)', num_ctrl_pts, lambda));
legend('曲线', '途径点', '控制点');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; axis equal;
view(120, 25);

subplot(1,2,2);
hold on;
plot3(curve_pts_uniform(1,:), curve_pts_uniform(2,:), curve_pts_uniform(3,:), 'b-', 'LineWidth', 2.5);
plot3(q_via(1,:), q_via(2,:), q_via(3,:), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
plot3(P_calc_T(1,:), P_calc_T(2,:), P_calc_T(3,:), 'k--s', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
for i = 1:num_via_pts
    text(q_via(1,i), q_via(2,i), q_via(3,i), sprintf('  q%d', i-1), 'FontSize', 10);
end
hold off;
title(sprintf('有理B样条曲线 (n=%d, \\lambda=%.0e)', num_ctrl_pts, lambda));
legend('曲线', '途径点', '控制点');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; axis equal;
view(120, 25);

% 计算速度
velocity_curve = zeros(num_plot_pts, 3);
for k = 1:num_plot_pts
    CK = CurveDerivs1(u_plot(k), U, p, P_calc_T, 1, num_ctrl_pts - 1);
    velocity_curve(k, :) = CK(2, :);
end

% subplot(1,2,2);
% plot(u_plot, vecnorm(velocity_curve, 2, 2), 'b-', 'LineWidth', 2);
% xlabel('参数 u'); ylabel('||dC/du||');
% title('速度模');
% grid on;

% 详细速度和加速度图
figure('Position', [100, 100, 1400, 800]);
for i = 1:3
    % 速度
    subplot(3,2,2*i-1);
    plot(u_plot, velocity_curve(:,i), 'LineWidth', 2);
    ylabel(sprintf('dC_%s/du', char('x'+i-1)));
    if i == 1
        title('速度分量');
    end
    grid on;
    % 加速度
    subplot(3,2,2*i);
    acc_curve = zeros(num_plot_pts, 1);
    for k = 1:num_plot_pts
        CK = CurveDerivs1(u_plot(k), U, p, P_calc_T, 2, num_ctrl_pts - 1);
        acc_curve(k) = CK(3, i);
    end
    plot(u_plot, acc_curve, 'LineWidth', 2);
    ylabel(sprintf('d²C_%s/du²', char('x'+i-1)));
    if i == 1
        title('加速度分量'); 
    end
    grid on;
end
subplot(3,2,5); xlabel('参数 u');
subplot(3,2,6); xlabel('参数 u');

% 标准的 FindSpan 函数（NURBS Book Algorithm A2.1）
function span = FindSpan(n, p, u, U)
    % 特殊情况
    if u >= U(n+2)  % u在最后一个节点
        span = n;
        return;
    end
    if u <= U(p+1)  % u在第一个节点
        span = p;
        return;
    end
    
    % 二分查找
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

% 标准的 BasisFuns 函数（NURBS Book Algorithm A2.2）
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

% 标准的 DersBasisFuns 函数（NURBS Book Algorithm A2.3）
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
    
    % 计算导数
    a = zeros(2, p+1);
    for r = 0:p
        s1 = 0; s2 = 1;
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
            j = s1; s1 = s2; s2 = j;
        end
    end
    
    % 乘以正确的因子
    r = p;
    for k = 1:n
        for j = 0:p
            ders(k+1, j+1) = ders(k+1, j+1) * r;
        end
        r = r * (p - k);
    end
end

% 曲线求值函数
function C = CurvePoint(n, p, U, P, u)
    span = FindSpan(n, p, u, U);
    N = BasisFuns(span, u, p, U);
    
    C = zeros(size(P, 1), 1);
    for i = 0:p
        C = C + N(i+1) * P(:, span-p+i+1);
    end
end

% 曲线导数函数
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

% 非均匀有理b样条
function C = NURBSCurvePoint(n, p, U, P, w, u)
    span = FindSpan(n, p, u, U);
    N = BasisFuns(span, u, p, U);
    
    C = zeros(size(P, 1), 1);
    denominator = 0.0;
    
    for i = 0:p
        idx = span - p + i + 1;
        C = C + N(i+1) * w(idx) * P(:, idx);
        denominator = denominator + N(i+1) * w(idx);
    end
    
    C = C / denominator;
end