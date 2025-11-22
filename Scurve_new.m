q0 = 0; q1 = 10; v0 = 1; v1 = 0;
vmax = 5; amax = 10; jmax = 30; jmin = -jmax;

if (vmax - v0) * jmax < amax^2
    Tj1 = sqrt((vmax - v0) / jmax);
    Ta = 2 * Tj1;
else
    Tj1 = amax / jmax;
    Ta = Tj1 + (vmax - v0) / amax;
end

if (vmax - v1) * jmax < amax^2
    Tj2 = sqrt((vmax - v1) / jmax);
    Td = 2 * Tj2;
else
    Tj2 = amax / jmax;
    Td = Tj2 + (vmax - v1) / amax;
end

Tv = (q1 - q0) / vmax - 1/2 * Ta * (1 + v0 / vmax) - 1/2 * Td * (1 + v1 / vmax);
if Tv > 0
    vlim = vmax;
    alim_a = amax;
    alim_d = amax;
elseif Tv < 0
    Tv = 0;
    Tj1 = amax / jmax;
    Tj2 = Tj1;
    delta = amax^4 / jmax^2 + 2 * (v0^2 + v1^2) + amax * (4 * (q1-q0) - 2 * (amax/jmax) * (v0 +v1));
    Ta = (amax^2/jmax - 2*v0 + sqrt(delta)) / (2 * amax);
    Td = (amax^2/jmax - 2*v1 + sqrt(delta)) / (2 * amax);
    if Ta < 0
        Ta = 0;
        Tj1 = 0;
        Td = 2 * (q1 - q0) / (v1 + v0);
        Tj2 = (jmax * (q1-q0) - sqrt(jmax*(jmax*(q1-q0)^2 + (v1+v0)^2 * (v1-v0)))) / (jmax * (v1 + v0));
    end
    if Td < 0
        Td = 0;
        Tj2 = 0;
        Ta = 2 * (q1 - q0) / (v1 + v0);
        Tj1 = (jmax * (q1-q0) - sqrt(jmax*(jmax*(q1-q0)^2 - (v1+v0)^2 * (v1-v0)))) / (jmax * (v1 + v0));
    end
    alim_a = jmax * Tj1;
    alim_d = -jmax * Tj2;
    vlim = v0 + (Ta - Tj1) * alim_a;
end

T = Ta + Td + Tv;
dt = 0.05;
t_segment = 0:dt:T;
n = length(t_segment);
q_all = zeros(1, n);
q_d_all = zeros(1, n);
q_dd_all = zeros(1, n);
q_ddd_all = zeros(1, n);
for i = 1:n
    t = t_segment(i);
    if t >= 0 && t <= Tj1
        q_all(i) = q0 + v0 * t + jmax * (t^3 / 6);
        q_d_all(i) = v0 + jmax * (t^2 / 2);
        q_dd_all(i) = jmax * t;
        q_ddd_all(i) = jmax;
    elseif t >= Tj1 && t <= Ta - Tj1
        q_all(i) = q0 + v0 * t + alim_a / 6 * (3*t^2 - 3*Tj1*t + Tj1^2);
        q_d_all(i) = v0 + alim_a * (t - Tj1 / 2);
        q_dd_all(i) = jmax * Tj1;
        q_ddd_all(i) = 0;
    elseif t >= Ta - Tj1 && t <= Ta
        q_all(i) = q0 + (vlim + v0) * (Ta / 2) - vlim * (Ta-t) - jmin * (Ta-t)^3 / 6;
        q_d_all(i) = vlim + jmin * (Ta-t)^2 / 2;
        q_dd_all(i) = -jmin * (Ta-t);
        q_ddd_all(i) = jmin;
    elseif t >= Ta && t <= Ta + Tv
        q_all(i) = q0 + (vlim+v0) * Ta /2 + vlim * (t-Ta);
        q_d_all(i) = vlim;
        q_dd_all(i) = 0;
        q_ddd_all(i) = 0;
    elseif t >= T - Td && t <= T - Td + Tj2
        q_all(i) = q1 - (vlim + v1) * Td / 2 + vlim * (t-T+Td) - jmax * (t-T+Td)^3 / 6;
        q_d_all(i) = vlim - jmax * (t-T+Td)^2 / 2;
        q_dd_all(i) = -jmax * (t-T+Td);
        q_ddd_all(i) = jmin;
    elseif t >= T - Td + Tj2 && t <= T - Tj2
        q_all(i) = q1 - (vlim+v1) * Td /2 + vlim * (t-T+Td) - (alim_d / 6) * (3*(t-T+Td)^2 - 3*Tj2*(t-T+Td) + Tj2^2);
        q_d_all(i) = vlim - alim_d * (t-T+Td-Tj2/2);
        q_dd_all(i) = -jmax * Tj2;
        q_ddd_all(i) = 0;
    elseif t >= T - Tj2 && t < T
        q_all(i) = q1 - v1 * (T - t) - jmax * (T-t)^3 / 6;
        q_d_all(i) = v1 + jmax * (T-t)^2 / 2;
        q_dd_all(i) = -jmax * (T-t);
        q_ddd_all(i) = jmax;
    end
end

figure;
subplot(4, 1, 1);
plot(t_segment, q_all, 'b-', 'LineWidth', 2);
grid on;
title('双S速度轨迹规划 - 位置曲线');
ylabel('位置 q(t)');
xlabel('时间 (s)');

subplot(4, 1, 2);
plot(t_segment, q_d_all, 'r-', 'LineWidth', 2);
grid on;
title('速度曲线');
ylabel('速度 v(t)');
xlabel('时间 (s)');

subplot(4, 1, 3);
plot(t_segment, q_dd_all, 'g-', 'LineWidth', 2);
grid on;
title('加速度曲线');
ylabel('加速度 a(t)');
xlabel('时间 (s)');

subplot(4, 1, 4);
plot(t_segment, q_ddd_all, 'm-', 'LineWidth', 2);
grid on;
title('加加速度曲线');
ylabel('加加速度 j(t)');
xlabel('时间 (s)');


