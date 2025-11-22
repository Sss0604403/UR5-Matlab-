clear;clc;close all;

% 参数设置
popsize = 30;       % 麻雀种群数量
max_iter = 1000;    % 最大迭代次数
dim = 5;            % 变量维度
lb = [0, 0.01, 0.01, 0, 0];    % 变量下界 
ub = [200, 100, 100, 2, 2];      % 变量上界 
Kp = 1;
Ki = 1;
Kd = 1;
lambda = 0.5;
mu = 0.5;
% 被控函数
num = [1];
den = [4.32, 19.1801, 1];

% 运行 SSA 和 ISSA 算法 
fprintf('正在运行标准 SSA 算法...');
[best_fitness_SSA, best_position_SSA, convergence_curve_SSA] = SSA(popsize, max_iter, lb, ub, dim, @ObjectiveFunction_FOPID);

fprintf('正在运行改进的 ISSA 算法...');
[best_fitness_ISSA, best_position_ISSA, convergence_curve_ISSA] = ISSA(popsize, max_iter, lb, ub, dim, @ObjectiveFunction_FOPID);

disp('----------------- ISSA 优化结果 -----------------');
disp(['最优控制器参数 [Kp, Ki, Kd, lambda, mu]: ', num2str(best_position_SSA)]);
disp(['最优性能指标 Q: ', num2str(best_fitness_SSA)]);
disp('-------------------------------------------------');

disp('----------------- ISSA 优化结果 -----------------');
disp(['最优控制器参数 [Kp, Ki, Kd, lambda, mu]: ', num2str(best_position_ISSA)]);
disp(['最优性能指标 Q: ', num2str(best_fitness_ISSA)]);
disp('-------------------------------------------------');

figure; 
% 绘制 SSA 的收敛曲线，使用红色实线
plot(1:max_iter, convergence_curve_SSA, 'r-', 'LineWidth', 2);
hold on; 
% 绘制 ISSA 的收敛曲线，使用蓝色虚线
plot(1:max_iter, convergence_curve_ISSA, 'b--', 'LineWidth', 2);
title('SSA vs ISSA 算法收敛曲线对比');
xlabel('迭代次数');
ylabel('最优适应度值');
legend('SSA','ISSA'); 
grid on;
hold off; 

fprintf('正在使用最优参数绘制阶跃响应曲线...');
model_name = 'ISSA_FOPID_Controller'; 
sim_time = 4; 

Kp     = best_position_ISSA(1);
Ki     = best_position_ISSA(2);
Kd     = best_position_ISSA(3);
lambda = best_position_ISSA(4);
mu     = best_position_ISSA(5);
% 运行仿真，并将输出保存到 simOut 变量中
simOut_ISSA = sim(model_name, 'SimulationMode', 'normal', 'StopTime', num2str(sim_time), 'SrcWorkspace', 'current');
% 从仿真结果中提取需要的数据
t_ISSA = simOut_ISSA.e.Time;      % 时间向量
e = simOut_ISSA.e.Data;      % 误差 e = r - y
r = ones(size(t_ISSA));      % 创建一个与时间向量同样大小的参考输入信号 r=1
y_ISSA = r - e;              % 通过误差计算出系统输出 y = r - e

Kp     = best_position_SSA(1);
Ki     = best_position_SSA(2);
Kd     = best_position_SSA(3);
lambda = best_position_SSA(4);
mu     = best_position_SSA(5);
% 运行仿真，并将输出保存到 simOut 变量中
simOut_SSA = sim(model_name, 'SimulationMode', 'normal', 'StopTime', num2str(sim_time), 'SrcWorkspace', 'current');
% 从仿真结果中提取需要的数据
t_SSA = simOut_SSA.e.Time;      % 时间向量
e = simOut_SSA.e.Data;      % 误差 e = r - y
r = ones(size(t_SSA));      % 创建一个与时间向量同样大小的参考输入信号 r=1
y_SSA = r - e;              % 通过误差计算出系统输出 y = r - e

figure; 
plot(t_SSA, y_SSA, 'b-', 'LineWidth', 2);  % SSA
hold on; 
plot(t_ISSA, y_ISSA, 'r--', 'LineWidth', 2);  % ISSA
plot(t_ISSA, ones(size(t_ISSA)), 'k:', 'LineWidth', 1.5);  % 期望输入
hold off;
title('最优FOPID控制器下的系统阶跃响应', 'FontSize', 14);
xlabel('时间 / s', 'FontSize', 12);
ylabel('振幅', 'FontSize', 12); 
legend('SSA优化结果', 'ISSA优化结果', '期望输入', 'Location', 'southeast');
grid on; 
axis([0, sim_time, 0, 1.4]); 
fprintf('所有任务完成！');