clear;
clc;
close all;

% 1. 参数设置
Num_runs = 30;          % 独立运行次数
popsize = 30;           % 种群数量
max_iter = 1000;        % 最大迭代次数

% 函数1的设置
% fobj = @SphereFunction;     
% dim = 2;               
% lb = -100 * ones(1, dim);               
% ub = 100 * ones(1, dim);                
% 函数2的设置
% fobj = @DixonPrice;     
% dim = 2;               
% lb = -10 * ones(1, dim);              
% ub = 10 * ones(1, dim);                
% 函数3的设置
% fobj = @MatyasFunction;     
% dim = 2;               
% lb = -10 * ones(1, dim);               
% ub = 10 * ones(1, dim);                
% 函数4的设置
% fobj = @LeonFunction;     
% dim = 2;               
% lb = -1.2 * ones(1, dim);               
% ub = 1.2 * ones(1, dim);                
% 函数5的设置
% fobj = @BealeFunction;     
% dim = 2;               
% lb = -4.5 * ones(1, dim);               
% ub = 4.5 * ones(1, dim);                
% 函数6的设置
fobj = @AlpineFunction;     
dim = 2;               
lb = -10 * ones(1, dim);               
ub = 10 * ones(1, dim);                
% 函数7的设置
% fobj = @EasonFunction;     
% dim = 2;               
% lb = -100 * ones(1, dim);               
% ub = 100 * ones(1, dim);                
% 函数8的设置
% fobj = @ShubertFunction;     
% dim = 2;               
% lb = -10 * ones(1, dim);               
% ub = 10 * ones(1, dim);                
% 函数9的设置
% fobj = @Cross_in_TrayFunction;     
% dim = 2;               
% lb = -10 * ones(1, dim);               
% ub = 10 * ones(1, dim);                
% 函数10的设置
% fobj = @ChenBirdFunction;     
% dim = 2;               
% lb = -500 * ones(1, dim);               
% ub = 500 * ones(1, dim);                

% 准备存储所有运行结果
convergence_curves_SSA = zeros(Num_runs, max_iter);
best_fitnesses_SSA = zeros(Num_runs, 1);

convergence_curves_ISSA = zeros(Num_runs, max_iter);
best_fitnesses_ISSA = zeros(Num_runs, 1);

convergence_curves_PSO = zeros(Num_runs, max_iter);
best_fitnesses_PSO = zeros(Num_runs, 1);

convergence_curves_GWO = zeros(Num_runs, max_iter);
best_fitnesses_GWO = zeros(Num_runs, 1);

convergence_curves_WOA = zeros(Num_runs, max_iter);
best_fitnesses_WOA = zeros(Num_runs, 1);


% 2. 循环运行算法以收集数据 
fprintf('开始对 %s 函数进行测试，共运行 %d 次...', func2str(fobj), Num_runs);

for i = 1:Num_runs
    fprintf('第 %d 次运行:', i);
    
    % 运行 SSA 
    fprintf('  正在运行标准 SSA 算法...');
    [best_fitness_SSA, ~, conv_curve_SSA] = SSA(popsize, max_iter, lb, ub, dim, fobj);
    convergence_curves_SSA(i, :) = conv_curve_SSA;
    best_fitnesses_SSA(i) = best_fitness_SSA;
    fprintf(' 完成. 最优值: %e', best_fitness_SSA);

    % 运行 ISSA 
    fprintf('  正在运行改进的 ISSA 算法...');
    [best_fitness_ISSA, ~, conv_curve_ISSA] = ISSA(popsize, max_iter, lb, ub, dim, fobj);
    convergence_curves_ISSA(i, :) = conv_curve_ISSA;
    best_fitnesses_ISSA(i) = best_fitness_ISSA;
    fprintf(' 完成. 最优值: %e', best_fitness_ISSA);

    % 运行 PSO 
    fprintf('  正在运行改进的 PSO 算法...');
    [best_fitness_PSO, ~, conv_curve_PSO] = PSO(popsize, max_iter, lb, ub, dim, fobj);
    convergence_curves_PSO(i, :) = conv_curve_PSO;
    best_fitnesses_PSO(i) = best_fitness_PSO;
    fprintf(' 完成. 最优值: %e', best_fitness_PSO);
    
    % 运行 GWO 
    fprintf('  正在运行改进的 GWO 算法...');
    [best_fitness_GWO, ~, conv_curve_GWO] = GWO(popsize, max_iter, lb, ub, dim, fobj);
    convergence_curves_GWO(i, :) = conv_curve_GWO;
    best_fitnesses_GWO(i) = best_fitness_GWO;
    fprintf(' 完成. 最优值: %e', best_fitness_GWO);

    % 运行 WOA 
    fprintf('  正在运行改进的 WOA 算法...');
    [best_fitness_WOA, ~, conv_curve_WOA] = WOA(popsize, max_iter, lb, ub, dim, fobj);
    convergence_curves_WOA(i, :) = conv_curve_WOA;
    best_fitnesses_WOA(i) = best_fitness_WOA;
    fprintf(' 完成. 最优值: %e', best_fitness_WOA);
end


% 3. 计算并显示统计结果 
fprintf('----------------- 统计结果 -----------------');
fprintf('算法	Mean		Std');
fprintf('\n');

mean_SSA = mean(best_fitnesses_SSA);
std_SSA = std(best_fitnesses_SSA);
fprintf('SSA	%e	%e;\n', mean_SSA, std_SSA);

mean_ISSA = mean(best_fitnesses_ISSA);
std_ISSA = std(best_fitnesses_ISSA);
fprintf('ISSA	%e	%e;\n', mean_ISSA, std_ISSA);

mean_PSO = mean(best_fitnesses_PSO);
std_PSO = std(best_fitnesses_PSO);
fprintf('PSO	%e	%e;\n', mean_PSO, std_PSO);

mean_GWO = mean(best_fitnesses_GWO);
std_GWO = std(best_fitnesses_GWO);
fprintf('GWO	%e	%e;\n', mean_GWO, std_GWO);

mean_WOA = mean(best_fitnesses_WOA);
std_WOA = std(best_fitnesses_WOA);
fprintf('WOA	%e	%e;\n', mean_WOA, std_WOA);

fprintf('----------------------------------------------');


% 4. 绘制收敛曲线 (图3 的内容)
figure;
% 计算平均收敛曲线
mean_convergence_SSA = mean(convergence_curves_SSA, 1);
mean_convergence_ISSA = mean(convergence_curves_ISSA, 1);
mean_convergence_PSO = mean(convergence_curves_PSO, 1);
mean_convergence_GWO = mean(convergence_curves_GWO, 1);
mean_convergence_WOA = mean(convergence_curves_WOA, 1);

% 使用 semilogy 对数坐标轴
semilogy(1:max_iter, mean_convergence_SSA, 'b--', 'LineWidth', 1.5);
hold on;
semilogy(1:max_iter, mean_convergence_ISSA, 'r-', 'LineWidth', 1.5);
semilogy(1:max_iter, mean_convergence_PSO, 'g-.', 'LineWidth', 1.5);
semilogy(1:max_iter, mean_convergence_GWO, 'm:', 'LineWidth', 1.5);
semilogy(1:max_iter, mean_convergence_WOA, 'k:', 'LineWidth', 1.5);

title(['收敛曲线对比 - ' func2str(fobj)]);
xlabel('迭代次数');
ylabel('平均最优适应度值 (对数刻度)');
legend('SSA', 'ISSA', 'PSO', 'GWO', 'WOA'); % 添加其他算法名 'SSA', 'ISSA', 'PSO'
grid on;
hold off;

fprintf('统计结果和收敛曲线已生成。');

% 5. 绘制二维参数空间图 (可选，仅在 dim=2 时有意义)
if dim == 2
    fprintf('正在为二维函数绘制参数空间图...');
    
    % 1. 创建一个二维网格
    % 从每个维度的下界到上界，均匀取100个点
    x1 = linspace(lb(1), ub(1), 100);
    x2 = linspace(lb(2), ub(2), 100);
    [X1, X2] = meshgrid(x1, x2);
    
    % 2. 计算网格上每个点的函数值
    Z = zeros(size(X1));
    for i = 1:size(X1, 1)
        for j = 1:size(X1, 2)
            % 将网格点 (X1(i,j), X2(i,j)) 作为输入传给目标函数
            Z(i, j) = fobj([X1(i, j), X2(i, j)]);
        end
    end
    
    % 3. 创建一个新的图形窗口并绘图
    figure;
    
    % 使用 surf 绘制三维曲面图，并设置漂亮的颜色
    surf(X1, X2, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    colormap('jet'); % 使用 'jet' 或 'parula' 等颜色映射
    % colorbar;       % 显示颜色条
    hold on;
    
    % 使用 contour 在底部绘制等高线图
    contour(X1, X2, Z, 20, 'LineWidth', 1, 'HandleVisibility', 'off');
    
    title(['参数空间 - ' func2str(fobj)]);
    xlabel('x_1');
    ylabel('x_2');
    zlabel(['f(' func2str(fobj) ')']);
    view(3); % 设置为三维视角
    grid on;
    axis tight; % 紧凑坐标轴
    
    fprintf('参数空间图绘制完成！');
else
    fprintf('注意: 当前维度为 %d，无法绘制二维参数空间图。', dim);
    fprintf('如需绘图，请在脚本开头设置 dim = 2 并重新运行。');
end

