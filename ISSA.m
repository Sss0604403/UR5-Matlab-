function [best_fitness, best_position, convergence_curve] = ISSA(popsize, max_iter, lb, ub, dim, fobj)
    
    % 初始化参数
    p_percent = 0.2;                        % 发现者占种群的比例
    p_num = round(popsize * p_percent);     % 发现者的数量
    
    sd_percent = 0.1;                       % 警戒者占种群的比例
    sd_num = round(popsize * sd_percent);   % 警戒者的数量

    ST = 0.8;                               % 安全阈值

    % Chebyshev映射初始化
    Y = zeros(popsize, dim); k = 4;
    y_seed = 2 * rand(1, dim) - 1;
    Y(1, :) = y_seed;
    for i = 2:popsize
        Y(i, :) = cos(k * acos(Y(i-1, :)));
    end
    positions = zeros(popsize, dim);
    for d = 1:dim
        positions(:, d) = lb(d) + (1 + Y(:, d)) .* (ub(d) - lb(d)) / 2;
    end

    fitness = zeros(popsize, 1);
    for i = 1:popsize
        fitness(i) = fobj(positions(i, :));
    end
    % 初始化全局最优
    [best_fitness, best_idx] = min(fitness);
    best_position = positions(best_idx, :);
    % 用于记录每次迭代的最优值，以绘制收敛曲线
    convergence_curve = zeros(max_iter, 1);

    % SSA 主循环
    for t = 1:max_iter
        % 1. 对所有麻雀的适应度值进行排序，找到当前最优和最差
        [~, sort_idx] = sort(fitness);
        % 记录当前全局最优和最差的位置
        best_pos_t = positions(sort_idx(1), :); % 当前迭代最好的位置
        worst_pos_t = positions(sort_idx(end), :); % 当前迭代最差的位置

        % 2. 更新“发现者”的位置
        % 循环处理前 p_num 个（即适应度最好的）麻雀
        for i = 1:p_num
            idx = sort_idx(i); % 获取当前处理的发现者在原始种群中的索引
            R2 = rand();       % 预警值
            
            if R2 < ST % 环境安全，进行广泛搜索
                alpha = rand();
                positions(idx, :) = positions(idx, :) .* exp(-i / (alpha * max_iter));
            else       % 有危险，飞到其他地方
                positions(idx, :) = positions(idx, :) + randn() * ones(1, dim);
            end
        end

        % 3. 更新“加入者”的位置
        % 循环处理从 p_num+1 到 popsize 的麻雀
        for i = (p_num + 1):popsize
            idx = sort_idx(i); % 获取当前处理的加入者在原始种群中的索引
            if i > (popsize + p_num) / 2  % 对应的是适应度差的加入者，它们很饿，需要飞到别处觅食
                positions(idx, :) = randn() * exp((worst_pos_t - positions(idx, :)) / (i^2));
            else % 对应的是适应度好的加入者，它们会飞向当前最优发现者的周围
                 % 生成一个每行元素都为1或-1的矩阵
                A = ones(1, dim);
                A_plus = 1 ./ (A' * (A * A')^(-1)); % A+ = A'*(A*A')^-1
                positions(idx, :) = best_pos_t + abs(positions(idx, :) - best_pos_t) .* A_plus' .* randn();
            end
        end

        % 4. 更新“警戒者”的位置
        % 从整个种群中随机选择 sd_num 只作为警戒者
        rand_indices = randperm(popsize, sd_num); 
        for i = 1:sd_num
            idx = rand_indices(i); % 获取当前处理的警戒者在原始种群中的索引
            if fitness(idx) > best_fitness % 如果当前警戒者不是最优的麻雀
                % 向当前全局最优位置靠近
                positions(idx, :) = best_position + randn() * abs(positions(idx, :) - best_position);
            else % 如果当前警戒者本身就是最优的麻雀
                % 为了避免陷入局部最优，它会向邻近的麻雀移动
                K = 2 * rand() - 1;
                positions(idx, :) = positions(idx, :) + K * (abs(positions(idx, :) - worst_pos_t) ./ (fitness(idx) - fitness(sort_idx(end)) + 1e-8));
            end
        end

        p = 0.2; 
        if rand < p
            % 对适应度差的麻雀自适应t分布扰动 
            start_index = floor(popsize / 2) + 1; 
            for i = start_index : popsize
                % 从自由度为 t 的 t 分布中生成随机数
                t_rand_vec = trnd(t, [1, dim]);                
                positions(i, :) = positions(i, :) + positions(i, :) .* t_rand_vec;
            end
        else
            % 萤火虫算法
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

        % 5. 边界处理和适应度更新
        for i = 1:popsize
            positions(i, :) = max(positions(i, :), lb);  % 检查下界
            positions(i, :) = min(positions(i, :), ub);  % 检查上界
            % 重新计算更新后位置的适应度
            fitness(i) = fobj(positions(i, :));
        end

        % 6. 更新全局最优解
        [current_best_fitness, current_best_idx] = min(fitness);
        if current_best_fitness < best_fitness
            best_fitness = current_best_fitness;
            best_position = positions(current_best_idx, :);
        end
        % 记录本次迭代的最优值
        convergence_curve(t) = best_fitness;
        
        % 在命令窗口显示迭代过程
        fprintf('迭代次数: %d, 最优适应度: %f\n', t, best_fitness);
    end
end


