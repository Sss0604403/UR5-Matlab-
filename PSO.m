function [gBest_score, gBest_position, convergence_curve] = PSO(popsize, max_iter, lb, ub, dim, fobj)
    w = 0.9;
    c1 = 2;
    c2 = 2;
    % 初始化
    positions = zeros(popsize, dim);
    for i = 1:popsize
        positions(i, :) = lb + (ub - lb) .* rand(1, dim);
    end
    velocities = zeros(popsize, dim);

    fitness = zeros(popsize, 1);
    for i = 1:popsize
        fitness(i) = fobj(positions(i, :));
    end

    pBest_position = positions;  % 粒子局部最优
    pBest_scores = fitness;
    [gBest_score, gBest_idx] = min(fitness);
    gBest_position = positions(gBest_idx, :);

    convergence_curve = zeros(1, max_iter);
    % 主循环
    for t = 1:max_iter
        for i = 1:popsize
            r1 = rand(1, dim);
            r2 = rand(1, dim);

            velocities(i, :) = w * velocities(i, :) + c1 * r1 .* (pBest_position(i, :) - positions(i, :)) ...
                                                    + c2 * r2 .* (gBest_position - positions(i, :));
            positions(i, :) = positions(i, :) + velocities(i, :);
            positions(i, :) = max(positions(i, :), lb);
            positions(i, :) = min(positions(i, :), ub);

            new_fitness = fobj(positions(i, :));
            if new_fitness < pBest_scores(i)
                pBest_scores(i) = new_fitness;
                pBest_position(i, :) = positions(i, :);
            end

            if new_fitness < gBest_score
                gBest_score = new_fitness;
                gBest_position = positions(i, :);
            end
        end
        convergence_curve(t) = gBest_score;

        fprintf('迭代次数: %d, 最优适应度: %f\n', t, gBest_score);
    end
end

