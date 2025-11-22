function [Leader_score, Leader_pos, convergence_curve] = WOA(popsize, max_iter, lb, ub, dim, fobj)
    Leader_pos = zeros(1, dim);
    Leader_score = inf;
    positions = zeros(popsize, dim);
    for i = 1:popsize
        positions(i, :) = lb + (ub - lb) .* rand(1, dim);
    end
    convergence_curve = zeros(1, max_iter);
    b = 1; % 螺旋形状常数

    for t = 1:max_iter
        for i = 1:popsize
            positions(i, :) = max(positions(i, :), lb);
            positions(i, :) = min(positions(i, :), ub);

            fitness = fobj(positions(i, :));
            if fitness < Leader_score  % 更新领导者
                Leader_score = fitness;
                Leader_pos = positions(i, :);
            end
        end
        % 更新鲸鱼位置
        a = 2 - t * (2 / max_iter);
        for i = 1:popsize
            r1 = rand();
            r2 = rand();
            A = 2 * a * r1 - a;
            C = 2 * r2;
            p = rand();
            l = (a - 1) * rand() + 1;

            if p < 0.5
                if abs(A) < 1 % 收缩包围或随机搜寻
                    D = abs(C * Leader_pos - positions(i, :));
                    positions(i, :) = Leader_pos - A * D;
                else
                    rand_whale_idx = floor(popsize * rand() + 1);
                    X_rand = positions(rand_whale_idx, :);
                    D = abs(C * X_rand - positions(i, :));
                    positions(i, :) = X_rand - A * D;
                end
            else
                D_prime = abs(Leader_pos - positions(i, :));
                positions(i, :) = D_prime .* exp(b * 1) .* cos(2 * pi* l) + Leader_pos;
            end
        end
        convergence_curve(t) = Leader_score;

        fprintf('迭代次数: %d, 最优适应度: %f\n', t, Leader_score);
    end
end

