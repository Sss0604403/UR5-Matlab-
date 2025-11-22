function [Alpha_score, Alpha_pos, convergence_curve] = GWO(popsize, max_iter, lb, ub, dim, fobj)
    Alpha_pos = zeros(1, dim);
    Alpha_score = inf;
    Beta_pos = zeros(1, dim);
    Beta_score = inf;
    Delta_pos = zeros(1, dim);
    Delta_score = inf;
    
    positions = zeros(popsize, dim);
    for i = 1:popsize
        positions(i, :) = lb + (ub - lb) .* rand(1, dim);
    end

    convergence_curve = zeros(1, max_iter);
    for t = 1:max_iter
        for i = 1:popsize
            positions(i, :) = max(positions(i, :), lb);
            positions(i, :) = min(positions(i, :), ub);

            fitness = fobj(positions(i, :));
            if fitness < Alpha_score
                % Alpha降级为Beta，Beta降级为Delta
                Delta_score = Beta_score;
                Delta_pos = Beta_pos;
                Beta_score = Alpha_score;
                Beta_pos = Alpha_pos;
                Alpha_score = fitness;
                Alpha_pos = positions(i, :);
            elseif fitness > Alpha_score && fitness < Beta_score
                % 当前狼介于Alpha和Beta之间
                Delta_score = Beta_score;
                Delta_pos = Beta_pos;
                Beta_score = fitness;
                Beta_pos = positions(i, :);
            elseif fitness > Beta_score && fitness < Delta_score
                % 当前狼介于Beta和Delta之间
                Delta_score = fitness;
                Delta_pos = positions(i, :);
            end
        end
        % 更新狼群位置
        a = 2 -t * (2 / max_iter);
        % 对每一只Omega狼
        for i = 1:popsize
            for j = 1:dim
                r1 = rand();
                r2 = rand();
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_alpha = abs(C1 * Alpha_pos(j) - positions(i, j));
                X1 = Alpha_pos(j) - A1 * D_alpha;

                r1 = rand();
                r2 = rand();
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_beta = abs(C2 * Beta_pos(j) - positions(i, j));
                X2 = Beta_pos(j) - A2 * D_beta;

                r1 = rand();
                r2 = rand();
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_delta = abs(C3 * Delta_pos(j) - positions(i, j));
                X3 = Delta_pos(j) - A3 * D_delta;
                
                positions(i, j) = (X1 + X2 + X3) / 3;
            end
        end
        convergence_curve(t) = Alpha_score;
        fprintf('迭代次数: %d, 最优适应度: %f\n', t, Alpha_score);
    end
end

