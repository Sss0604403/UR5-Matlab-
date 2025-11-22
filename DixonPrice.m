% 基准函数2
function f_2 = DixonPrice(x)
    % (f48)
    d = length(x); % 获取维度
    term1 = (x(1) - 1)^2;
    sum_term = 0;
    for i = 2:d
        sum_term = sum_term + i * (2*x(i)^2 - x(i-1))^2;
    end
    
    f_2 = term1 + sum_term;
end

