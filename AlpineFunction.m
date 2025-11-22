% 基准函数6
function f_6 = AlpineFunction(x)
    % f6
    d = length(x);
    sum = 0.0;
    for i = 1:d
        sum = sum + abs(x(i) * sin(x(i)) + 0.1 * x(i));
    end
    f_6 = sum;
end

