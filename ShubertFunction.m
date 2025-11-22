% 基准函数8
function f_8 = ShubertFunction(x)
    % f132
    d = length(x);
    multi = 1.0;
    for i = 1:d
        sum = 0.0;
        for j = 1:5
            sum = sum + cos((j+1) * x(i) + j);
        end
        multi = multi * sum;
    end
    f_8 = multi;
end

