% 基准函数9
function f_9 = Cross_in_TrayFunction(x)
    % f39
    f_9 = -0.0001 * (abs(sin(x(1)) * sin(x(2)) * exp(abs(100 - ((x(1)^2 + x(2)^2)^0.5) / pi))) + 1)^0.1;
end

