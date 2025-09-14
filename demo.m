theta0 = [1,1,1,1,1,1];
T = forwardmotion(theta0);
disp("正运动学末端输出矩阵：");
disp(T);  

theta = inversemotion(T);
disp("逆运动学反求关节角度: ");
disp(theta);

T = forwardmotion(theta);
disp("验证正逆运动学算法的有效性");
disp(T);

