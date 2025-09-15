theta0 = [1,1,1,1,1,1];
T1 = forwardmotion(theta0);
disp("正运动学末端输出矩阵：");
disp(T1);  

theta_all = inversemotion(T1);
disp("逆运动学反求关节角度: ");
disp(theta_all);

theta_1 = theta_all(1,:);
T2 = forwardmotion(theta_1);
disp("验证正逆运动学算法的有效性");
disp(T2);

