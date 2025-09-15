function T = forwardmotion(theta)
    %DH参数
    a = [0,-0.425,-0.39225,0,0,0];
    d = [0.089159,0,0,0.10915,0.09465,0.08230];
    alpha = [pi/2,0,0,pi/2,-pi/2,0];

    %每隔一个自由度的变换矩阵
    T01 = T_para(theta(1),d(1),a(1),alpha(1));
    T12 = T_para(theta(2),d(2),a(2),alpha(2));
    T23 = T_para(theta(3),d(3),a(3),alpha(3));
    T34 = T_para(theta(4),d(4),a(4),alpha(4));
    T45 = T_para(theta(5),d(5),a(5),alpha(5));
    T56 = T_para(theta(6),d(6),a(6),alpha(6));
    
    T = T01 * T12 * T23 * T34 * T45 * T56;%末端位置
end

%计算i-1关节与i关节之间的转换关系
function T = T_para(theta,d,a,alpha)
    Ttheta = [cos(theta),-sin(theta),0,0;
        sin(theta),cos(theta),0,0;
        0,0,1,0;
        0,0,0,1];
    Td = [1,0,0,0;
          0,1,0,0;
          0,0,1,d;
          0,0,0,1];
    Talpha = [1,0,0,0;
        0,cos(alpha),-sin(alpha),0;
        0,sin(alpha),cos(alpha),0;
        0,0,0,1];
    Ta = [1,0,0,a;
          0,1,0,0;
          0,0,1,0;
          0,0,0,1];

    T = Ttheta * Td * Talpha * Ta;
end

