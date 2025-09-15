function theta= inversemotion(T)

    %DH参数
    a = [0,-0.425,-0.39225,0,0,0];
    d = [0.089159,0,0,0.10915,0.09465,0.08230];
    alpha = [pi/2,0,0,pi/2,-pi/2,0];

    %将末端矩阵转为为noa形式
    nx = T(1,1); ny = T(2,1);nz=T(3,1);
    ox = T(1,2); oy = T(2,2);oz=T(3,2);
    ax = T(1,3); ay = T(2,3);az=T(3,3);
    px = T(1,4); py = T(2,4);pz=T(3,4);

    %求解关节角1
    m1 = d(6) * ay - py;
    n1 = ax * d(6) - px;
    theta1(1,1) = atan2(m1,n1) - atan2(d(4),sqrt(m1^2+n1^2-(d(4))^2));
    theta1(1,2) = atan2(m1,n1) - atan2(d(4),-sqrt(m1^2+n1^2-(d(4))^2));

    %求解关节角5
    theta5(1,1:2) = acos(ax * sin(theta1)-ay * cos(theta1));
    theta5(2,1:2) = -acos(ax * sin(theta1)-ay * cos(theta1));

    %求解关节角6
    m6 = nx * sin(theta1) - ny * cos(theta1);
    n6 = ox * sin(theta1) - oy * cos(theta1);
    theta6 = atan2(m6,n6) - atan2(sin(theta5),0);

    % %求解关节角3
    % m3 = d(5) * (sin(theta6) .* (nx * cos(theta1) + ny * sin(theta1)) + cos(theta6) .* (ox * cos(theta1) + oy * sin(theta1)) %少一个括号- d(6) * (ax * cos(theta1) + ay * sin(theta1)) + px * cos(theta1) + py * sin(theta1))%这多一个括号;
    % n3 = pz - d(1) - az * d(6) + d(5) * (oz * cos(theta6) + nz * sin(theta6));
    % theta3(1:2,:) = acos(m3.^2 + n3.^2 - a(2)^2 - a(3)^2) / (2 * a(2) * a(3));
    % theta3(3:4,:) = -acos(m3.^2 + n3.^2 - a(2)^2 - a(3)^2) / (2 * a(2) * a(3));

    %求解关节角3
    m3 = d(5) * (sin(theta6) .* (nx * cos(theta1) + ny * sin(theta1)) + cos(theta6) .* (ox * cos(theta1) + oy * sin(theta1))) - d(6) * (ax * cos(theta1) + ay * sin(theta1)) + px * cos(theta1) + py * sin(theta1);
    n3 = pz - d(1) - az * d(6) + d(5) * (oz * cos(theta6) + nz * sin(theta6));
    acos_arg = (m3.^2 + n3.^2 - (a(2))^2 - (a(3))^2) / (2 * a(2) * a(3));
    acos_arg = max(min(acos_arg, 1.0), -1.0);  % 强制限制在合法范围
    theta3(1:2,:) = acos(acos_arg);
    theta3(3:4,:) = -acos(acos_arg);


    %求解关节角2
    m2(1:2,:) = m3;
    m2(3:4,:) = m3;
    n2(1:2,:) = n3;
    n2(3:4,:) = n3;
    s2 = ((a(3) * cos(theta3) + a(2)) .* n2 - a(3) * sin(theta3) .* m2) ./ ((a(2))^2 + (a(3))^2 + 2 * a(2) * a(3) * cos(theta3));
    c2 = (m2 + a(3) * sin(theta3) .* s2) ./ (a(3) * cos(theta3) + a(2));
    % c2 = max(min(c2, 1.0), -1.0);  % 强制限制在[-1,1]
    % s2 = sqrt(1 - c2^2);            % 此时s2必为实数
    theta2 = atan2(s2,c2); 

    %整理关节角1 5 6 3 2
    theta(1:4,1) = theta1(1,1);
    theta(5:8,1) = theta1(1,2);
    theta(:,2) = [theta2(1,1),theta2(3,1),theta2(2,1),theta2(4,1),theta2(1,2),theta2(3,2),theta2(2,2),theta2(4,2)];
    theta(:,3) = [theta3(1,1),theta3(3,1),theta3(2,1),theta3(4,1),theta3(1,2),theta3(3,2),theta3(2,2),theta3(4,2)];
    theta(1:2,5) = theta5(1,1);
    theta(3:4,5) = theta5(2,1);
    theta(5:6,5) = theta5(1,2);
    theta(7:8,5) = theta5(2,2);
    theta(1:2,6) = theta6(1,1);
    theta(3:4,6) = theta6(2,1);
    theta(5:6,6) = theta6(1,2);
    theta(7:8,6) = theta6(2,2);

    %求解关节角4
    theta(:,4) = atan2(-sin(theta(:,6)) .* (nx * cos(theta(:,1)) + ny * sin(theta(:,1))) - cos(theta(:,6)) .* (ox * cos(theta(:,1)) + oy * sin(theta(:,1))),oz * cos(theta(:,6)) + nz * sin(theta(:,6))) - theta(:,2) - theta(:,3);
end

