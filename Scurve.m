clear; clc; close all;

a = [0, -0.42500, -0.39225, 0, 0, 0];
d = [0.089159, 0, 0, 0.10915, 0.09465, 0.08230];
alpha = [pi/2, 0, 0, pi/2, -pi/2, 0];

L1 = Link('d', d(1), 'a', a(1), 'alpha', alpha(1), 'standard');
L2 = Link('d', d(2), 'a', a(2), 'alpha', alpha(2), 'standard');
L3 = Link('d', d(3), 'a', a(3), 'alpha', alpha(3), 'standard');
L4 = Link('d', d(4), 'a', a(4), 'alpha', alpha(4), 'standard');
L5 = Link('d', d(5), 'a', a(5), 'alpha', alpha(5), 'standard');
L6 = Link('d', d(6), 'a', a(6), 'alpha', alpha(6), 'standard');

UR5 = SerialLink([L1, L2, L3, L4, L5, L6], 'name', 'UR5');

joint_limits = struct();
joint_limits.vmax = [180, 180, 180, 360, 360, 360] * pi/180;  
joint_limits.amax = [300, 300, 300, 600, 600, 600] * pi/180;  
joint_limits.jmax = [900, 900, 900, 1800, 1800, 1800] * pi/180;  

waypoints_deg = [
    0,    0,   0,   0,   0,   0;      
    30,   -45,  60,  0,   90,   0;    
    60,   -60,  90,  -45,  60,  45;   
    90,   -30,  45,  -90,  30,  90;   
    45,   -45,  30,  -45,  45,  45    
];
waypoints = waypoints_deg * pi/180;

[t_total, q_total, qd_total, qdd_total, qddd_total] = plan_joint_trajectory(waypoints, joint_limits);

visualize_trajectory(t_total, q_total, qd_total, qdd_total, qddd_total, waypoints);

animate_ur5_robot(UR5, t_total, q_total, waypoints, true);

plot_end_effector_trajectory(UR5, q_total);

function visualize_trajectory(t, q, qd, qdd, qddd, waypoints)
    figure('Position', [50, 50, 1400, 900], 'Name', 'UR5关节空间S型曲线轨迹');
    
    joint_names = {'关节1', '关节2', '关节3', '关节4', '关节5', '关节6'};
    colors = lines(6);
    n_joints = size(q, 2);
    
    q_deg = q * 180/pi;
    qd_deg = qd * 180/pi;
    qdd_deg = qdd * 180/pi;
    qddd_deg = qddd * 180/pi;
    waypoints_deg = waypoints * 180/pi;
    
    % 位置
    subplot(4, 1, 1);
    hold on;
    for j = 1:n_joints
        plot(t, q_deg(:, j), 'LineWidth', 2, 'Color', colors(j,:), 'DisplayName', joint_names{j});
    end
    grid on;
    xlabel('时间 (s)', 'FontSize', 11);
    ylabel('关节角度 (°)', 'FontSize', 11);
    title('关节位置曲线', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'FontSize', 9);
    
    % 速度
    subplot(4, 1, 2);
    hold on;
    for j = 1:n_joints
        plot(t, qd_deg(:, j), 'LineWidth', 2, 'Color', colors(j,:), 'DisplayName', joint_names{j});
    end
    grid on;
    xlabel('时间 (s)', 'FontSize', 11);
    ylabel('关节速度 (°/s)', 'FontSize', 11);
    title('关节速度曲线', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'FontSize', 9);
    
    % 加速度
    subplot(4, 1, 3);
    hold on;
    for j = 1:n_joints
        plot(t, qdd_deg(:, j), 'LineWidth', 2, 'Color', colors(j,:), 'DisplayName', joint_names{j});
    end
    grid on;
    xlabel('时间 (s)', 'FontSize', 11);
    ylabel('关节加速度 (°/s²)', 'FontSize', 11);
    title('关节加速度曲线', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'FontSize', 9);
    
    % 加加速度
    subplot(4, 1, 4);
    hold on;
    for j = 1:n_joints
        plot(t, qddd_deg(:, j), 'LineWidth', 2, 'Color', colors(j,:), 'DisplayName', joint_names{j});
    end
    grid on;
    xlabel('时间 (s)', 'FontSize', 11);
    ylabel('关节加加速度 (°/s³)', 'FontSize', 11);
    title('关节加加速度曲线 (Jerk)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'eastoutside', 'FontSize', 9);
end

function animate_ur5_robot(robot, t, q, waypoints, save_video)
    if nargin < 5
        save_video = false;
    end
    fprintf('\n开始UR5机器人动画演示...\n');
    
    fig = figure('Position', [100, 100, 1000, 800], 'Name', 'UR5机器人S型曲线运动');
    
    set(fig, 'Color', 'w');  % 白色背景
    
    % % ========== 创建视频写入对象 ==========
    % if save_video
    %     video_filename = sprintf('UR5_trajectory_%s.mp4', datestr(now, 'yyyymmdd_HHMMSS'));
    %     v = VideoWriter(video_filename, 'MPEG-4');
    %     v.FrameRate = 30;  % 帧率
    %     v.Quality = 95;    % 视频质量 (0-100)
    %     open(v);
    %     fprintf('正在录制视频: %s', video_filename);
    % end
    % 绘制初始姿态
    robot.plot(q(1,:), 'workspace', [-0.8 0.8 -0.8 0.8 0 1.0], 'trail', 'r-', 'fps', 30);
    % 标记路径点
    hold on;
    for i = 1:size(waypoints, 1)
        T = robot.fkine(waypoints(i,:));
        pos = transl(T);
        plot3(pos(1), pos(2), pos(3), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'LineWidth', 2);
        text(pos(1), pos(2), pos(3)+0.05, sprintf('P%d', i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'g');
    end
    title('UR5机器人S型曲线轨迹运动', 'FontSize', 14, 'FontWeight', 'bold');
     
    % 动画循环
    skip = 5;  % 每5帧显示一次
    for i = 1:skip:length(t)
        robot.animate(q(i,:));
        title(sprintf('UR5机器人运动 - 时间: %.2f s / %.2f s', t(i), t(end)), 'FontSize', 14, 'FontWeight', 'bold');
        drawnow;
        if save_video
            frame = getframe(fig);
            writeVideo(v, frame);
        end
        pause(0.01);
    end
    
    % 最终姿态
    robot.animate(q(end,:));
    title(sprintf('UR5机器人运动完成 - 总时间: %.2f s', t(end)), 'FontSize', 14, 'FontWeight', 'bold');
    if save_video
        % 在最后一帧停留1秒
        for i = 1:30
            frame = getframe(fig);
            writeVideo(v, frame);
        end
        close(v);
        fprintf('视频保存成功: %s', video_filename);
        fprintf('视频时长: %.2f 秒', length(t)/v.FrameRate);
    end
    fprintf('动画演示完成！\n');
end

function plot_end_effector_trajectory(robot, q)
    fprintf('\n计算末端执行器轨迹...\n');
    
    n = size(q, 1);
    ee_pos = zeros(n, 3);
    
    for i = 1:n
        T = robot.fkine(q(i,:));
        ee_pos(i,:) = transl(T);
    end
    
    figure('Position', [150, 150, 1000, 800], 'Name', '末端执行器轨迹');
    
    % 3D轨迹
    subplot(2, 2, [1, 3]);
    plot3(ee_pos(:,1), ee_pos(:,2), ee_pos(:,3), 'b-', 'LineWidth', 2);
    hold on;
    plot3(ee_pos(1,1), ee_pos(1,2), ee_pos(1,3), 'go', ...
          'MarkerSize', 15, 'MarkerFaceColor', 'g');
    plot3(ee_pos(end,1), ee_pos(end,2), ee_pos(end,3), 'ro', ...
          'MarkerSize', 15, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('X (m)', 'FontSize', 11);
    ylabel('Y (m)', 'FontSize', 11);
    zlabel('Z (m)', 'FontSize', 11);
    title('末端执行器三维轨迹', 'FontSize', 12, 'FontWeight', 'bold');
    legend('轨迹', '起点', '终点', 'Location', 'best');
    view(3);
    axis equal;
    
    % XY平面投影
    subplot(2, 2, 2);
    plot(ee_pos(:,1), ee_pos(:,2), 'b-', 'LineWidth', 2);
    hold on;
    plot(ee_pos(1,1), ee_pos(1,2), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(ee_pos(end,1), ee_pos(end,2), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('X (m)', 'FontSize', 11);
    ylabel('Y (m)', 'FontSize', 11);
    title('XY平面投影', 'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    
    % XZ平面投影
    subplot(2, 2, 4);
    plot(ee_pos(:,1), ee_pos(:,3), 'b-', 'LineWidth', 2);
    hold on;
    plot(ee_pos(1,1), ee_pos(1,3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(ee_pos(end,1), ee_pos(end,3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('X (m)', 'FontSize', 11);
    ylabel('Z (m)', 'FontSize', 11);
    title('XZ平面投影', 'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    
    fprintf('末端执行器轨迹绘制完成！\n');
    fprintf('轨迹长度: %.3f m\n', sum(sqrt(sum(diff(ee_pos).^2, 2))));
end