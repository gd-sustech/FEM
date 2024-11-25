function sliding_motion_rk4()
    % 参数设置
    m = 0.7;     % 质量 (kg)
    k = 0.01;    % 空气阻力系数 (kg/m)
    g = 9.81;    % 重力加速度 (m/s^2)
    mu = 0.5;    % 摩擦系数

    % 初始条件
    x0 = 0;      % 初始位置 (m)
    v0 = 3;      % 初始速度 (m/s)
    init_conditions = [x0; v0];

    % 时间设置
    dt = 0.001;  % 减小时间步长 (s)
    t_end = 10;  % 结束时间 (s)
    t = 0:dt:t_end; % 时间数组
    n = length(t); % 时间步数

    % 初始化数组
    x = zeros(1, n);
    v = zeros(1, n);
    x(1) = x0;   % 设置初始位置
    v(1) = v0;   % 设置初始速度

    % 四阶龙格-库塔法迭代
    for i = 1:n-1
        k1 = odefun(t(i), [x(i); v(i)], m, k, g, mu);
        k2 = odefun(t(i) + dt/2, [x(i) + dt/2 * k1(1); v(i) + dt/2 * k1(2)], m, k, g, mu);
        k3 = odefun(t(i) + dt/2, [x(i) + dt/2 * k2(1); v(i) + dt/2 * k2(2)], m, k, g, mu);
        k4 = odefun(t(i) + dt, [x(i) + dt * k3(1); v(i) + dt * k3(2)], m, k, g, mu);
        
        % 更新位置和速度
        x(i+1) = x(i) + dt/6 * (k1(1) + 2*k2(1) + 2*k3(1) + k4(1));
        v(i+1) = v(i) + dt/6 * (k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
        
        % 检查速度是否过小，防止数值不稳定
        if abs(v(i+1)) < 1e-5
            v(i+1) = 0; % 如果速度太小则设为0
        end
    end

    % 绘图
    figure;
    plot(t, x, 'LineWidth', 2);
    xlabel('时间 (s)');
    ylabel('位置 x (m)');
    title('物体滑行运动 (四阶龙格-库塔法)');
    grid on;
end

function dydt = odefun(t, y, m, k, g, mu)
    x = y(1);
    v = y(2);
    
    % 计算加速度
    a = -(k/m) * v^2 - g * mu;

    % 返回微分方程的右边
    dydt = [v; a];
end
