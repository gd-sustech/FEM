% 参数设置
m = 1;          % 质量
k = 0.1;       % 空气阻力常数
mu = 0.1;      % 摩擦系数
g = 9.81;      % 重力加速度
v0 = 10;       % 初始速度

% 时间设置
tspan = [0, 10]; % 时间范围
initial_conditions = [v0]; % 初始速度

% 微分方程定义
odefun = @(t, v) - (k/m) * v^2 - (mu/m) * g;

% 求解速度
[t, v] = ode45(odefun, tspan, initial_conditions);

% 计算位移
x = cumtrapz(t, v); % 使用梯形法则积分得到位移

% 绘图
figure;
subplot(2,1,1);
plot(t, v);
xlabel('时间 (s)');
ylabel('速度 (m/s)');
title('速度随时间变化图');
grid on;

subplot(2,1,2);
plot(t, x);
xlabel('时间 (s)');
ylabel('位移 (m)');
title('位移随时间变化图');
grid on;
