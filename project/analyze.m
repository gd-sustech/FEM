% 参数定义
sigma_inf = 1e7; % 远场拉伸应力 (Pa)
r_hole = 0.5;    % 孔半径 (m)
E = 1e10;        % 弹性模量 (Pa)
nu = 0.3;        % 泊松比

% 极坐标网格生成
theta = linspace(0, pi/2, 100); % 四分之一平板，角度从 0 到 pi/2
r = linspace(r_hole, 2, 100);   % 半径从孔边开始
[Theta, R] = meshgrid(theta, r);

% 解析解计算
sigma_rr = sigma_inf * (1 - (r_hole^2 ./ R.^2)) + ...
           sigma_inf * (1 + 3 * (r_hole^4 ./ R.^4)) .* cos(2 * Theta);

sigma_tt = sigma_inf * (1 + (r_hole^2 ./ R.^2)) - ...
           sigma_inf * (1 + 3 * (r_hole^4 ./ R.^4)) .* cos(2 * Theta);

tau_rt = -sigma_inf * (1 - 3 * (r_hole^4 ./ R.^4)) .* sin(2 * Theta);

% 等效应力 (Von Mises)
sigma_eq = sqrt(sigma_rr.^2 + sigma_tt.^2 - sigma_rr .* sigma_tt + 3 * tau_rt.^2);

% 可视化等效应力分布
figure;
polarplot = polaraxes;
hold on;
pcolor(R .* cos(Theta), R .* sin(Theta), sigma_eq); % 将极坐标转换为笛卡尔坐标
shading interp; % 平滑颜色
colorbar;
title('解析解的等效应力分布 (Von Mises 应力)');
xlabel('X 坐标 (m)');
ylabel('Y 坐标 (m)');
axis equal;
