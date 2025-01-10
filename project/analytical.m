% 材料参数与远场应力
E = 10e9; % 杨氏模量 (Pa)
nu = 0.3; % 泊松比
sigma_inf = 1e4; % 远场拉伸应力 (Pa)
r_hole = 0.5; % 孔半径 (m)
xc = -1; % 孔中心 x 坐标
yc = -1; % 孔中心 y 坐标
x_coor = msh.POS(:, 1); % x 坐标
y_coor = msh.POS(:, 2); % y 坐标
IEN = msh.TRIANGLES(:, 1:3); % 三角形单元 (节点索引)
% 节点坐标 (四分之一平板的有限元网格)
x_nodes = x_coor; % x 坐标
y_nodes = y_coor; % y 坐标

% 初始化应力分布
sigma_x = zeros(size(x_nodes)); % x方向应力
sigma_y = zeros(size(y_nodes)); % y方向应力
tau_xy = zeros(size(x_nodes)); % 剪切应力
von_mises = zeros(size(x_nodes)); % 总应力 (Von Mises 应力)

% 逐节点计算解析解
for i = 1:length(x_nodes)
    % 将坐标转换为极坐标 (r, theta)
    x = x_nodes(i);
    y = y_nodes(i);
    r = sqrt((x - xc)^2 + (y - yc)^2); % 距离孔中心的径向距离
    theta = atan2(y - yc, x - xc); % 极角

    % 判断是否在四分之一平板的有效区域
    if r < r_hole || x < -1 || y < -1 || x > 1 || y > 1
        sigma_x(i) = NaN; % 忽略无效区域或孔内区域
        sigma_y(i) = NaN;
        tau_xy(i) = NaN;
        von_mises(i) = NaN;
        continue;
    end

    % 解析解公式 (极坐标下的应力分量)
    sigma_rr = sigma_inf * (1 - r_hole^2 / r^2) + ...
               sigma_inf * (1 + 3 * r_hole^4 / r^4) * cos(2 * theta);
    sigma_tt = sigma_inf * (1 + r_hole^2 / r^2) - ...
               sigma_inf * (1 + 3 * r_hole^4 / r^4) * cos(2 * theta);
    tau_rt = -sigma_inf * (1 - 3 * r_hole^4 / r^4) * sin(2 * theta);

    % 转换为笛卡尔坐标下的应力分量
    sigma_x(i) = sigma_rr * cos(theta)^2 + sigma_tt * sin(theta)^2 + 2 * tau_rt * sin(theta) * cos(theta);
    sigma_y(i) = sigma_rr * sin(theta)^2 + sigma_tt * cos(theta)^2 - 2 * tau_rt * sin(theta) * cos(theta);
    tau_xy(i) = (sigma_tt - sigma_rr) * sin(theta) * cos(theta) + tau_rt * (cos(theta)^2 - sin(theta)^2);

    % 计算 Von Mises 应力
    von_mises(i) = sqrt(sigma_x(i)^2 + sigma_y(i)^2 - sigma_x(i) * sigma_y(i) + 3 * tau_xy(i)^2);
end

num_elements = size(IEN, 1); % 单元总数
element_von_mises = zeros(num_elements, 1); % 每个单元的平均 Von Mises 应力

% 逐单元计算平均应力
for elem = 1:num_elements
    % 获取当前单元的节点索引
    nodes = IEN(elem, :);
    
    % 获取当前单元的 Von Mises 应力
    vm_stresses = von_mises(nodes);
    
    % 计算平均 Von Mises 应力
    element_von_mises(elem) = mean(vm_stresses, 'omitnan');
end

% 输出平均单元应力结果
% fprintf('每个单元的平均 Von Mises 应力:\n');
% disp(element_von_mises);

% 可视化单元平均应力分布
figure;
patch('Faces', IEN, 'Vertices', [x_coor, y_coor], ...
      'FaceVertexCData', element_von_mises, ... % 应力数据
      'FaceColor', 'flat', ...                 % 填充颜色
      'EdgeColor', 'k', ...                    % 网格线颜色 (黑色)
      'LineWidth', 0.5);                       % 网格线宽度
colorbar;
title('四分之一带孔平板的平均单元应力分布 (Von Mises 应力)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;
xlim([-1 1]);
ylim([-1 1]);
save('Stress_exact', 'von_mises');


