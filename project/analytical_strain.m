%%%这段代码是用解析解来算应变u，v
% 材料参数与远场应力
E = 10e9; % 杨氏模量 (Pa)
nu = 0.3; % 泊松比
G = E / (2 * (1 + nu)); % 剪切模量
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

end

% 计算应变分量
epsilon_x = (1 / (E / (1 - nu^2))) * (sigma_x - nu * sigma_y);
epsilon_y = (1 / (E / (1 - nu^2))) * (sigma_y - nu * sigma_x);
gamma_xy = (1 / G) * tau_xy;

% 计算位移梯度
du_dx = epsilon_x;  % 位移在x方向的梯度
dv_dy = epsilon_y;  % 位移在y方向的梯度
du_dy = gamma_xy / 2; % 位移在y方向的梯度（剪切应变）
dv_dx = gamma_xy / 2; % 位移在x方向的梯度（剪切应变）

% 积分得到位移
u = zeros(size(x_nodes)); % x方向位移
v = zeros(size(y_nodes)); % y方向位移

% 积分计算位移
for i = 2:length(x_nodes)
    u(i) = u(i-1) + du_dx(i) * (x_nodes(i) - x_nodes(i-1)); % 积分计算u
    v(i) = v(i-1) + dv_dy(i) * (y_nodes(i) - y_nodes(i-1)); % 积分计算v
end
u_exact=u;
v_exact=v;

% 保存位移结果
save('Displacement_exact', 'u_exact', 'v_exact', 'x_coor', 'y_coor');

