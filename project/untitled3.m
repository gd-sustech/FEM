% 计算单元应变
n_elements = size(IEN, 1); % 单元数量
strain = zeros(n_elements, 3); % 每个单元的平均应变 [εx, εy, γxy]

for ee = 1:n_elements
    % 获取单元节点坐标和位移
    nodes = IEN(ee, :);
    x_ele = x_coor(nodes);
    y_ele = y_coor(nodes);
    d_ele = [u(nodes), v(nodes)]'; % 位移向量，维度 [6x1]
    d_ele = d_ele(:); % 转为列向量

    % 高斯积分点
    n_int = 2; % 每方向积分点数
    [xi, eta, weight] = Gauss2D(n_int, n_int);

    % 初始化单元应变
    strain_ele = zeros(3, 1);

    % 计算单元内的平均应变
    for ll = 1:n_int^2
        % 形函数及其导数
        [N, dN_dxi, dN_deta] = Tri3ShapeFunctions(xi(ll), eta(ll));
        J = [dN_dxi * x_ele, dN_dxi * y_ele; dN_deta * x_ele, dN_deta * y_ele];
        detJ = det(J);
        invJ = inv(J);
        dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
        dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;

        % 构建应变-位移矩阵 B
        B = zeros(3, 6);
        for i = 1:3
            B(:, (i-1)*2+1:(i-1)*2+2) = [dN_dx(i), 0; 0, dN_dy(i); dN_dy(i), dN_dx(i)];
        end

        % 计算高斯点应变
        strain_gauss = B * d_ele;
        strain_ele = strain_ele + strain_gauss * weight(ll);
    end

    % 平均应变
    strain(ee, :) = strain_ele / sum(weight);
end

% 可视化应变
% 选择应变分量 εx
strain_x = strain(:, 1);

% 计算单元中心点
x_centroid = mean(x_coor(IEN), 2);
y_centroid = mean(y_coor(IEN), 2);

% 绘制应变分布图
figure;
scatter(x_centroid, y_centroid, 30, strain_x, 'filled');
colorbar;
title('四分之一带孔平板的应变分布 (εx)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;

% 可选：绘制等效应变分布
strain_eq = sqrt(strain(:, 1).^2 + strain(:, 2).^2 + 0.5 * strain(:, 3).^2);
figure;
scatter(x_centroid, y_centroid, 30, strain_eq, 'filled');
colorbar;
title('四分之一带孔平板的等效应变分布');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;
