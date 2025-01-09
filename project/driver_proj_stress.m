
% 计算单元应力
stress = zeros(n_elements, 3); % 每个单元的平均应力 [σx, σy, τxy]

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

    % 初始化单元应力
    stress_ele = zeros(3, 1);

    % 计算单元内的平均应力
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

        % 计算高斯点应力
        stress_gauss = C * strain_gauss;
        stress_ele = stress_ele + stress_gauss * weight(ll);
    end

    % 平均应力
    stress(ee, :) = stress_ele / sum(weight);
end

% 计算等效应力 (Von Mises 应力)
stress_eq = sqrt(stress(:, 1).^2 + stress(:, 2).^2 - stress(:, 1) .* stress(:, 2) + 3 * stress(:, 3).^2);

% 绘制等效应力分布到网格
figure;
patch('Faces', IEN, 'Vertices', [x_coor, y_coor], 'FaceVertexCData', stress_eq, ...
    'FaceColor', 'flat', 'EdgeColor', 'k');
colorbar;
title('四分之一带孔平板的应力分布 (Von Mises 应力)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;
