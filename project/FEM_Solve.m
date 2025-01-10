function [u, v, stress, stress_eq] = FEM_Solve(x_coor, y_coor, IEN, E, nu, sigma_inf, r_hole, xc, yc)
    % 平面应力本构矩阵
    C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    
    % 自由度编号
    n_nodes = size(x_coor, 1);
    ID = zeros(n_nodes, 2);
    counter = 0;
    for i = 1:n_nodes
        counter = counter + 1;
        ID(i, 1) = counter;
        counter = counter + 1;
        ID(i, 2) = counter;
    end
    n_eq = max(ID(:));
    
    % LM 数组
    LM = zeros(size(IEN, 1), 6);
    for ee = 1:size(IEN, 1)
        for aa = 1:3
            node = IEN(ee, aa);
            LM(ee, (aa - 1) * 2 + 1) = ID(node, 1);
            LM(ee, (aa - 1) * 2 + 2) = ID(node, 2);
        end
    end
    
    % 全局刚度矩阵和力向量
    K = spalloc(n_eq, n_eq, 36 * size(IEN, 1));
    F = zeros(n_eq, 1);

    % 高斯积分点
    n_int = 3;
    [xi, eta, weight] = Gauss2D(n_int, n_int);

    % 单元刚度矩阵组装
    for ee = 1:size(IEN, 1)
        x_ele = x_coor(IEN(ee, :));
        y_ele = y_coor(IEN(ee, :));
        k_ele = zeros(6, 6);
        for ll = 1:n_int^2
            [N, dN_dxi, dN_deta] = Tri3ShapeFunctions(xi(ll), eta(ll));
            J = [dN_dxi * x_ele, dN_dxi * y_ele; dN_deta * x_ele, dN_deta * y_ele];
            detJ = det(J);
            invJ = inv(J);
            dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
            dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;
            B = zeros(3, 6);
            for i = 1:3
                B(:, (i-1)*2+1:(i-1)*2+2) = [dN_dx(i), 0; 0, dN_dy(i); dN_dy(i), dN_dx(i)];
            end
            k_ele = k_ele + B' * C * B * detJ * weight(ll);
        end
        for i = 1:6
            for j = 1:6
                if LM(ee, i) > 0 && LM(ee, j) > 0
                    K(LM(ee, i), LM(ee, j)) = K(LM(ee, i), LM(ee, j)) + k_ele(i, j);
                end
            end
        end
    end

    % 边界条件
    fixed_u = find(abs(x_coor - (-1)) < 1e-5);
    for i = 1:length(fixed_u)
        dof = ID(fixed_u(i), 1);
        K(dof, :) = 0;
        K(dof, dof) = 1;
        F(dof) = 0;
    end
    fixed_v = find(abs(y_coor - (-1)) < 1e-5);
    for i = 1:length(fixed_v)
        dof = ID(fixed_v(i), 2);
        K(dof, :) = 0;
        K(dof, dof) = 1;
        F(dof) = 0;
    end
    right_edge = find(abs(x_coor - 1) < 1e-5);
    force_per_node = sigma_inf;
    for i = 1:length(right_edge)
        dof = ID(right_edge(i), 1);
        F(dof) = F(dof) + force_per_node;
    end

    % 求解位移
    d = K \ F;
    u = d(1:2:end);
    v = d(2:2:end);

    % Von Mises 应力（从应力回收）
    stress = zeros(size(IEN, 1), 3);
    stress_eq = zeros(size(IEN, 1), 1);
    % (可以添加计算应力部分)
end
