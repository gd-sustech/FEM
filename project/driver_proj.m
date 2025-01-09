
% % 加载 Gmsh 网格
% load('mesh.m'); % 确保 mesh 文件在当前目录下

% 提取节点和单元信息
x_coor = msh.POS(:, 1); % x 坐标
y_coor = msh.POS(:, 2); % y 坐标
IEN = msh.TRIANGLES(:, 1:3); % 三角形单元 (节点索引)

% 材料属性
E = 10e9;  % 弹性模量 (Pa)
nu = 0.3;  % 泊松比

% 平面应力本构矩阵
C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

% 自由度编号
n_np = size(x_coor, 1); % 总节点数
ID = zeros(n_np, 2); % 每个节点两个自由度 (u, v)
counter = 0;
for i = 1:n_np
    counter = counter + 1;
    ID(i, 1) = counter;   % u 自由度
    counter = counter + 1;
    ID(i, 2) = counter;   % v 自由度
end
n_eq = max(ID(:)); % 总方程数

% LM 数组 (局部到全局自由度映射)
LM = zeros(size(IEN, 1), 6); % 3 节点 * 2 DOFs/节点
for ee = 1:size(IEN, 1)
    for aa = 1:3
        node = IEN(ee, aa);
        LM(ee, (aa-1)*2 + 1) = ID(node, 1);
        LM(ee, (aa-1)*2 + 2) = ID(node, 2);
    end
end

% 全局刚度矩阵和力向量
K = spalloc(n_eq, n_eq, 36 * size(IEN, 1));
F = zeros(n_eq, 1);

% 高斯积分点
n_int = 2; % 每个方向的积分点数
[xi, eta, weight] = Gauss2D(n_int, n_int);

% 单元刚度矩阵组装
for ee = 1:size(IEN, 1)
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));
    
    k_ele = zeros(6, 6); % 单元刚度矩阵
    for ll = 1:n_int^2
        % 映射到物理域
        [N, dN_dxi, dN_deta] = Tri3ShapeFunctions(xi(ll), eta(ll));
        J = [dN_dxi * x_ele, dN_dxi * y_ele; dN_deta * x_ele, dN_deta * y_ele];
        detJ = det(J);
        invJ = inv(J);
        dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
        dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;
        
        B = zeros(3, 6); % 应变-位移矩阵
        for i = 1:3
            B(:, (i-1)*2+1:(i-1)*2+2) = [dN_dx(i), 0; 0, dN_dy(i); dN_dy(i), dN_dx(i)];
        end
        
        % 计算单元刚度矩阵
        k_ele = k_ele + B' * C * B * detJ * weight(ll);
    end
    
    % 组装到全局刚度矩阵
    for i = 1:6
        for j = 1:6
            if LM(ee, i) > 0 && LM(ee, j) > 0
                K(LM(ee, i), LM(ee, j)) = K(LM(ee, i), LM(ee, j)) + k_ele(i, j);
            end
        end
    end
end

% 施加边界条件
% 左边界 u = 0
fixed_u = find(abs(x_coor - (-1)) < 1e-5);
for i = 1:length(fixed_u)
    dof = ID(fixed_u(i), 1);
    K(dof, :) = 0;
    K(dof, dof) = 1;
    F(dof) = 0;
end

% 下边界 v = 0
fixed_v = find(abs(y_coor - (-1)) < 1e-5);
for i = 1:length(fixed_v)
    dof = ID(fixed_v(i), 2);
    K(dof, :) = 0;
    K(dof, dof) = 1;
    F(dof) = 0;
end

% 右边界施加拉伸力
right_edge = find(abs(x_coor - 1) < 1e-5);
force_per_node = 10e3 * 0.1; % 10 KPa
for i = 1:length(right_edge)
    dof = ID(right_edge(i), 1);
    F(dof) = F(dof) + force_per_node;
end

% 求解位移
d = K \ F;

% 位移后处理
u = d(1:2:end);
v = d(2:2:end);

% 绘制变形图
figure;
trimesh(IEN, x_coor, y_coor, zeros(size(x_coor)), 'EdgeColor', 'k');
hold on;
quiver(x_coor, y_coor, u, v, 'r');
title('四分之一带孔平板应力分析');
xlabel('X');
ylabel('Y');
axis equal;
