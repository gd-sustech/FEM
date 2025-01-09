% 加载网格数据
% load('mesh.m'); % 确保文件在当前目录下

% 提取节点和单元信息
nodes = msh.POS(:, 1:2); % 节点坐标 (x, y)
elements = msh.TRIANGLES(:, 1:3); % 三角形单元 (节点索引)

% 绘制网格
figure;
trimesh(elements, nodes(:, 1), nodes(:, 2), zeros(size(nodes, 1), 1), 'EdgeColor', 'k');
axis equal;
title('四分之一带孔平板网格');
xlabel('X');
ylabel('Y');
