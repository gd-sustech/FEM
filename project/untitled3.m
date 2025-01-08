% 加载结果数据
load('StressAnalysis.mat', 'u', 'v', 'x_coor', 'y_coor', 'IEN');

% 放大位移，以便可视化
u_scaled = 10 * u; % 放大10倍
v_scaled = 10 * v;

% 绘制变形后的网格
figure;
hold on;
for ee = 1:size(IEN, 1)
    % 获取当前单元的节点坐标
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));
    
    % 根据位移计算变形后的节点位置
    x_ele_deformed = x_ele + u_scaled(IEN(ee, :));
    y_ele_deformed = y_ele + v_scaled(IEN(ee, :));
    
    % 绘制变形后的单元
    fill(x_ele_deformed, y_ele_deformed, 'w', 'EdgeColor', 'k');
end

% 设置图形
title('Deformed Mesh');
axis equal;
xlabel('X (m)');
ylabel('Y (m)');

