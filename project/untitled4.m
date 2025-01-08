% 绘制位移的等高线图
figure;
contourf(x_coor, y_coor, u, 20, 'LineStyle', 'none');
colorbar;
title('Displacement in X Direction');
xlabel('X (m)');
ylabel('Y (m)');
axis equal;
