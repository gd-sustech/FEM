% Load the results from StressAnalysis.mat
load('StressAnalysis.mat', 'u', 'v', 'x_coor', 'y_coor', 'IEN');

% Remove invalid nodes (those set to NaN)
valid_nodes = ~isnan(x_coor) & ~isnan(y_coor);
x_coor_valid = x_coor(valid_nodes);
y_coor_valid = y_coor(valid_nodes);
u_valid = u(valid_nodes);
v_valid = v(valid_nodes);

% Plot the deformed shape
scale_factor = 1e3; % Scale factor for visualization
x_deformed = x_coor_valid + scale_factor * u_valid;
y_deformed = y_coor_valid + scale_factor * v_valid;

figure;
subplot(1, 2, 1);
scatter(x_coor_valid, y_coor_valid, 10, 'filled');
title('Original Mesh');
xlabel('x (m)'); ylabel('y (m)');
axis equal;

tiledlayout(1,2);
nexttile;
plot(x_deformed,);...
