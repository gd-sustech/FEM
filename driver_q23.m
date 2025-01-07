clear all; clc;

% Problem setup
kappa = 1.0; % conductivity

% Manufactured solution and derivatives
u_exact = @(x, y) x .* (1 - x) .* y .* (1 - y);       % u (exact solution)
u_exact_x = @(x, y) (1 - 2 * x) .* y .* (1 - y);      % du/dx
u_exact_y = @(x, y) x .* (1 - x) .* (1 - 2 * y);      % du/dy
f = @(x, y) 2.0 * kappa * (x .* (1 - x) + y .* (1 - y)); % source term

% Sobolev norms (integrate error function)
compute_error = @(u_h, x_coor, y_coor, IEN, n_el, n_en, quad_pts, w) ...
    calculate_error(u_h, u_exact, u_exact_x, u_exact_y, x_coor, y_coor, IEN, n_el, n_en, quad_pts, w);

% Convergence analysis
n_el_x_list = [10, 20, 40, 80]; % Varying mesh sizes
errors_L2 = zeros(length(n_el_x_list), 1);
errors_H1 = zeros(length(n_el_x_list), 1);

for i = 1:length(n_el_x_list)
    n_el_x = n_el_x_list(i);
    n_el_y = n_el_x;
    
    % Generate mesh for quadrilateral elements
    [x_coor, y_coor, IEN, n_el, n_np, n_en, n_eq, ID, LM] = generate_mesh_quad(n_el_x, n_el_y);
    
    % Solve the problem
    u_h = solve_fem(x_coor, y_coor, IEN, n_el, n_en, f, kappa, n_eq, ID, LM);
    
    % Calculate errors
    [L2_error, H1_error] = compute_error(u_h, x_coor, y_coor, IEN, n_el, n_en, 3, 1/3);
    errors_L2(i) = L2_error;
    errors_H1(i) = H1_error;
end

% Log-log plot
figure;
loglog(1 ./ n_el_x_list, errors_L2, '-o', 'DisplayName', 'L2 Error');
hold on;
loglog(1 ./ n_el_x_list, errors_H1, '-x', 'DisplayName', 'H1 Error');
xlabel('Mesh size (h)');
ylabel('Error');
legend;
title('Convergence rate analysis for FEM');
grid on;
