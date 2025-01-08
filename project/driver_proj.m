clear all; clc;

% Material properties
E = 10e9;  % Elastic modulus (Pa)
nu = 0.3;  % Poisson's ratio

% Plane stress constitutive matrix
C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

% Geometric parameters
L = 4;   % Length of the square plate (m)
r_hole = 0.5; % Radius of the circular hole (m)

% Mesh parameters
n_el_x = 60;  % Number of elements in x-direction
n_el_y = 60;  % Number of elements in y-direction
n_el = n_el_x * n_el_y; % Total number of elements

n_np_x = n_el_x + 1;      % Number of nodal points in x-direction
n_np_y = n_el_y + 1;      % Number of nodal points in y-direction
n_np = n_np_x * n_np_y;   % Total number of nodal points

hx = L / n_el_x;          % Element size in x-direction
hy = L / n_el_y;          % Element size in y-direction

% Generate nodal coordinates
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);
for ny = 1:n_np_y
    for nx = 1:n_np_x
        index = (ny - 1) * n_np_x + nx;
        x_coor(index) = (nx - 1) * hx;
        y_coor(index) = (ny - 1) * hy;
    end
end

% Remove nodes inside the hole
nodes_to_remove = find((x_coor - L/2).^2 + (y_coor - L/2).^2 < r_hole^2);
x_coor(nodes_to_remove) = NaN;
y_coor(nodes_to_remove) = NaN;

% Adjust IEN array (element connectivity)
IEN = zeros(n_el, 4);
for ex = 1:n_el_x
    for ey = 1:n_el_y
        ee = (ey - 1) * n_el_x + ex; % Element index
        n1 = (ey - 1) * n_np_x + ex;
        n2 = n1 + 1;
        n3 = n2 + n_np_x;
        n4 = n1 + n_np_x;
        
        if any(ismember([n1, n2, n3, n4], nodes_to_remove))
            IEN(ee, :) = NaN; % Mark as invalid element
        else
            IEN(ee, :) = [n1, n2, n3, n4];
        end
    end
end
IEN = IEN(~isnan(IEN(:, 1)), :); % Remove invalid elements

% ID array (degrees of freedom)
ID = zeros(n_np, 2); % Two DOFs per node (u, v)
counter = 0;
for i = 1:n_np
    if ~ismember(i, nodes_to_remove) % Only include valid nodes
        counter = counter + 1;
        ID(i, 1) = counter;   % u DOF
        counter = counter + 1;
        ID(i, 2) = counter;   % v DOF
    end
end
n_eq = max(ID(:)); % Total number of equations

% LM array (local to global mapping)
LM = zeros(size(IEN, 1), 8); % 4 nodes * 2 DOFs per node
for ee = 1:size(IEN, 1)
    for aa = 1:4
        node = IEN(ee, aa);
        LM(ee, (aa-1)*2 + 1) = ID(node, 1);
        LM(ee, (aa-1)*2 + 2) = ID(node, 2);
    end
end

% Allocate global stiffness matrix and force vector
K = spalloc(n_eq, n_eq, 64 * size(IEN, 1));
F = zeros(n_eq, 1);

% Quadrature rule
n_int = 2; % Number of integration points per direction
[xi, eta, weight] = Gauss2D(n_int, n_int);

% Assembly process
for ee = 1:size(IEN, 1)
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));
    
    k_ele = zeros(8, 8); % Element stiffness matrix
    for ll = 1:n_int^2
        % Map Gauss points to physical domain
        [N, dN_dxi, dN_deta] = Quad4ShapeFunctions(xi(ll), eta(ll));
        J = [dN_dxi * x_ele, dN_dxi * y_ele; dN_deta * x_ele, dN_deta * y_ele];
        detJ = det(J);
        invJ = inv(J);
        dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
        dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;
        
        B = zeros(3, 8); % Strain-displacement matrix
        for i = 1:4
            B(:, (i-1)*2+1:(i-1)*2+2) = [dN_dx(i), 0; 0, dN_dy(i); dN_dy(i), dN_dx(i)];
        end
        
        % Element stiffness matrix
        k_ele = k_ele + B' * C * B * detJ * weight(ll);
    end
    
    % Assemble global stiffness matrix
    for i = 1:8
        for j = 1:8
            if LM(ee, i) > 0 && LM(ee, j) > 0
                K(LM(ee, i), LM(ee, j)) = K(LM(ee, i), LM(ee, j)) + k_ele(i, j);
            end
        end
    end
end

% Apply boundary conditions and external forces
% Left edge: u = 0 (fixed in x-direction)
fixed_u = find(x_coor == 0 & ~isnan(x_coor));
for i = 1:length(fixed_u)
    dof = ID(fixed_u(i), 1);
    K(dof, :) = 0;
    K(dof, dof) = 1;
    F(dof) = 0;
end

% Bottom edge: v = 0 (fixed in y-direction)
fixed_v = find(y_coor == 0 & ~isnan(y_coor));
for i = 1:length(fixed_v)
    dof = ID(fixed_v(i), 2);
    K(dof, :) = 0;
    K(dof, dof) = 1;
    F(dof) = 0;
end

% Right edge: Apply uniform tension in x-direction
right_edge = find(x_coor == L & ~isnan(x_coor));
force_per_node = 10e3* hy; % 10 KPa tension
for i = 1:length(right_edge)
    dof = ID(right_edge(i), 1);
    F(dof) = F(dof) + force_per_node;
end

% Solve for displacements
d = K \ F;

% Post-processing: Insert displacements back into full vector
u = zeros(n_np, 1); % Displacement in x-direction
v = zeros(n_np, 1); % Displacement in y-direction
for i = 1:n_np
    if ID(i, 1) > 0
        u(i) = d(ID(i, 1));
    end
    if ID(i, 2) > 0
        v(i) = d(ID(i, 2));
    end
end

% Save results
save('StressAnalysis', 'u', 'v', 'x_coor', 'y_coor', 'IEN',"n_el_x", "n_el_y");
