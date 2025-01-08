clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% Error norms
L2_errors = [];
H1_errors = [];
mesh_sizes = [];

% quadrature rule
n_int = 3; % use a 3-point Gaussian quadrature rule for triangles
[xi, eta, weight] = GaussTriangle(n_int);

% mesh generation
n_en   = 3;               % number of nodes in a triangle element
for n_el_x = [20, 40, 60, 80, 100]
    n_el_y = n_el_x;
    n_quad = n_el_x * n_el_y; % total number of quadrilateral elements

n_el = n_quad * 2;        % each quadrilateral is split into 2 triangles
n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir
mesh_sizes = [mesh_sizes, hx];
% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array for triangular elements
IEN = zeros(n_el, n_en);
tri_idx = 0;

for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = (ey-1) * n_el_x + ex; % quadrilateral element index

    % Nodes of the quadrilateral
    n1 = (ey-1) * n_np_x + ex;
    n2 = n1 + 1;
    n3 = n1 + n_np_x + 1;
    n4 = n1 + n_np_x;

    % Subdivide quadrilateral into two triangles
    tri_idx = tri_idx + 1;
    IEN(tri_idx, :) = [n1, n2, n3]; % first triangle

    tri_idx = tri_idx + 1;
    IEN(tri_idx, :) = [n1, n3, n4]; % second triangle
  end
end

% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele = zeros(n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dksi = 0.0; dx_deta = 0.0;
    dy_dksi = 0.0; dy_deta = 0.0;
    
    % Loop over nodes to compute Jacobian and shape function contributions
    for aa = 1 : n_en
      [N, dN_dksi, dN_deta] = TriangleShapeFunctions(xi(ll), eta(ll));
      x_l = x_l + x_ele(aa) * N(aa);
      y_l = y_l + y_ele(aa) * N(aa);
      dx_dksi  = dx_dksi  + x_ele(aa) * dN_dksi(aa);
      dx_deta  = dx_deta  + x_ele(aa) * dN_deta(aa);
      dy_dksi  = dy_dksi  + y_ele(aa) * dN_dksi(aa);
      dy_deta  = dy_deta  + y_ele(aa) * dN_deta(aa);
    end
    
    detJ = dx_dksi * dy_deta - dx_deta * dy_dksi;
    
    for aa = 1 : n_en
      [N, dN_dksi, dN_deta] = TriangleShapeFunctions(xi(ll), eta(ll));
      Na_x = (dN_dksi(aa) * dy_deta - dN_deta(aa) * dy_dksi) / detJ;
      Na_y = (-dN_dksi(aa) * dx_deta + dN_deta(aa) * dx_dksi) / detJ;
      
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * N(aa);
      
      for bb = 1 : n_en
        Nb_x = (dN_dksi(bb) * dy_deta - dN_deta(bb) * dy_dksi) / detJ;
        Nb_y = (-dN_dksi(bb) * dx_deta + dN_deta(bb) * dx_dksi) / detJ;
        k_ele(aa, bb) = k_ele(aa, bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end
    end
  end
  
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        end
      end
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% Full displacement vector
disp = zeros(n_np, 1);
for ii = 1:n_np
    if ID(ii) > 0
        disp(ii) = dn(ID(ii));
    end
end

% Error calculation
L2_error = 0.0;
H1_error = 0.0;

for ee = 1:n_el
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));
    d_ele = disp(IEN(ee, :));

    for ll = 1:n_int
        x_l = 0.0; y_l = 0.0;
        dx_dksi = 0.0; dx_deta = 0.0;
        dy_dksi = 0.0; dy_deta = 0.0;

        for aa = 1:n_en
            [N, dN_dksi, dN_deta] = TriangleShapeFunctions(xi(ll), eta(ll));
            x_l = x_l + x_ele(aa) * N(aa);
            y_l = y_l + y_ele(aa) * N(aa);
            dx_dksi = dx_dksi + x_ele(aa) * dN_dksi(aa);
            dx_deta = dx_deta + x_ele(aa) * dN_deta(aa);
            dy_dksi = dy_dksi + y_ele(aa) * dN_dksi(aa);
            dy_deta = dy_deta + y_ele(aa) * dN_deta(aa);
        end

        detJ = dx_dksi * dy_deta - dx_deta * dy_dksi;

        u_h = 0.0; u_hx = 0.0; u_hy = 0.0;
        for aa = 1:n_en
            [N, dN_dksi, dN_deta] = TriangleShapeFunctions(xi(ll), eta(ll));
            Na_x = (dN_dksi(aa) * dy_deta - dN_deta(aa) * dy_dksi) / detJ;
            Na_y = (-dN_dksi(aa) * dx_deta + dN_deta(aa) * dx_dksi) / detJ;
            u_h = u_h + d_ele(aa) * N(aa);
            u_hx = u_hx + d_ele(aa) * Na_x;
            u_hy = u_hy + d_ele(aa) * Na_y;
        end

        u_exact = exact(x_l, y_l);
        u_exact_x = exact_x(x_l, y_l);
        u_exact_y = exact_y(x_l, y_l);

        L2_error = L2_error + weight(ll) * detJ * (u_h - u_exact)^2;
        H1_error = H1_error + weight(ll) * detJ * ((u_h - u_exact)^2 + (u_hx - u_exact_x)^2 + (u_hy - u_exact_y)^2);
    end
end

L2_errors = [L2_errors, sqrt(L2_error)];
H1_errors = [H1_errors, sqrt(H1_error)];
end
% Plot convergence
figure;
loglog(mesh_sizes, L2_errors, '-o', 'DisplayName', 'L2 Norm');
hold on;
loglog(mesh_sizes, H1_errors, '-x', 'DisplayName', 'H1 Norm');
xlabel('Mesh Size (h)');
ylabel('Error');
legend show;
title('Convergence of FEM Solution (Triangular)');
grid on;

