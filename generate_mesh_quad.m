function [x_coor, y_coor, IEN, n_el, n_np, n_en, n_eq, ID, LM] = generate_mesh_quad(n_el_x, n_el_y)
    % INPUT:
    % n_el_x, n_el_y: number of elements in x and y directions
    % 
    % OUTPUT:
    % x_coor, y_coor: coordinates of all nodes
    % IEN: connectivity matrix (element to node mapping)
    % n_el: total number of elements
    % n_np: total number of nodes
    % n_en: number of nodes per element
    % n_eq: number of unknowns
    % ID: global node-to-equation mapping
    % LM: local element-to-global equation mapping
    
    % Number of nodes in x and y directions
    n_np_x = n_el_x + 1; % nodes in x direction
    n_np_y = n_el_y + 1; % nodes in y direction
    n_np = n_np_x * n_np_y; % total number of nodes
    
    % Total number of elements
    n_el = n_el_x * n_el_y;
    n_en = 4; % each quadrilateral element has 4 nodes

    % Generate node coordinates
    x_coor = zeros(n_np, 1); % x-coordinates of nodes
    y_coor = zeros(n_np, 1); % y-coordinates of nodes
    
    hx = 1.0 / n_el_x; % element size in x direction
    hy = 1.0 / n_el_y; % element size in y direction
    
    % Assign coordinates to each node
    for ny = 1:n_np_y
        for nx = 1:n_np_x
            index = (ny - 1) * n_np_x + nx; % node index
            x_coor(index) = (nx - 1) * hx; % x-coordinate
            y_coor(index) = (ny - 1) * hy; % y-coordinate
        end
    end
    
    % Generate IEN array (element to node mapping)
    IEN = zeros(n_el, n_en); % connectivity matrix
    for ex = 1:n_el_x
        for ey = 1:n_el_y
            ee = (ey - 1) * n_el_x + ex; % element index
            n1 = (ey - 1) * n_np_x + ex; % bottom-left node
            n2 = n1 + 1;                 % bottom-right node
            n3 = n1 + n_np_x + 1;        % top-right node
            n4 = n1 + n_np_x;            % top-left node
            IEN(ee, :) = [n1, n2, n3, n4];
        end
    end
    
    % Generate ID array (global node-to-equation mapping)
    ID = zeros(n_np, 1); % initialize
    counter = 0;
    for ny = 2:n_np_y - 1
        for nx = 2:n_np_x - 1
            index = (ny - 1) * n_np_x + nx; % node index
            counter = counter + 1;          % equation number
            ID(index) = counter;  
        end
    end
    n_eq = counter; % total number of equations

    % Generate LM array (local element-to-global equation mapping)
    LM = ID(IEN); % map local element nodes to global equations
end
