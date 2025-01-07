function [L2_error, H1_error] = calculate_error(u_h, u_exact, u_exact_x, u_exact_y, x_coor, y_coor, IEN, n_el, n_en, quad_pts, weight)
    L2_error = 0;
    H1_error = 0;
    
    for ee = 1:n_el
        % Element coordinates
        x_ele = x_coor(IEN(ee, :));
        y_ele = y_coor(IEN(ee, :));
        
        for qp = 1:length(quad_pts)
            [xi, eta, w] = GaussTriangle(quad_pts); % Get quadrature data
            [N, dN_dxi, dN_deta] = TriangleShapeFunctions(xi(qp), eta(qp));
            
            % Compute derivatives
            dx_dxi = 0; dx_deta = 0;
            dy_dxi = 0; dy_deta = 0;
            for aa = 1:n_en
                dx_dxi = dx_dxi + x_ele(aa) * dN_dxi(aa);
                dx_deta = dx_deta + x_ele(aa) * dN_deta(aa);
                dy_dxi = dy_dxi + y_ele(aa) * dN_dxi(aa);
                dy_deta = dy_deta + y_ele(aa) * dN_deta(aa);
            end
            
            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi; % Jacobian determinant
            
            % Map shape function derivatives to physical domain
            dN_dx = (dN_dxi * dy_deta - dN_deta * dy_dxi) / detJ;
            dN_dy = (-dN_dxi * dx_deta + dN_deta * dx_dxi) / detJ;
            
            % Evaluate solution and exact derivatives
            u_h_val = 0;
            du_h_dx = 0;
            du_h_dy = 0;
            for aa = 1:n_en
                u_h_val = u_h_val + u_h(IEN(ee, aa)) * N(aa);
                du_h_dx = du_h_dx + u_h(IEN(ee, aa)) * dN_dx(aa);
                du_h_dy = du_h_dy + u_h(IEN(ee, aa)) * dN_dy(aa);
            end
            
            u_ex = u_exact(xi(qp), eta(qp));
            u_ex_x = u_exact_x(xi(qp), eta(qp));
            u_ex_y = u_exact_y(xi(qp), eta(qp));
            
            % Accumulate L2 and H1 errors
            L2_error = L2_error + w * detJ * (u_h_val - u_ex)^2;
            H1_error = H1_error + w * detJ * ((du_h_dx - u_ex_x)^2 + (du_h_dy - u_ex_y)^2);
        end
    end
    
    L2_error = sqrt(L2_error);
    H1_error = sqrt(H1_error);
end
