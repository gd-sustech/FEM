function [N, dN_dxi, dN_deta] = Tri3ShapeFunctions(xi, eta)
    % Tri3ShapeFunctions computes the shape functions and their derivatives
    % for a 3-node triangular element in the natural coordinate system.
    %
    % Inputs:
    %   xi  - Natural coordinate xi (ranging from 0 to 1 for a triangle)
    %   eta - Natural coordinate eta (ranging from 0 to 1 for a triangle)
    %
    % Outputs:
    %   N       - Shape function values at (xi, eta) [1x3]
    %   dN_dxi  - Derivative of shape functions w.r.t xi [1x3]
    %   dN_deta - Derivative of shape functions w.r.t eta [1x3]

    % Shape functions for the 3-node triangle
    N = [1 - xi - eta, xi, eta];

    % Derivatives of shape functions w.r.t. xi
    dN_dxi = [-1, 1, 0];

    % Derivatives of shape functions w.r.t. eta
    dN_deta = [-1, 0, 1];
end
