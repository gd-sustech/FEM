function [N, dN_dxi, dN_deta] = TriangleShapeFunctions(xi, eta)
    % TriangleShapeFunctions computes the shape functions and their derivatives 
    % for linear triangular elements.
    %
    % Input:
    %   xi  - Natural coordinate in the xi direction
    %   eta - Natural coordinate in the eta direction
    %
    % Output:
    %   N        - Shape functions [N1, N2, N3]
    %   dN_dxi   - Derivative of N with respect to xi [dN1/dxi, dN2/dxi, dN3/dxi]
    %   dN_deta  - Derivative of N with respect to eta [dN1/deta, dN2/deta, dN3/deta]

    % Shape functions for a linear triangular element
    N = [1 - xi - eta; xi; eta];

    % Derivatives of shape functions with respect to xi and eta
    dN_dxi = [-1; 1; 0];   % dN1/dxi, dN2/dxi, dN3/dxi
    dN_deta = [-1; 0; 1];  % dN1/deta, dN2/deta, dN3/deta
end
