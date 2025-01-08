function [N, dN_dxi, dN_deta] = Quad4ShapeFunctions(xi, eta)
% Quad4ShapeFunctions: Compute shape functions and their derivatives
% for a 4-node quadrilateral element at a given (xi, eta) point.
% 
% INPUTS:
%   xi, eta - Coordinates of the Gauss point in the natural coordinate system.
% 
% OUTPUTS:
%   N        - Shape functions (1x4 array)
%   dN_dxi   - Derivatives of shape functions w.r.t xi (1x4 array)
%   dN_deta  - Derivatives of shape functions w.r.t eta (1x4 array)

% Shape functions
N = [(1 - xi) * (1 - eta), ...  % N1
     (1 + xi) * (1 - eta), ...  % N2
     (1 + xi) * (1 + eta), ...  % N3
     (1 - xi) * (1 + eta)] / 4;

% Derivatives of shape functions w.r.t xi
dN_dxi = [-(1 - eta), ...  % dN1/dxi
           (1 - eta), ...  % dN2/dxi
           (1 + eta), ...  % dN3/dxi
          -(1 + eta)] / 4;

% Derivatives of shape functions w.r.t eta
dN_deta = [-(1 - xi), ...  % dN1/deta
           -(1 + xi), ...  % dN2/deta
            (1 + xi), ...  % dN3/deta
            (1 - xi)] / 4;
end
