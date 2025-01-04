function [xi, eta, weight] = GaussTriangle(n_int)
    % GaussTriangle generates Gauss quadrature points and weights for 
    % integration over a triangular element.
    %
    % Input:
    %   n_int - Number of integration points (1, 3, or 4 supported here)
    %
    % Output:
    %   xi     - Natural coordinates in xi direction
    %   eta    - Natural coordinates in eta direction
    %   weight - Corresponding weights for each integration point

    if n_int == 1
        % Single integration point (degree of precision = 1)
        xi = 1/3;
        eta = 1/3;
        weight = 1.0;

    elseif n_int == 3
        % Three integration points (degree of precision = 2)
        xi = [1/6; 2/3; 1/6];
        eta = [1/6; 1/6; 2/3];
        weight = [1/3; 1/3; 1/3];

    elseif n_int == 4
        % Four integration points (degree of precision = 3)
        xi = [1/3; 0.6; 0.2; 0.2];
        eta = [1/3; 0.2; 0.6; 0.2];
        weight = [-27/48; 25/48; 25/48; 25/48];

    else
        error('Unsupported number of integration points. Use 1, 3, or 4.');
    end
end
