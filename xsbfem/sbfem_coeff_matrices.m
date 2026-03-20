function [E0, E1, E2] = sbfem_coeff_matrices(nodes, sc, D, ngp)
% SBFEM_COEFF_MATRICES - Compute SBFEM coefficient matrices E0, E1, E2
%
% This function computes the coefficient matrices for a polygon subdomain
% in the Scaled Boundary Finite Element Method (SBFEM).
%
% Input:
%   nodes - n x 2 matrix of boundary node coordinates (ordered CCW)
%   sc    - 1 x 2 scaling center coordinates
%   D     - 3 x 3 material constitutive matrix
%   ngp   - number of Gauss points per line element (default: 2)
%
% Output:
%   E0, E1, E2 - SBFEM coefficient matrices
%

if nargin < 4
    ngp = 2;
end

n = size(nodes, 1);     % Number of boundary nodes
ndof = 2 * n;           % Total DOFs

% Shift coordinates relative to scaling center
xy = nodes - repmat(sc, n, 1);

% Initialize coefficient matrices
E0 = zeros(ndof, ndof);
E1 = zeros(ndof, ndof);
E2 = zeros(ndof, ndof);

% Gauss quadrature points and weights
[gp, gw] = gauss_quadrature(ngp);

% Loop over line elements on the boundary
% Each element connects node e to node e+1 (with wrap-around)
for e = 1:n
    % Node indices for this line element
    i1 = e;
    i2 = mod(e, n) + 1;
    
    % Coordinates of element nodes (relative to scaling center)
    x1 = xy(i1, 1); y1 = xy(i1, 2);
    x2 = xy(i2, 1); y2 = xy(i2, 2);
    
    % Global DOF indices for this element
    dof_elem = [2*i1-1, 2*i1, 2*i2-1, 2*i2];
    
    % Gauss integration
    for g = 1:ngp
        eta = gp(g);
        w   = gw(g);
        
        % Shape functions for 2-node line element
        N1 = (1 - eta) / 2;
        N2 = (1 + eta) / 2;
        
        % Shape function derivatives w.r.t. eta
        dN1 = -1/2;
        dN2 =  1/2;
        
        % Coordinate interpolation on boundary
        x_eta = N1*x1 + N2*x2;
        y_eta = N1*y1 + N2*y2;
        
        % Derivatives of coordinates w.r.t. eta
        x_eta_d = dN1*x1 + dN2*x2;  % dx/d_eta
        y_eta_d = dN1*y1 + dN2*y2;  % dy/d_eta
        
        % Jacobian determinant |J| (Eq. 10)
        Jdet = x_eta * y_eta_d - x_eta_d * y_eta;
        absJ = abs(Jdet);
        
        if absJ < 1e-15
            warning('Near-zero Jacobian determinant in element %d', e);
            continue;
        end
        
        % SBFEM strain-displacement matrices b1, b2 (Eq. 9)
        % b1(eta) = (1/|J|) * [y_eta,eta    0      ]
        %                      [0           -x_eta,eta]
        %                      [-x_eta,eta   y_eta,eta]
        b1 = (1/absJ) * [ y_eta_d,       0;
                           0,            -x_eta_d;
                          -x_eta_d,       y_eta_d];
        
        % b2(eta) = (1/|J|) * [-y_eta    0    ]
        %                      [0         x_eta]
        %                      [x_eta    -y_eta]
        b2 = (1/absJ) * [-y_eta,    0;
                           0,        x_eta;
                           x_eta,   -y_eta];
        
        % Shape function matrix N (2 x 4) for displacement interpolation
        Nmat = [N1, 0,  N2, 0;
                0,  N1, 0,  N2];
        
        % Derivative of shape function matrix (2 x 4)
        dNmat = [dN1, 0,   dN2, 0;
                 0,   dN1, 0,   dN2];
        
        % B1 = b1 * N, B2 = b2 * N,eta (Eq. 12)
        B1_elem = b1 * Nmat;     % 3 x 4
        B2_elem = b2 * dNmat;    % 3 x 4
        
        % Accumulate into global coefficient matrices (Eq. 15)
        E0(dof_elem, dof_elem) = E0(dof_elem, dof_elem) + ...
            B1_elem' * D * B1_elem * absJ * w;
        
        E1(dof_elem, dof_elem) = E1(dof_elem, dof_elem) + ...
            B2_elem' * D * B1_elem * absJ * w;
        
        E2(dof_elem, dof_elem) = E2(dof_elem, dof_elem) + ...
            B2_elem' * D * B2_elem * absJ * w;
    end
end

end

function [pts, wts] = gauss_quadrature(n)
% Gauss-Legendre quadrature points and weights on [-1, 1]
switch n
    case 1
        pts = 0;
        wts = 2;
    case 2
        pts = [-1/sqrt(3), 1/sqrt(3)];
        wts = [1, 1];
    case 3
        pts = [-sqrt(3/5), 0, sqrt(3/5)];
        wts = [5/9, 8/9, 5/9];
    case 4
        pts = [-0.861136311594953, -0.339981043584856, ...
                0.339981043584856,  0.861136311594953];
        wts = [0.347854845137454, 0.652145154862546, ...
               0.652145154862546, 0.347854845137454];
    otherwise
        error('Gauss quadrature not implemented for n = %d', n);
end
end
