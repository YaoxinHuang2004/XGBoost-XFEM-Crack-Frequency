function [E0, E1, E2] = sbfem_coeff_matrices_crack(nodes, sc, D, ngp)
% SBFEM_COEFF_MATRICES_CRACK - Coefficient matrices for crack subdomain
%
% For a crack subdomain, the boundary is OPEN: nodes 1..n define n-1 line
% elements (node n does NOT connect back to node 1). Nodes 1 and n are the
% upper and lower crack mouth nodes at the same physical location.
% The scaling center is at the crack tip.
%
% Input:
%   nodes - n x 2 boundary node coordinates (ordered: upper_mouth -> CCW -> lower_mouth)
%   sc    - 1 x 2 scaling center (crack tip position)
%   D     - 3 x 3 material constitutive matrix
%   ngp   - Gauss points per element (default: 3)
%
% Output:
%   E0, E1, E2 - coefficient matrices (2n x 2n)

if nargin < 4
    ngp = 3;
end

n = size(nodes, 1);
ndof = 2 * n;

xy = nodes - repmat(sc, n, 1);

E0 = zeros(ndof, ndof);
E1 = zeros(ndof, ndof);
E2 = zeros(ndof, ndof);

[gp, gw] = gauss_quadrature_local(ngp);

% Open boundary: n-1 elements (no wrap-around from node n to node 1)
n_elem = n - 1;

for e = 1:n_elem
    i1 = e;
    i2 = e + 1;
    
    x1 = xy(i1,1); y1 = xy(i1,2);
    x2 = xy(i2,1); y2 = xy(i2,2);
    
    dof_elem = [2*i1-1, 2*i1, 2*i2-1, 2*i2];
    
    for g = 1:ngp
        eta = gp(g);
        w   = gw(g);
        
        N1 = (1 - eta)/2;
        N2 = (1 + eta)/2;
        dN1 = -1/2;
        dN2 =  1/2;
        
        x_eta   = N1*x1 + N2*x2;
        y_eta   = N1*y1 + N2*y2;
        x_eta_d = dN1*x1 + dN2*x2;
        y_eta_d = dN1*y1 + dN2*y2;
        
        Jdet = x_eta * y_eta_d - x_eta_d * y_eta;
        absJ = abs(Jdet);
        
        if absJ < 1e-15
            continue;
        end
        
        b1 = (1/absJ) * [ y_eta_d,   0;
                           0,         -x_eta_d;
                          -x_eta_d,    y_eta_d];
        
        b2 = (1/absJ) * [-y_eta,   0;
                           0,       x_eta;
                           x_eta,  -y_eta];
        
        Nmat  = [N1, 0,  N2, 0;  0, N1, 0, N2];
        dNmat = [dN1, 0, dN2, 0; 0, dN1, 0, dN2];
        
        B1_e = b1 * Nmat;
        B2_e = b2 * dNmat;
        
        E0(dof_elem, dof_elem) = E0(dof_elem, dof_elem) + B1_e' * D * B1_e * absJ * w;
        E1(dof_elem, dof_elem) = E1(dof_elem, dof_elem) + B2_e' * D * B1_e * absJ * w;
        E2(dof_elem, dof_elem) = E2(dof_elem, dof_elem) + B2_e' * D * B2_e * absJ * w;
    end
end

end

function [pts, wts] = gauss_quadrature_local(n)
switch n
    case 1
        pts = 0; wts = 2;
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
    case 5
        pts = [-0.906179845938664, -0.538469310105683, 0, ...
                0.538469310105683,  0.906179845938664];
        wts = [0.236926885056189, 0.478628670499366, 0.568888888888889, ...
               0.478628670499366, 0.236926885056189];
    otherwise
        error('Not implemented for n=%d', n);
end
end
