function [K, data] = sbfem_stiffness(nodes, sc, D, ngp)
% SBFEM_STIFFNESS - Compute stiffness matrix of an SBFEM subdomain
%
% Computes the stiffness matrix for a polygon subdomain using the 
% Scaled Boundary Finite Element Method via Schur decomposition.
%
% Input:
%   nodes - n x 2 matrix of boundary node coordinates (ordered CCW)
%   sc    - 1 x 2 scaling center coordinates
%   D     - 3 x 3 material constitutive matrix
%   ngp   - number of Gauss points (default: 2)
%
% Output:
%   K    - 2n x 2n stiffness matrix
%   data - structure with intermediate results for SIF computation

if nargin < 4
    ngp = 2;
end

n = size(nodes, 1);
ndof = 2 * n;

% Step 1: Compute coefficient matrices E0, E1, E2
[E0, E1, E2] = sbfem_coeff_matrices(nodes, sc, D, ngp);

% Step 2: Form Hamiltonian matrix Z (Eq. 17)
% Z = [E0^{-1} E1^T,        -E0^{-1}           ]
%     [E1 E0^{-1} E1^T - E2, -E1 E0^{-1}        ]
E0inv = E0 \ eye(ndof);

Z = zeros(2*ndof, 2*ndof);
Z(1:ndof, 1:ndof)           =  E0inv * E1';
Z(1:ndof, ndof+1:2*ndof)    = -E0inv;
Z(ndof+1:2*ndof, 1:ndof)    =  E1 * E0inv * E1' - E2;
Z(ndof+1:2*ndof, ndof+1:2*ndof) = -E1 * E0inv;

% Step 3: Schur decomposition (Eq. 18): ZV = VS
[V, S] = schur(Z, 'complex');

% Sort eigenvalues: negative real parts first (bounded domain modes)
eigvals = diag(S);
[~, idx] = sort(real(eigvals));
idx_neg = idx(real(eigvals(idx)) < 0);
idx_pos = idx(real(eigvals(idx)) >= 0);
reorder = [idx_neg; idx_pos];

% Reorder Schur form
select = false(2*ndof, 1);
select(idx_neg) = true;
try
    [V, S] = ordschur(V, S, select);
catch
    % Manual reordering if ordschur fails
    V = V(:, reorder);
    S = V \ (Z * V);
end

% Step 4: Extract submatrices (Eq. 20)
Vu = V(1:ndof, 1:ndof);         % Modal displacements (bounded)
Vq = V(ndof+1:2*ndof, 1:ndof);  % Modal forces (bounded)
Sn = S(1:ndof, 1:ndof);         % Eigenvalue block (negative real parts)

% Step 5: Compute stiffness matrix K = Vq * Vu^{-1} (Eq. 24)
K = real(Vq / Vu);

% Symmetrize (should be symmetric in theory)
K = (K + K') / 2;

% Store data for SIF computation
data.Vu = Vu;
data.Vq = Vq;
data.Sn = Sn;
data.E0 = E0;
data.E1 = E1;
data.E2 = E2;
data.Z  = Z;
data.V  = V;
data.S  = S;
data.nodes = nodes;
data.sc = sc;
data.D  = D;

end
