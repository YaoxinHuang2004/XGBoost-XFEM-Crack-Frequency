function [K, KI, KII, data] = sbfem_crack_tip(nodes, sc, D, ub, E_mod, nu, plane_type, ngp)
% SBFEM_CRACK_TIP - SBFEM stiffness and SIF for crack tip subdomain
%
% Input:
%   nodes      - n x 2 boundary nodes (ordered: lower_mouth -> CCW -> upper_mouth)
%   sc         - 1 x 2 crack tip position (scaling center)
%   D          - 3 x 3 constitutive matrix
%   ub         - 2n x 1 boundary displacements ([] to skip SIF)
%   E_mod, nu  - elastic constants
%   plane_type - 'strain' or 'stress'
%   ngp        - Gauss points (default: 3)

if nargin < 8, ngp = 3; end
if nargin < 7, plane_type = 'strain'; end
if nargin < 6, nu = 0.3; end
if nargin < 5, E_mod = 1; end
if nargin < 4, ub = []; end

n = size(nodes, 1);
ndof = 2 * n;

%% Step 1: Coefficient matrices (open boundary: n-1 line elements)
[E0, E1, E2] = sbfem_coeff_matrices_crack(nodes, sc, D, ngp);

%% Step 2: Regularize E0
reg = 1e-12 * max(abs(diag(E0)));
if reg < eps, reg = 1e-12; end
E0reg = E0 + reg * eye(ndof);

%% Step 3: Hamiltonian Z (Eq. 17)
E0inv = E0reg \ eye(ndof);
Z = [ E0inv * E1',              -E0inv;
      E1 * E0inv * E1' - E2,    -E1 * E0inv ];

%% Step 4: Eigendecomposition
[Phi, Lambda] = eig(Z);
eigvals_all = diag(Lambda);

% Sort by real part
[~, idx_sort] = sort(real(eigvals_all));
eigvals_all = eigvals_all(idx_sort);
Phi = Phi(:, idx_sort);

% Select N eigenvalues with most negative real parts (bounded domain)
eigvals_bounded = eigvals_all(1:ndof);
Phi_bounded = Phi(:, 1:ndof);

% Modal matrices
Vu = Phi_bounded(1:ndof, :);
Vq = Phi_bounded(ndof+1:end, :);

%% Step 5: Stiffness K = Vq * Vu^{-1} (Eq. 24)
K = real(Vq / Vu);
K = (K + K') / 2;

%% Step 6: Identify singular eigenvalues: -1 < Re(lambda) < 0
sing_mask = (real(eigvals_bounded) > -0.95) & (real(eigvals_bounded) < -0.05);
sing_idx = find(sing_mask);

% Store data
data.Vu = Vu;
data.Vq = Vq;
data.eigvals_bounded = eigvals_bounded;
data.sing_idx = sing_idx;
data.sing_eigvals = eigvals_bounded(sing_idx);
data.nodes = nodes;
data.sc = sc;

%% Step 7: SIF computation
KI = NaN; KII = NaN;
data.KI_disp = NaN; data.KII_disp = NaN;
data.KI_stress = NaN; data.KII_stress = NaN;

if isempty(ub) || isempty(sing_idx)
    return;
end

% Coordinates relative to crack tip
xy = nodes - repmat(sc, n, 1);

% ---- Determine local crack coordinate system ----
mouth_vec = (xy(1,:) + xy(n,:)) / 2;
mouth_dir = mouth_vec / max(norm(mouth_vec), 1e-15);
ext_dir = -mouth_dir;                        % Crack extension direction
crack_normal = [-ext_dir(2), ext_dir(1)];    % Perpendicular

% Ensure node n (upper face) is on + normal side
if dot(xy(n,:), crack_normal) < dot(xy(1,:), crack_normal)
    crack_normal = -crack_normal;
end

data.ext_dir = ext_dir;
data.crack_normal = crack_normal;

% Integration constants: c = Vu^{-1} * ub (Eq. 22)
c_all = Vu \ ub;

% Singular modes
Vu_s = Vu(:, sing_idx);
c_s  = c_all(sing_idx);
lambda_s = eigvals_bounded(sing_idx);
Lam_s = diag(lambda_s);

% Singular displacement at boundary (Eq. 28)
u_sing = real(Vu_s * c_s);

% ---- Displacement-based SIF (Eq. 26) ----
delta_ux = u_sing(2*n-1) - u_sing(1);
delta_uy = u_sing(2*n)   - u_sing(2);

delta_opening = dot([delta_ux, delta_uy], crack_normal);
delta_sliding = dot([delta_ux, delta_uy], ext_dir);

r0 = (norm(xy(1,:)) + norm(xy(n,:))) / 2;

G = E_mod / (2*(1+nu));
if strcmpi(plane_type, 'strain')
    kappa = 3 - 4*nu;
else
    kappa = (3-nu) / (1+nu);
end

factor = G / (kappa + 1) * sqrt(2*pi / r0);
KI  = factor * delta_opening;
KII = factor * delta_sliding;

data.KI_disp = KI;
data.KII_disp = KII;
data.r0 = r0;
data.delta_opening = delta_opening;
data.delta_sliding = delta_sliding;
data.c_all = c_all;
data.c_s = c_s;
data.u_sing = u_sing;

% ---- Stress-based SIF (Eq. 32) ----
% Key improvement: directly evaluate at theta=0 without interpolation
[KI_s, KII_s] = stress_based_sif_direct(xy, D, Vu_s, Lam_s, c_s, ext_dir, crack_normal);

data.KI_stress = KI_s;
data.KII_stress = KII_s;

end


function [KI, KII] = stress_based_sif_direct(xy, D, Vu_s, Lam_s, c_s, ext_dir, crack_normal)
% STRESS_BASED_SIF_DIRECT - Evaluate stress mode directly at theta=0
%
% Instead of sampling at Gauss points and interpolating, this function:
% 1. Computes theta at each boundary NODE
% 2. Finds the element that contains theta=0 (crack front direction)
% 3. Uses bisection to find the exact eta within that element where theta=0
% 4. Evaluates the stress mode at that exact point
%
% This avoids all interpolation errors.

n = size(xy, 1);
n_elem = n - 1;
KI = NaN; KII = NaN;

% Compute angle of each boundary node in local crack coordinates
theta_nodes = zeros(n, 1);
for i = 1:n
    px = dot(xy(i,:), ext_dir);
    py = dot(xy(i,:), crack_normal);
    theta_nodes(i) = atan2(py, px);
end

% Find element where theta crosses 0 (from negative to positive)
target_elem = [];
for e = 1:n_elem
    th1 = theta_nodes(e);
    th2 = theta_nodes(e+1);
    % theta=0 is inside this element if th1 and th2 straddle 0
    if (th1 <= 0 && th2 >= 0) || (th1 >= 0 && th2 <= 0)
        % Avoid elements near crack face (|theta| near pi)
        if abs(th1) < 0.9*pi && abs(th2) < 0.9*pi
            target_elem = e;
            break;
        end
    end
end

if isempty(target_elem)
    % Fallback: find element with node angles closest to 0
    [~, closest_node] = min(abs(theta_nodes(2:end-1)));
    closest_node = closest_node + 1;  % Offset for skipping node 1
    target_elem = max(1, closest_node - 1);
end

% Bisection to find eta where theta(eta) = 0 within the target element
i1 = target_elem;
i2 = target_elem + 1;
x1 = xy(i1,1); y1 = xy(i1,2);
x2 = xy(i2,1); y2 = xy(i2,2);
dof_e = [2*i1-1, 2*i1, 2*i2-1, 2*i2];

eta_lo = -1; eta_hi = 1;
for iter = 1:50
    eta_mid = (eta_lo + eta_hi) / 2;
    N1 = (1-eta_mid)/2; N2 = (1+eta_mid)/2;
    xp = N1*x1 + N2*x2;
    yp = N1*y1 + N2*y2;
    
    px = dot([xp, yp], ext_dir);
    py = dot([xp, yp], crack_normal);
    th_mid = atan2(py, px);
    
    if abs(th_mid) < 1e-12
        break;
    end
    
    % Determine which half contains theta=0
    N1_lo = (1-eta_lo)/2; N2_lo = (1+eta_lo)/2;
    xp_lo = N1_lo*x1 + N2_lo*x2;
    yp_lo = N1_lo*y1 + N2_lo*y2;
    th_lo = atan2(dot([xp_lo,yp_lo], crack_normal), dot([xp_lo,yp_lo], ext_dir));
    
    if th_lo * th_mid < 0
        eta_hi = eta_mid;
    else
        eta_lo = eta_mid;
    end
end

% Evaluate stress mode at the found eta (theta ≈ 0)
eta_star = (eta_lo + eta_hi) / 2;
N1 = (1-eta_star)/2; N2 = (1+eta_star)/2;
dN1 = -1/2; dN2 = 1/2;

xp = N1*x1 + N2*x2;
yp = N1*y1 + N2*y2;
xd = dN1*x1 + dN2*x2;
yd = dN1*y1 + dN2*y2;

Jdet = xp*yd - xd*yp;
absJ = abs(Jdet);

if absJ < 1e-15
    return;
end

b1 = (1/absJ)*[yd, 0; 0, -xd; -xd, yd];
b2 = (1/absJ)*[-yp, 0; 0, xp; xp, -yp];
Nm  = [N1, 0, N2, 0; 0, N1, 0, N2];
dNm = [dN1, 0, dN2, 0; 0, dN1, 0, dN2];
B1e = b1 * Nm;
B2e = b2 * dNm;

% Stress mode (Eq. 31): Psi_s = D[B1*Vu_s*(-Lambda_s) + B2*Vu_s]
Vu_s_e = Vu_s(dof_e, :);
Psi = D * (B1e * Vu_s_e * (-Lam_s) + B2e * Vu_s_e);

% Singular stress at theta=0, xi=1 (boundary): sigma_s = Psi * c_s
stress_s = real(Psi * c_s);  % [sigma_xx; sigma_yy; sigma_xy] in global

% Distance from crack tip to this point
L0 = norm([xp, yp]);

% Transform stress to local crack coordinates
ca = ext_dir(1); sa = ext_dir(2);
T_s = [ca^2, sa^2, 2*ca*sa;
       sa^2, ca^2, -2*ca*sa;
       -ca*sa, ca*sa, ca^2-sa^2];
stress_local = T_s * stress_s;

% Eq. 32: {KI; KII} = sqrt(2*pi*L0) * {sigma_yy_local; sigma_xy_local}
KI  = sqrt(2*pi*L0) * stress_local(2);
KII = sqrt(2*pi*L0) * stress_local(3);

end
