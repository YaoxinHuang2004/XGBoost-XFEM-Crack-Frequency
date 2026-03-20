function [KI, KII] = sbfem_compute_sif(data, ub, crack_angle, E, nu, plane_type)
% SBFEM_COMPUTE_SIF - Compute stress intensity factors from SBFEM solution
%
% Computes mode I and II SIFs using both displacement-based and stress-based
% methods from the SBFEM solution of a crack tip subdomain.
%
% Input:
%   data        - SBFEM solution data structure (from sbfem_stiffness)
%   ub          - boundary displacement vector (2n x 1)
%   crack_angle - crack inclination angle in radians (angle of crack line w.r.t. x-axis)
%   E, nu       - elastic constants
%   plane_type  - 'strain' or 'stress'
%
% Output:
%   KI, KII - Mode I and II stress intensity factors
%

% Extract SBFEM data
Vu   = data.Vu;
Sn   = data.Sn;
nodes = data.nodes;
sc    = data.sc;
D     = data.D;

n    = size(nodes, 1);
ndof = 2 * n;

%% Step 1: Compute integration constants c = Vu^{-1} * ub (Eq. 22)
c = Vu \ ub;

%% Step 2: Identify singular eigenvalue block
% Eigenvalues of Sn (upper triangular from Schur decomposition)
eigvals = diag(Sn);
% Singular modes: -1 < Re(lambda) < 0 (for r^{-1/2} singularity, Re(lambda) ≈ -0.5)
tol = 0.05;
sing_idx = find(real(eigvals) > -1 + tol & real(eigvals) < 0 - tol);

if isempty(sing_idx)
    warning('No singular eigenvalues found. SIFs set to zero.');
    KI = 0; KII = 0;
    return;
end

%% Method 1: Displacement-based SIF computation (Eq. 26)
% Find crack face nodes on the boundary
% The crack faces are at theta = +/-pi in local crack coordinates
% In global coordinates, the crack direction from the tip is at angle (crack_angle + pi)

% Coordinates relative to scaling center
xy_rel = nodes - repmat(sc, n, 1);

% Transform to local crack coordinate system
% Local x-axis: along crack extension direction (opposite to crack direction)
% The crack comes INTO the tip from direction crack_angle + pi
% The crack front (theta=0) is at direction crack_angle (extending the crack)
% Wait, need to be more careful.
% For the right crack tip: crack extends to the left (toward center), 
% so the crack direction is crack_angle + pi, and theta=0 is at crack_angle
% For the left crack tip: crack extends to the right (toward center),
% so the crack direction is crack_angle, and theta=0 is at crack_angle + pi

% The local coordinate system at the crack tip:
% The x_local axis points along the crack extension (theta=0 direction)
% For simplicity, we'll determine which direction the crack comes from
% by looking at the geometry

% Compute angles of boundary nodes in local crack coordinate system
cos_a = cos(crack_angle);
sin_a = sin(crack_angle);

% Rotation matrix from global to local crack coordinate system
R = [cos_a, sin_a; -sin_a, cos_a];

% Transform boundary coordinates to local crack system
xy_local = (R * xy_rel')';  % n x 2

% Find the two boundary nodes closest to theta = +pi and theta = -pi
% (these are the crack mouth nodes on the upper and lower crack faces)
angles = atan2(xy_local(:,2), xy_local(:,1));

% Find nodes near theta = pi (upper crack face) and theta = -pi (lower crack face)
% These nodes should be very close in position but on opposite sides of the crack
upper_face_idx = [];
lower_face_idx = [];

for i = 1:n
    if abs(angles(i) - pi) < 0.15 || abs(angles(i) + pi) < 0.15
        if xy_local(i,2) >= 0
            upper_face_idx = [upper_face_idx, i];
        else
            lower_face_idx = [lower_face_idx, i];
        end
    end
end

% Shear modulus and kappa
G = E / (2*(1+nu));
if strcmpi(plane_type, 'strain')
    kappa = 3 - 4*nu;
else
    kappa = (3 - nu)/(1 + nu);
end

% Try displacement-based method if crack face nodes found
if ~isempty(upper_face_idx) && ~isempty(lower_face_idx)
    % Use singular modes for displacement
    Vu_s = Vu(:, sing_idx);
    c_s  = c(sing_idx);
    
    % Singular displacement at boundary (xi=1)
    u_s = real(Vu_s * c_s);  % 2n x 1
    
    % Get displacements at crack face nodes
    % Take the node pair closest to the crack line
    [~, ui] = min(abs(abs(angles(upper_face_idx)) - pi));
    [~, li] = min(abs(abs(angles(lower_face_idx)) - pi));
    
    iup = upper_face_idx(ui);
    ilo = lower_face_idx(li);
    
    % Global displacements at these nodes
    ux_up = u_s(2*iup-1);  uy_up = u_s(2*iup);
    ux_lo = u_s(2*ilo-1);  uy_lo = u_s(2*ilo);
    
    % Transform to local crack coordinates
    u_up_local = R * [ux_up; uy_up];
    u_lo_local = R * [ux_lo; uy_lo];
    
    % Distance from crack tip to crack mouth node
    r0_up = norm(xy_rel(iup,:));
    r0_lo = norm(xy_rel(ilo,:));
    r0 = (r0_up + r0_lo) / 2;
    
    % SIFs from displacement method (Eq. 26)
    factor = G / (kappa + 1) * sqrt(2*pi/r0);
    KI_disp  = factor * (u_up_local(2) - u_lo_local(2));
    KII_disp = factor * (u_up_local(1) - u_lo_local(1));
else
    KI_disp = NaN;
    KII_disp = NaN;
end

%% Method 2: Stress-based SIF computation (Eq. 32)
% Compute stress modes at boundary points and interpolate to theta=0

% Compute stress mode at Gauss points along boundary elements
n_sample = 100;  % Total sample points around boundary
theta_samples = zeros(n_sample, 1);
sigma_yy_s = zeros(n_sample, 1);
sigma_xy_s = zeros(n_sample, 1);
L0_samples = zeros(n_sample, 1);

xy = nodes - repmat(sc, n, 1);
sample_count = 0;

pts_per_elem = max(2, ceil(n_sample / n));

for e = 1:n
    i1 = e;
    i2 = mod(e, n) + 1;
    
    x1 = xy(i1,1); y1 = xy(i1,2);
    x2 = xy(i2,1); y2 = xy(i2,2);
    
    dof_elem = [2*i1-1, 2*i1, 2*i2-1, 2*i2];
    
    for p = 1:pts_per_elem
        eta = -1 + 2*(p-0.5)/pts_per_elem;
        
        % Shape functions
        N1 = (1-eta)/2; N2 = (1+eta)/2;
        dN1 = -1/2;     dN2 = 1/2;
        
        x_pt = N1*x1 + N2*x2;
        y_pt = N1*y1 + N2*y2;
        
        x_d = dN1*x1 + dN2*x2;
        y_d = dN1*y1 + dN2*y2;
        
        Jdet = x_pt*y_d - x_d*y_pt;
        absJ = abs(Jdet);
        
        if absJ < 1e-15, continue; end
        
        % b1 and b2
        b1 = (1/absJ) * [y_d, 0; 0, -x_d; -x_d, y_d];
        b2 = (1/absJ) * [-y_pt, 0; 0, x_pt; x_pt, -y_pt];
        
        Nmat = [N1, 0, N2, 0; 0, N1, 0, N2];
        dNmat = [dN1, 0, dN2, 0; 0, dN1, 0, dN2];
        
        B1_elem = b1 * Nmat;
        B2_elem = b2 * dNmat;
        
        % Extract singular modes
        Vu_s_elem = zeros(4, length(sing_idx));
        for si = 1:length(sing_idx)
            Vu_s_elem(:, si) = Vu(dof_elem, sing_idx(si));
        end
        
        % Stress mode (Eq. 31): Psi_s = D * [B1 * Vu_s * (-Sn_s) + B2 * Vu_s]
        Sn_s = Sn(sing_idx, sing_idx);
        psi_s = D * (B1_elem * Vu_s_elem * (-Sn_s) + B2_elem * Vu_s_elem);
        
        % Compute stress at this point
        c_s = c(sing_idx);
        stress_s = real(psi_s * c_s);
        
        % Transform point to local crack coordinates
        pt_local = R * [x_pt; y_pt];
        theta_pt = atan2(pt_local(2), pt_local(1));
        L0_pt = norm([x_pt, y_pt]);
        
        sample_count = sample_count + 1;
        if sample_count > n_sample
            break;
        end
        theta_samples(sample_count) = theta_pt;
        
        % Transform stress to local crack coordinates
        % sigma_local = T * sigma_global where T is stress transformation
        c2 = cos_a^2; s2 = sin_a^2; cs = cos_a*sin_a;
        T_stress = [c2, s2, 2*cs; s2, c2, -2*cs; -cs, cs, c2-s2];
        stress_local = T_stress * stress_s;
        
        sigma_yy_s(sample_count) = stress_local(2);  % sigma_yy in local
        sigma_xy_s(sample_count) = stress_local(3);  % sigma_xy in local
        L0_samples(sample_count) = L0_pt;
    end
    if sample_count >= n_sample, break; end
end

% Trim arrays
theta_samples = theta_samples(1:sample_count);
sigma_yy_s = sigma_yy_s(1:sample_count);
sigma_xy_s = sigma_xy_s(1:sample_count);
L0_samples = L0_samples(1:sample_count);

% Sort by angle
[theta_sorted, sort_idx] = sort(theta_samples);
sigma_yy_sorted = sigma_yy_s(sort_idx);
sigma_xy_sorted = sigma_xy_s(sort_idx);
L0_sorted = L0_samples(sort_idx);

% Interpolate to theta = 0 (crack front)
% Remove duplicates and NaN
valid = ~isnan(sigma_yy_sorted) & ~isinf(sigma_yy_sorted);
theta_v = theta_sorted(valid);
syy_v = sigma_yy_sorted(valid);
sxy_v = sigma_xy_sorted(valid);
L0_v = L0_sorted(valid);

if length(theta_v) >= 3
    try
        syy_0 = interp1(theta_v, syy_v, 0, 'spline');
        sxy_0 = interp1(theta_v, sxy_v, 0, 'spline');
        L0_0  = interp1(theta_v, L0_v,  0, 'spline');
        
        KI_stress  = sqrt(2*pi*L0_0) * syy_0;
        KII_stress = sqrt(2*pi*L0_0) * sxy_0;
    catch
        KI_stress = NaN;
        KII_stress = NaN;
    end
else
    KI_stress = NaN;
    KII_stress = NaN;
end

%% Choose best result
if ~isnan(KI_stress)
    KI = KI_stress;
    KII = KII_stress;
elseif ~isnan(KI_disp)
    KI = KI_disp;
    KII = KII_disp;
else
    KI = 0;
    KII = 0;
end

end
