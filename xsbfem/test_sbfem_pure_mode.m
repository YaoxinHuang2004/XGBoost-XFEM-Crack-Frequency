%% test_sbfem_pure_mode.m
% Validate SBFEM SIF computation for pure mode I, II, and mixed mode
clear; clc; close all;
addpath(pwd);

fprintf('=== SBFEM SIF Validation ===\n\n');

E_mod = 100; nu = 0.3; plane_type = 'strain';
G = E_mod / (2*(1+nu));
kappa = 3 - 4*nu;
D = get_material_matrix(E_mod, nu, plane_type);
a = 0.5;

R = 0.2;        % Boundary radius
n_pts = 48;     % Boundary nodes

eps_a = pi / n_pts;
theta = linspace(-(pi - eps_a), (pi - eps_a), n_pts)';
bnd = R * [cos(theta), sin(theta)];
sc = [0, 0];

fprintf('Parameters: E=%.0f, nu=%.1f, G=%.4f, kappa=%.2f\n', E_mod, nu, G, kappa);
fprintf('Boundary: R=%.2f, n_pts=%d\n\n', R, n_pts);

%% Test 1: Pure Mode I (KI=1, KII=0)
fprintf('--- Test 1: Pure Mode I ---\n');
KI_app = 1.0; KII_app = 0.0;
ub = williams_displacement(theta, R, KI_app, KII_app, G, kappa);

[~, KI_c, KII_c, data] = sbfem_crack_tip(bnd, sc, D, ub, E_mod, nu, plane_type, 3);

print_result('Mode I', KI_app, KII_app, KI_c, KII_c, data);

%% Test 2: Pure Mode II (KI=0, KII=1)
fprintf('\n--- Test 2: Pure Mode II ---\n');
KI_app = 0.0; KII_app = 1.0;
ub = williams_displacement(theta, R, KI_app, KII_app, G, kappa);

[~, KI_c, KII_c, data] = sbfem_crack_tip(bnd, sc, D, ub, E_mod, nu, plane_type, 3);

print_result('Mode II', KI_app, KII_app, KI_c, KII_c, data);

%% Test 3: Mixed mode (KI=0.7, KII=0.5)
fprintf('\n--- Test 3: Mixed Mode ---\n');
KI_app = 0.7; KII_app = 0.5;
ub = williams_displacement(theta, R, KI_app, KII_app, G, kappa);

[~, KI_c, KII_c, data] = sbfem_crack_tip(bnd, sc, D, ub, E_mod, nu, plane_type, 3);

print_result('Mixed', KI_app, KII_app, KI_c, KII_c, data);

%% Test 4: Inclined crack (alpha=45 deg)
fprintf('\n--- Test 4: Inclined crack (45 deg) ---\n');
alpha = pi/4;
sigma_test = 1.0;
KI_exact = sigma_test * sqrt(pi*a) * cos(alpha)^2;
KII_exact = sigma_test * sqrt(pi*a) * cos(alpha) * sin(alpha);

tip = [a*cos(alpha), a*sin(alpha)];
Rmat = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];

bnd_g = zeros(n_pts, 2);
for i = 1:n_pts
    pt = Rmat * [R*cos(theta(i)); R*sin(theta(i))];
    bnd_g(i,:) = tip + pt';
end

ub_g = zeros(2*n_pts, 1);
for i = 1:n_pts
    th = theta(i); fac = 1/(2*G) * sqrt(R/(2*pi));
    c2 = cos(th/2); s2 = sin(th/2);
    ux_l = KI_exact*fac*c2*(kappa-1+2*s2^2) + KII_exact*fac*s2*(kappa+1+2*c2^2);
    uy_l = KI_exact*fac*s2*(kappa+1-2*c2^2) - KII_exact*fac*c2*(kappa-1-2*s2^2);
    u_g = Rmat * [ux_l; uy_l];
    ub_g(2*i-1) = u_g(1); ub_g(2*i) = u_g(2);
end

[~, KI_c, KII_c, data] = sbfem_crack_tip(bnd_g, tip, D, ub_g, E_mod, nu, plane_type, 3);

fprintf('  Exact:    KI=%.6f, KII=%.6f\n', KI_exact, KII_exact);
fprintf('  Computed: KI=%.6f, KII=%.6f\n', KI_c, KII_c);
fprintf('  Error:    KI=%.2f%%, KII=%.2f%%\n', ...
    abs(KI_c-KI_exact)/max(abs(KI_exact),1e-15)*100, ...
    abs(KII_c-KII_exact)/max(abs(KII_exact),1e-15)*100);

%% Test 5: Convergence study
fprintf('\n--- Test 5: Convergence (Mode I, KI=1) ---\n');
n_list = [16, 24, 32, 48, 64, 80, 96];
KI_conv = zeros(length(n_list),1);

for ip = 1:length(n_list)
    np = n_list(ip);
    ea = pi/np;
    th_c = linspace(-(pi-ea), (pi-ea), np)';
    bnd_c = R * [cos(th_c), sin(th_c)];
    ub_c = williams_displacement(th_c, R, 1.0, 0.0, G, kappa);
    [~, ki, ~, ~] = sbfem_crack_tip(bnd_c, [0,0], D, ub_c, E_mod, nu, plane_type, 3);
    KI_conv(ip) = ki;
    fprintf('  n=%3d: KI=%.8f, err=%.4f%%\n', np, ki, abs(ki-1)*100);
end

%% Plot convergence
figure;
semilogy(n_list, abs(KI_conv - 1)*100, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Number of boundary nodes'); ylabel('Relative error (%)');
title('SBFEM Mode I SIF convergence'); grid on;
saveas(gcf, 'test_convergence.png');

fprintf('\nAll tests completed.\n');

%% ====================================================================
function ub = williams_displacement(theta, R, KI, KII, G, kappa)
% Williams expansion displacement at r=R for given theta angles
n = length(theta);
ub = zeros(2*n, 1);
fac = 1/(2*G) * sqrt(R/(2*pi));
for i = 1:n
    th = theta(i);
    c2 = cos(th/2); s2 = sin(th/2);
    ub(2*i-1) = KI*fac*c2*(kappa-1+2*s2^2) + KII*fac*s2*(kappa+1+2*c2^2);
    ub(2*i)   = KI*fac*s2*(kappa+1-2*c2^2) - KII*fac*c2*(kappa-1-2*s2^2);
end
end

function print_result(name, KI_app, KII_app, KI_c, KII_c, data)
fprintf('  Applied:  KI=%.4f, KII=%.4f\n', KI_app, KII_app);
fprintf('  Computed: KI=%.6f, KII=%.6f\n', KI_c, KII_c);
if KI_app ~= 0
    fprintf('  KI error: %.4f%%\n', abs(KI_c-KI_app)/abs(KI_app)*100);
else
    fprintf('  KI (abs): %.6f\n', abs(KI_c));
end
if KII_app ~= 0
    fprintf('  KII error: %.4f%%\n', abs(KII_c-KII_app)/abs(KII_app)*100);
else
    fprintf('  KII (abs): %.6f\n', abs(KII_c));
end

% Diagnostics
fprintf('  Singular eigenvalues (%d found):', length(data.sing_idx));
for j = 1:length(data.sing_eigvals)
    fprintf(' %.4f+%.4fi', real(data.sing_eigvals(j)), imag(data.sing_eigvals(j)));
end
fprintf('\n');

% Show first 8 bounded eigenvalues
ev = data.eigvals_bounded;
fprintf('  First 8 eigenvalues (bounded): ');
for j = 1:min(8, length(ev))
    fprintf('%.3f ', real(ev(j)));
end
fprintf('...\n');

if isfield(data, 'KI_stress') && ~isnan(data.KI_stress)
    fprintf('  Stress method: KI=%.6f, KII=%.6f\n', data.KI_stress, data.KII_stress);
end
end
