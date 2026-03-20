%% main_inclined_crack.m
% =======================================================================
% Reproduces Inclined centre crack under uniaxial tension
% Yaoxin Huang
% =======================================================================
clear; clc; close all;
addpath(pwd);

fprintf('=== XSBFEM: Inclined Centre Crack Under Uniaxial Tension ===\n\n');

%% Problem Parameters (Section 5.6)
b     = 16.0;        % Plate dimension (m)
a     = 0.5;         % Half crack length (m)
sigma = 1.0e-3;      % Applied stress (MPa) = 1 kPa
E_mod = 1.0;         % Young's modulus (MPa)
nu    = 0.2;         % Poisson's ratio
plane_type = 'strain';

G     = E_mod / (2*(1+nu));
kappa = 3 - 4*nu;
D     = get_material_matrix(E_mod, nu, plane_type);

angles_deg = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180];

fprintf('G = %.4f, kappa = %.1f\n', G, kappa);

%% Diagnostic: Pure Mode I (straight crack)
fprintf('\n--- Diagnostic: Pure Mode I ---\n');

R_test = 0.3; n_test = 48;
eps_a = pi / n_test;
theta_test = linspace(-(pi - eps_a), (pi - eps_a), n_test)';
bnd_test = R_test * [cos(theta_test), sin(theta_test)];

ub_test = williams_disp(theta_test, R_test, 1.0, 0.0, G, kappa);
[~, KI_t, KII_t, dt] = sbfem_crack_tip(bnd_test, [0,0], D, ub_test, E_mod, nu, plane_type, 3);

fprintf('  KI:  exact=1.0000, disp=%.6f (%.2f%%), stress=%.6f (%.2f%%)\n', ...
    dt.KI_disp, abs(dt.KI_disp-1)*100, dt.KI_stress, abs(dt.KI_stress-1)*100);
fprintf('  KII: exact=0.0000, disp=%.6f,              stress=%.6f\n', ...
    dt.KII_disp, dt.KII_stress);

%% Diagnostic: Pure Mode II
fprintf('\n--- Diagnostic: Pure Mode II ---\n');
ub_test2 = williams_disp(theta_test, R_test, 0.0, 1.0, G, kappa);
[~, ~, ~, dt2] = sbfem_crack_tip(bnd_test, [0,0], D, ub_test2, E_mod, nu, plane_type, 3);

fprintf('  KI:  exact=0.0000, disp=%.6f,              stress=%.6f\n', ...
    dt2.KI_disp, dt2.KI_stress);
fprintf('  KII: exact=1.0000, disp=%.6f (%.2f%%), stress=%.6f (%.2f%%)\n', ...
    dt2.KII_disp, abs(dt2.KII_disp-1)*100, dt2.KII_stress, abs(dt2.KII_stress-1)*100);

%% Part A: Direct SBFEM Validation for all angles
fprintf('\n--- Part A: Direct SBFEM for all angles ---\n');

R_boundary = 0.3;
n_boundary = 48;

n_angles = length(angles_deg);
KI_analytical   = zeros(n_angles, 1);
KII_analytical  = zeros(n_angles, 1);
KI_disp   = zeros(n_angles, 1);
KII_disp  = zeros(n_angles, 1);
KI_stress  = zeros(n_angles, 1);
KII_stress = zeros(n_angles, 1);

for ia = 1:n_angles
    alpha_deg = angles_deg(ia);
    alpha = alpha_deg * pi / 180;
    
    % Analytical SIFs (Eq. 51)
    KI_exact  = sigma * sqrt(pi*a) * cos(alpha)^2;
    KII_exact = sigma * sqrt(pi*a) * cos(alpha) * sin(alpha);
    
    KI_analytical(ia)  = KI_exact;
    KII_analytical(ia) = KII_exact;
    
    % Right crack tip
    tip = [a*cos(alpha), a*sin(alpha)];
    
    % Create circular boundary in local crack coordinates, transform to global
    eps_angle = pi / n_boundary;
    theta_local = linspace(-(pi - eps_angle), (pi - eps_angle), n_boundary)';
    
    cos_a = cos(alpha); sin_a = sin(alpha);
    R_mat = [cos_a, -sin_a; sin_a, cos_a];
    
    boundary_coords = zeros(n_boundary, 2);
    for i = 1:n_boundary
        pt = R_mat * [R_boundary*cos(theta_local(i)); R_boundary*sin(theta_local(i))];
        boundary_coords(i,:) = tip + pt';
    end
    
    % Williams expansion boundary displacement
    ub_exact = zeros(2*n_boundary, 1);
    for i = 1:n_boundary
        th = theta_local(i);
        fac = 1/(2*G) * sqrt(R_boundary/(2*pi));
        c2 = cos(th/2); s2 = sin(th/2);
        
        ux_l = KI_exact*fac*c2*(kappa-1+2*s2^2) + KII_exact*fac*s2*(kappa+1+2*c2^2);
        uy_l = KI_exact*fac*s2*(kappa+1-2*c2^2) - KII_exact*fac*c2*(kappa-1-2*s2^2);
        
        u_g = R_mat * [ux_l; uy_l];
        ub_exact(2*i-1) = u_g(1);
        ub_exact(2*i)   = u_g(2);
    end
    
    % SBFEM computation
    [~, ~, ~, se_data] = sbfem_crack_tip(boundary_coords, tip, D, ub_exact, ...
                                          E_mod, nu, plane_type, 3);
    
    KI_disp(ia)    = se_data.KI_disp;
    KII_disp(ia)   = se_data.KII_disp;
    KI_stress(ia)  = se_data.KI_stress;
    KII_stress(ia) = se_data.KII_stress;
    
    fprintf('  theta=%2d: KI [exact %.6f | disp %.6f (%.1f%%) | stress %.6f (%.1f%%)]\n', ...
        alpha_deg, KI_exact, KI_disp(ia), ...
        abs(KI_disp(ia)-KI_exact)/max(abs(KI_exact),1e-15)*100, ...
        KI_stress(ia), abs(KI_stress(ia)-KI_exact)/max(abs(KI_exact),1e-15)*100);
    fprintf('           KII[exact %.6f | disp %.6f (%.1f%%) | stress %.6f (%.1f%%)]\n', ...
        KII_exact, KII_disp(ia), ...
        abs(KII_disp(ia)-KII_exact)/max(abs(KII_exact),1e-15)*100, ...
        KII_stress(ia), abs(KII_stress(ia)-KII_exact)/max(abs(KII_exact),1e-15)*100);
end

%% Plotting
fprintf('\n--- Generating Plots (Figure 18) ---\n');

fig = figure('Position', [100 100 1200 500]);

% --- Mode I SIF ---
subplot(1,2,1);
plot(angles_deg, KI_analytical, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Analytical result');
hold on;
plot(angles_deg, KI_stress, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Stress based method');
plot(angles_deg, KI_disp, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Displacement based Method');
hold off;
xlabel('Angle \theta / (°)', 'FontSize', 12);
ylabel('Mode I SIF / kPa \cdot m^{1/2}', 'FontSize', 12);
title('(a) Mode I SIF', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10);
grid on; xlim([0 90]);
set(gca, 'FontSize', 11);

% --- Mode II SIF ---
subplot(1,2,2);
plot(angles_deg, KII_analytical, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Analytical result');
hold on;
plot(angles_deg, KII_stress, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Stress based method');
plot(angles_deg, KII_disp, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 8, ...
    'DisplayName', 'Displacement based method');
hold off;
xlabel('Angle \theta / (°)', 'FontSize', 12);
ylabel('Mode II SIF / kPa \cdot m^{1/2}', 'FontSize', 12);
title('(b) Mode II SIF', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10);
grid on; xlim([0 90]);
ylim_max = max(abs(KII_analytical)) * 1.3;
ylim([0, ylim_max]);
set(gca, 'FontSize', 11);

sgtitle('SIFs for Inclined Centre Crack Under Different Angles', ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(fig, 'sif_comparison.png');
fprintf('Figure saved.\n');

%% Summary table
fprintf('\n=== Summary ===\n');
fprintf('%5s | %10s %10s %7s %10s %7s | %10s %10s %7s %10s %7s\n', ...
    'Angle', 'KI_exact', 'KI_disp', 'err%', 'KI_str', 'err%', ...
    'KII_exact', 'KII_disp', 'err%', 'KII_str', 'err%');
fprintf('%s\n', repmat('-', 1, 115));
for ia = 1:n_angles
    eI_d  = abs(KI_disp(ia)-KI_analytical(ia))/max(abs(KI_analytical(ia)),1e-15)*100;
    eI_s  = abs(KI_stress(ia)-KI_analytical(ia))/max(abs(KI_analytical(ia)),1e-15)*100;
    eII_d = abs(KII_disp(ia)-KII_analytical(ia))/max(abs(KII_analytical(ia)),1e-15)*100;
    eII_s = abs(KII_stress(ia)-KII_analytical(ia))/max(abs(KII_analytical(ia)),1e-15)*100;
    fprintf('%5d | %10.6f %10.6f %6.2f %10.6f %6.2f | %10.6f %10.6f %6.2f %10.6f %6.2f\n', ...
        angles_deg(ia), KI_analytical(ia), KI_disp(ia), eI_d, KI_stress(ia), eI_s, ...
        KII_analytical(ia), KII_disp(ia), eII_d, KII_stress(ia), eII_s);
end

fprintf('\nDone.\n');

%% ====================================================================
function ub = williams_disp(theta, R, KI, KII, G, kappa)
% Williams expansion displacement at r=R
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
