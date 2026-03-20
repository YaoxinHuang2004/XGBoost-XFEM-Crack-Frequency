function [ke_enriched, phys_dofs, virt_dofs] = ...
    heaviside_element_stiffness(elem_coords, phi_elem, D)
% HEAVISIDE_ELEMENT_STIFFNESS - Stiffness of element completely cut by crack
%
% Splits a quad element into two polygons along the crack and computes
% the enriched stiffness matrix using SBFEM for each subdomain.
%
% Input:
%   elem_coords - 4 x 2 nodal coordinates (CCW ordering)
%   phi_elem    - 4 x 1 level set values at element nodes
%   D           - 3 x 3 material constitutive matrix
%
% Output:
%   ke_enriched - 12 x 12 enriched stiffness matrix
%                 (8 physical DOFs + 4 virtual DOFs)
%   phys_dofs   - indices of physical DOFs in the element
%   virt_dofs   - indices of virtual DOFs (at intersection points)

% Find intersection points of crack with element edges
n = 4;
intersections = [];
inter_edges = [];

for i = 1:n
    j = mod(i, n) + 1;
    if phi_elem(i) * phi_elem(j) < 0
        % Edge is crossed by crack
        t = phi_elem(i) / (phi_elem(i) - phi_elem(j));
        inter_pt = elem_coords(i,:) + t * (elem_coords(j,:) - elem_coords(i,:));
        intersections = [intersections; inter_pt];
        inter_edges = [inter_edges; i, j];
    end
end

if size(intersections, 1) ~= 2
    % Element is not cleanly cut - fall back to standard Q4
    ke_enriched = zeros(12, 12);
    ke_q4 = q4_stiffness(elem_coords, D);
    ke_enriched(1:8, 1:8) = ke_q4;
    phys_dofs = 1:8;
    virt_dofs = 9:12;
    return;
end

% Identify nodes on each side of the crack
pos_nodes = find(phi_elem > 0);  % Nodes with phi > 0
neg_nodes = find(phi_elem < 0);  % Nodes with phi < 0

% Build two polygons (subdomains)
% Subdomain 1: nodes with phi > 0 + intersection points
% Subdomain 2: nodes with phi < 0 + intersection points

% Order polygon nodes CCW
poly1_nodes = order_polygon_ccw([elem_coords(pos_nodes,:); intersections]);
poly2_nodes = order_polygon_ccw([elem_coords(neg_nodes,:); intersections]);

% Compute SBFEM stiffness for each subdomain
sc1 = mean(poly1_nodes);
sc2 = mean(poly2_nodes);

[K_sub1, ~] = sbfem_stiffness(poly1_nodes, sc1, D);
[K_sub2, ~] = sbfem_stiffness(poly2_nodes, sc2, D);

% Map subdomain DOFs to element DOFs
% Physical DOFs: 1:8 (2 per physical node, 4 nodes)
% Virtual DOFs: 9:12 (2 per virtual node, 2 intersection points)

% Initialize enriched stiffness
ke_enriched = zeros(12, 12);

% Create DOF mapping for subdomain 1
% Sub1 has: positive nodes + 2 intersection points
n_sub1 = size(poly1_nodes, 1);
sub1_global_dofs = zeros(2*n_sub1, 1);

for i = 1:length(pos_nodes)
    nd = pos_nodes(i);
    % Find this node in poly1
    for pp = 1:n_sub1
        if norm(poly1_nodes(pp,:) - elem_coords(nd,:)) < 1e-10
            sub1_global_dofs(2*pp-1) = 2*nd-1;  % ux
            sub1_global_dofs(2*pp)   = 2*nd;    % uy
            break;
        end
    end
end

% Map intersection points to virtual DOFs (9-12)
virt_count = 0;
for pp = 1:n_sub1
    if sub1_global_dofs(2*pp-1) == 0
        % This is an intersection point
        virt_count = virt_count + 1;
        sub1_global_dofs(2*pp-1) = 8 + 2*virt_count - 1;
        sub1_global_dofs(2*pp)   = 8 + 2*virt_count;
    end
end

% Accumulate sub1 stiffness
for i = 1:2*n_sub1
    for j = 1:2*n_sub1
        gi = sub1_global_dofs(i);
        gj = sub1_global_dofs(j);
        if gi > 0 && gj > 0
            ke_enriched(gi, gj) = ke_enriched(gi, gj) + K_sub1(i, j);
        end
    end
end

% Create DOF mapping for subdomain 2
n_sub2 = size(poly2_nodes, 1);
sub2_global_dofs = zeros(2*n_sub2, 1);

for i = 1:length(neg_nodes)
    nd = neg_nodes(i);
    for pp = 1:n_sub2
        if norm(poly2_nodes(pp,:) - elem_coords(nd,:)) < 1e-10
            sub2_global_dofs(2*pp-1) = 2*nd-1;
            sub2_global_dofs(2*pp)   = 2*nd;
            break;
        end
    end
end

virt_count2 = 0;
for pp = 1:n_sub2
    if sub2_global_dofs(2*pp-1) == 0
        virt_count2 = virt_count2 + 1;
        sub2_global_dofs(2*pp-1) = 8 + 2*virt_count2 - 1;
        sub2_global_dofs(2*pp)   = 8 + 2*virt_count2;
    end
end

% Accumulate sub2 stiffness
for i = 1:2*n_sub2
    for j = 1:2*n_sub2
        gi = sub2_global_dofs(i);
        gj = sub2_global_dofs(j);
        if gi > 0 && gj > 0
            ke_enriched(gi, gj) = ke_enriched(gi, gj) + K_sub2(i, j);
        end
    end
end

phys_dofs = 1:8;
virt_dofs = 9:12;

end

function ordered = order_polygon_ccw(pts)
% Order polygon vertices counterclockwise around centroid
c = mean(pts);
angles = atan2(pts(:,2) - c(2), pts(:,1) - c(1));
[~, idx] = sort(angles);
ordered = pts(idx, :);
end
