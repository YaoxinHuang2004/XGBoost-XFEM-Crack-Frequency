function [elem_type, phi, psi, crack_tips, super_elem_info] = ...
    classify_elements(coords, connect, crack_center, crack_half_len, crack_angle, n_layers)
% CLASSIFY_ELEMENTS - Classify elements by their relationship to the crack
%
% Uses level set functions to classify elements as:
%   0 - Standard FE element (not affected by crack)
%   1 - Heaviside enriched (completely cut by crack)
%   2 - Crack tip super element (contains or surrounds crack tip)
%
% Input:
%   coords          - node coordinates (nnodes x 2)
%   connect         - element connectivity (nel x 4)
%   crack_center    - center of crack [x, y]
%   crack_half_len  - half length of crack
%   crack_angle     - inclination angle in radians
%   n_layers        - number of element layers around crack tip for super element
%
% Output:
%   elem_type       - classification of each element (nel x 1)
%   phi             - normal level set at nodes (nnodes x 1)
%   psi             - tangential level set at nodes (nnodes x 1)
%   crack_tips      - [tip1_x, tip1_y; tip2_x, tip2_y]
%   super_elem_info - structure with super element information

if nargin < 6
    n_layers = 3;
end

nnodes = size(coords, 1);
nel = size(connect, 1);
a = crack_half_len;
theta = crack_angle;

% Crack tip coordinates
tip1 = crack_center + a * [cos(theta), sin(theta)];   % Right tip
tip2 = crack_center - a * [cos(theta), sin(theta)];   % Left tip
crack_tips = [tip1; tip2];

% Crack direction and normal vectors
n_crack = [cos(theta), sin(theta)];    % Along crack
n_normal = [-sin(theta), cos(theta)];  % Perpendicular to crack

% Compute level set functions at all nodes
% phi: signed distance to crack line (normal direction)
% psi: signed distance along crack direction (positive outside crack tips)
phi = zeros(nnodes, 1);
psi = zeros(nnodes, 1);

for i = 1:nnodes
    dx = coords(i,:) - crack_center;
    phi(i) = dot(dx, n_normal);            % Normal distance
    tang = dot(dx, n_crack);               % Tangential projection
    psi(i) = abs(tang) - a;                % Positive outside crack tips
end

% Step 1: Find elements containing crack tips
elem_type = zeros(nel, 1);

% For each element, check if it contains a crack tip
tip_elements = zeros(2, 1);  % Element containing each tip

for tip_id = 1:2
    tip_pos = crack_tips(tip_id, :);
    min_dist = inf;
    
    for e = 1:nel
        enodes = connect(e, :);
        elem_coords = coords(enodes, :);
        
        % Check if point is inside the quad element
        if point_in_quad(tip_pos, elem_coords)
            tip_elements(tip_id) = e;
            break;
        end
        
        % Fallback: find element with centroid closest to tip
        centroid = mean(elem_coords);
        d = norm(centroid - tip_pos);
        if d < min_dist
            min_dist = d;
            tip_elements(tip_id) = e;
        end
    end
end

% Step 2: Build super elements around crack tips
% For each tip, expand n_layers layers of elements
super_elem_sets = cell(2, 1);

for tip_id = 1:2
    % Start with the tip element
    seed = tip_elements(tip_id);
    layer_elems = seed;
    all_elems = seed;
    
    for lay = 1:n_layers-1
        % Find all nodes in current layer elements
        layer_nodes = unique(connect(layer_elems, :));
        
        % Find all elements sharing these nodes
        new_elems = [];
        for e = 1:nel
            if any(ismember(connect(e,:), layer_nodes)) && ~ismember(e, all_elems)
                new_elems = [new_elems; e];
            end
        end
        
        layer_elems = new_elems;
        all_elems = [all_elems; new_elems];
    end
    
    super_elem_sets{tip_id} = unique(all_elems);
end

% Mark super element members
for tip_id = 1:2
    elem_type(super_elem_sets{tip_id}) = 2;
end

% Step 3: Classify remaining elements
for e = 1:nel
    if elem_type(e) == 2
        continue;  % Already classified as super element
    end
    
    enodes = connect(e, :);
    phi_e = phi(enodes);
    psi_e = psi(enodes);
    
    % Check if element is cut by crack
    phi_cross = (max(phi_e) > 0 && min(phi_e) < 0);  % Crack line crosses element
    psi_inside = any(psi_e < 0);  % Some part is between crack tips
    
    if phi_cross && psi_inside
        elem_type(e) = 1;  % Heaviside enriched
    end
end

% Step 4: Build super element boundary information
super_elem_info = struct();

for tip_id = 1:2
    se_elems = super_elem_sets{tip_id};
    
    % Collect all nodes in the super element
    se_nodes_all = unique(connect(se_elems, :));
    
    % Find boundary nodes (nodes on the outer boundary of the super element)
    % A boundary node is shared with elements outside the super element
    boundary_nodes = [];
    interior_nodes = [];
    
    for ni = 1:length(se_nodes_all)
        nd = se_nodes_all(ni);
        % Find all elements containing this node
        elems_with_node = find(any(connect == nd, 2));
        % Check if any of these elements are outside the super element
        if any(~ismember(elems_with_node, se_elems))
            boundary_nodes = [boundary_nodes; nd];
        else
            interior_nodes = [interior_nodes; nd];
        end
    end
    
    % Order boundary nodes counterclockwise around the crack tip
    tip_pos = crack_tips(tip_id, :);
    angles = atan2(coords(boundary_nodes,2) - tip_pos(2), ...
                   coords(boundary_nodes,1) - tip_pos(1));
    [~, order] = sort(angles);
    boundary_nodes = boundary_nodes(order);
    
    % Find crack mouth nodes (boundary nodes on the crack face)
    % These are boundary nodes where phi ≈ 0 and psi < 0
    crack_mouth_idx = [];
    for ni = 1:length(boundary_nodes)
        nd = boundary_nodes(ni);
        if abs(phi(nd)) < max(abs(phi))*0.1 && psi(nd) < 0
            crack_mouth_idx = [crack_mouth_idx; ni];
        end
    end
    
    super_elem_info(tip_id).elements = se_elems;
    super_elem_info(tip_id).boundary_nodes = boundary_nodes;
    super_elem_info(tip_id).interior_nodes = interior_nodes;
    super_elem_info(tip_id).all_nodes = se_nodes_all;
    super_elem_info(tip_id).tip_position = tip_pos;
    super_elem_info(tip_id).crack_mouth_idx = crack_mouth_idx;
end

end

function inside = point_in_quad(pt, quad_coords)
% Check if a point is inside a convex quadrilateral using cross products
inside = true;
n = size(quad_coords, 1);
for i = 1:n
    j = mod(i, n) + 1;
    edge = quad_coords(j,:) - quad_coords(i,:);
    to_pt = pt - quad_coords(i,:);
    cross_z = edge(1)*to_pt(2) - edge(2)*to_pt(1);
    if cross_z < -1e-10
        inside = false;
        return;
    end
end
end
