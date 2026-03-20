function [coords, connect, nel, nnodes] = create_graded_mesh(Lx, Ly, nx, ny, grade_power)
% CREATE_GRADED_MESH - Create a graded quadrilateral mesh
%
% Creates a structured Q4 mesh on [-Lx/2, Lx/2] x [-Ly/2, Ly/2] with
% node grading that concentrates elements near the center.
%
% Input:
%   Lx, Ly      - domain dimensions
%   nx, ny      - number of elements in x and y directions
%   grade_power - grading power (1=uniform, >1=concentrated at center)
%
% Output:
%   coords  - (nnodes x 2) nodal coordinates
%   connect - (nel x 4) element connectivity (CCW node ordering)
%   nel     - number of elements
%   nnodes  - number of nodes

if nargin < 5
    grade_power = 1;
end

% Generate 1D node positions with grading
% Map t in [-1,1] to x in [-L/2, L/2] with concentration at t=0
t_x = linspace(-1, 1, nx+1)';
t_y = linspace(-1, 1, ny+1)';

if grade_power == 1
    x_nodes = t_x * Lx/2;
    y_nodes = t_y * Ly/2;
else
    x_nodes = sign(t_x) .* abs(t_x).^grade_power * Lx/2;
    y_nodes = sign(t_y) .* abs(t_y).^grade_power * Ly/2;
end

% Generate 2D node coordinates (tensor product)
nnodes = (nx+1) * (ny+1);
coords = zeros(nnodes, 2);

node_id = 0;
for j = 1:ny+1
    for i = 1:nx+1
        node_id = node_id + 1;
        coords(node_id, :) = [x_nodes(i), y_nodes(j)];
    end
end

% Generate element connectivity (CCW ordering)
nel = nx * ny;
connect = zeros(nel, 4);

elem_id = 0;
for j = 1:ny
    for i = 1:nx
        elem_id = elem_id + 1;
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);
        connect(elem_id, :) = [n1, n2, n3, n4];
    end
end

end
