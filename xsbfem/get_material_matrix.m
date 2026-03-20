function D = get_material_matrix(E, nu, plane_type)
% GET_MATERIAL_MATRIX - Compute 2D plane stress/strain constitutive matrix
%
% Input:
%   E          - Young's modulus
%   nu         - Poisson's ratio
%   plane_type - 'strain' or 'stress'
%
% Output:
%   D - 3x3 constitutive matrix [sigma_xx; sigma_yy; sigma_xy] = D * [eps_xx; eps_yy; gamma_xy]

if strcmpi(plane_type, 'strain')
    factor = E / ((1 + nu) * (1 - 2*nu));
    D = factor * [1-nu,  nu,    0;
                  nu,    1-nu,  0;
                  0,     0,    (1-2*nu)/2];
elseif strcmpi(plane_type, 'stress')
    factor = E / (1 - nu^2);
    D = factor * [1,   nu,  0;
                  nu,  1,   0;
                  0,   0,  (1-nu)/2];
else
    error('plane_type must be ''strain'' or ''stress''');
end
end
