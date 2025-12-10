% NAME_OF_FUNCTION: <description>
%
%   CALL: 
%
%  INPUT:
%
% OUTPUT:
%
%
% ProjectName - STIGA
% Copyright (C) 2025 Alen Kushova, Gabriele Loli
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See <https://www.gnu.org/licenses/> for more details.

function [geo, msh, space] = parametric_domain (problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

rdim  = numel(trial_degree); % dimention of the st_problem d+1.
switch rdim 
    case 2
        geo   = cell(1,rdim+1); 
        msh   = cell(1,rdim+1); 
        space = cell(2,rdim+1);
        geo{1}  = geo_load('geo_square.txt');
    case 3
        geo   = cell(1,rdim+1); 
        msh   = cell(1,rdim+1); 
        space = cell(2,rdim+1);
        geo{1}  = geo_load('geo_cube.txt');
    case 4
        geo   = cell(1,rdim+1); 
        msh   = cell(1,rdim+1); 
        space = cell(2,rdim+1);
        geo{1}= 'No ST geometries for 3D+1D case';
end

if (rdim<4) % here we have the space-time geometries
    % space time knots for trial functions
    [knots, zeta]    = kntrefine(geo{1}.nurbs.knots, nsub-1, trial_degree, trial_regularity);
    knots   = kntunclamp(knots, trial_degree, trial_regularity, []);
    % space time knots for test functions
    [test_knots, ~]    = kntrefine(geo{1}.nurbs.knots, nsub-1, test_degree, test_regularity);
    test_knots   = kntunclamp(test_knots, test_degree, test_regularity, []);
    % define quadrature rules 
    rule     = msh_gauss_nodes (nquad);
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh{1}   = msh_cartesian (zeta, qn, qw, geo{1});
    space{1,1} = sp_bspline (knots, trial_degree, msh{1});
    space{2,1} = sp_bspline (test_knots, test_degree, msh{1});
end


for i = 2:rdim+1
    j = i-1;
    geo{i}   = geo_load(nrbline ([0 0], [1 0]));
    [knots, zeta] = kntrefine(geo{i}.nurbs.knots, nsub(j)-1, trial_degree(j), trial_regularity(j));
    knots   = kntunclamp(knots, trial_degree(j), trial_regularity(j), []);
    [test_knots, ~] = kntrefine(geo{i}.nurbs.knots, nsub(j)-1, test_degree(j), test_regularity(j));
    test_knots = kntunclamp(test_knots, test_degree(j), test_regularity(j), []);
    rule     = msh_gauss_nodes (nquad(j));
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh{i}   = msh_cartesian (zeta, qn, qw, geo{i});
    space{1,i} = sp_bspline (knots, trial_degree(j), msh{i}); 
    space{2,i} = sp_bspline (test_knots, test_degree(j), msh{i}); 
end
end