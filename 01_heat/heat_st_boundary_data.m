% HEAT_ST_BOUNDARY_DATA: <description>
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
%
function [vel_drchlt, drchlt_dofs, vel_iniz] = ... 
    heat_st_boundary_data(spx, spt, xmsh, tmsh, h, drchlt_sides, isitweak)
  arguments
    spx
    spt
    xmsh
    tmsh
    h
    drchlt_sides
    isitweak = 'no'
  end
% projection of initial data.
dim = numel(spx.degree);
switch dim
 case 1 % 1D in space
  vel0 = op_f_v_tp (spx, xmsh, @(x) h(x, 0)); % valutazione di \int_\Omega u0 * v dx
 case 2 % 2D in spazio
  vel0 = op_f_v_tp (spx, xmsh, @(x,y) h(x,y,0)); % valutazione di \int_\Omega u0 * v dx
 case 3 % 3D in spazio
  vel0 = op_f_v_tp (spx, xmsh, @(x,y,z) h(x,y,z,0)); % valutazione di \int_\Omega u0 * v dx
end
Mx = op_u_v_tp (spx, spx, xmsh); % L2 repr. matrix to project the data
fprintf('Projecting initial data... \n\n')
vel_iniz = Mx\vel0; % L2 projection of initial data

if dim == 1
  % THIS CASE MUST BE FIXED
  fprintf('Projecting Dirichlet boundary data... \n\n')
  fprintf('THIS CASE IS NOT WORKING! \n\n')
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_st (spx, spt, xmsh, tmsh, h, drchlt_sides, isitweak)
else
  % DIRICHLET DATA IN SPACE-TIME
  fprintf('Projecting Dirichlet boundary data... \n\n')
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_st (spx, spt, xmsh, tmsh, h, drchlt_sides, isitweak);
end

end