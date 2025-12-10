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
function [u_drchlt, drchlt_dofs, u_iniz] = ... 
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
dim = xmsh.ndim;
switch dim
 case 1 % 1D in space
  u0 = op_f_v_tp (spx, xmsh, @(x) h(x, 0)); % valutazione di \int_\Omega u0 * v dx
 case 2 % 2D in spazio
  u0 = op_f_v_tp (spx, xmsh, @(x,y) h(x,y,0)); % valutazione di \int_\Omega u0 * v dx
 case 3 % 3D in spazio
  u0 = op_f_v_tp (spx, xmsh, @(x,y,z) h(x,y,z,0)); % valutazione di \int_\Omega u0 * v dx
end
Mx = op_u_v_tp (spx, spx, xmsh); % L2 repr. matrix to project the data
fprintf('Projecting initial data... \n\n')
u_iniz = Mx\u0; % L2 projection of initial data

if dim == 1
  % The 1d space domain hase 0d-boundaries. We treat it separately.  
  fprintf('Projecting Dirichlet boundary data... \n\n')
  Mt = op_u_v_tp (spt, spt, tmsh); % L2 repr. matrix to project the data
  u_l_t = op_f_v_tp (spt, tmsh, @(t) h(xmsh.map(xmsh.breaks{:}(1)), t));
  u_r_t = op_f_v_tp (spt, tmsh, @(t) h(xmsh.map(xmsh.breaks{:}(end)), t));
  % Compute the projections on {a} x [0,T] and {b} x [0,T] 
  % where [a,b] is the space interval domain.
  u_l_drchlet = reshape(Mt\u_l_t, 1, []); % L2 projection of dirichlet data at {a} x [0,T]
  u_r_drchlet = reshape(Mt\u_r_t, 1, []); % L2 projection of dirichlet data at {b} x [0,T]
  u_drchlt = cat(1, u_l_drchlet, u_r_drchlet);
  u_drchlt = u_drchlt(:);

  % Find the corresponding dofs in space-time
  dofs_st = [];
  switch isitweak
      case 'no'
          start = 2;  % 2:spt.ndof if the initial data are strongly imposed
          % and drop the first two dofs at the initial boundaries.
          u_drchlt = u_drchlt(3:end);
      case 'yes'
          start = 1;  % 1:spt.ndof if the initial data are weakly imposed
  end
  % first take boundary dofs in space
  dofs = [spx.boundary(1).dofs; spx.boundary(2).dofs]; 
  % then construct the corresponding dofs in space-time
  for i = start:spt.ndof
    dofs_st = cat(1, dofs_st, dofs + (i-1)*(spx.ndof));
  end 
  drchlt_dofs = dofs_st;
else
  % DIRICHLET DATA IN SPACE-TIME FOR DIM => 2
  fprintf('Projecting Dirichlet boundary data... \n\n')
  [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_st (spx, spt, xmsh, tmsh, h, drchlt_sides, isitweak);
end

end