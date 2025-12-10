% SP_DRCHLT_L2_PROJ_ST: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj_st (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    object defining the space of discrete functions (see sp_vector)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%  h:     function handle to compute the Dirichlet condition
%  sides: boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
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


function [u, dofs_st] = sp_drchlt_l2_proj_st (sp, spt, msh, msht, h, sides, isitweak)
  arguments
    sp
    spt
    msh
    msht
    h
    sides
    isitweak = 'no'
  end

  rhs  = zeros (sp.ndof, spt.ndof);

  dofs = [];

  nent = 0;
  for iside = sides
    nent = nent + msh.boundary(iside).nel * sp.boundary(iside).nsh_max^2;
    dofs = union (dofs, sp.boundary(iside).dofs);
  end
  
  dofs_st = [];
  switch isitweak
      case 'no'
          start = 2;  % 2:spt.ndof if the initial data are strongly imposed
      case 'yes'
          start = 1;  % 1:spt.ndof if the initial data are weakly imposed
  end

  for i = start:spt.ndof
    dofs_st = cat(1, dofs_st, dofs + (i-1)*(sp.ndof));
  end

  rows = zeros (nent, 1);
  cols = zeros (nent, 1);
  vals = zeros (nent, 1);
  
% then: kron( dofs.time, dofs ) which are the drchlet dofs in space-time

  ncounter = 0;
  for iside = sides
% Restrict the function handle to the specified side, in any dimension, hside = @(x,y) h(x,y,iside)
    hside = @(varargin) h(varargin{:},iside);
    [rs,cs,vs] = op_u_v_tp (sp.boundary(iside), sp.boundary(iside), msh.boundary(iside));

    bnd_dofs = sp.boundary(iside).dofs;
    
    rows(ncounter+(1:numel(rs))) = bnd_dofs(rs);
    cols(ncounter+(1:numel(rs))) = bnd_dofs(cs);
    vals(ncounter+(1:numel(rs))) = vs;
    ncounter = ncounter + numel (rs);

    rhs(bnd_dofs,:) = rhs(bnd_dofs,:) + ...
         reshape( op_f_v_st_tp (sp.boundary(iside), spt, msh.boundary(iside), msht, hside), ...
                    numel(bnd_dofs), []);
  end

  rhs = rhs(:);

  Ms = sparse (rows(1:ncounter), cols(1:ncounter), vals(1:ncounter));
  Mt = op_u_v_tp (spt, spt, msht);
  M = kron(Mt,Ms);
  u = M(dofs_st, dofs_st) \ rhs(dofs_st, 1); % could be more efficient
 
end
