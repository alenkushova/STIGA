% ST_SG_ERROR_TP: Evaluate the error in 'Schrödinger Graph' (SG)-norm.
%
%   [errGraph,errl2,errGs] = st_SG_error_tp (spaceS,spaceT,mshS,mshT,u,uex,rhs)
%
% INPUT:
%
%    space:  object defining the space of discrete functions (see sp_scalar)
%    msh:    object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    src:     function handle to evaluate the source term, meaning: 
%             'source f' = 'Schrödinger op.' of the exact solution uex
%
% OUTPUT:
%
%     errGraph:  error in SG norm
%     errl2:  error in L^2 norm
%     errGs: error in SG seminorm (L^2 norm of the residual) 
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

function [errGraph,errl2,errGs] = st_SG_error_tp  (space, msh, u, uex, rhs)

  if (numel(u) ~= space.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  errl2 = 0; errGs = 0;
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', true, 'hessian', true);
    
    [~, err_l2, err_Gs] = st_SG_error(sp_col, msh_col, u, uex, rhs);

    errGs = errGs + err_Gs.^2;
    errl2 = errl2 + err_l2.^2;

  end
  
  errGraph = sqrt (errl2 + errGs);
  errl2 = sqrt (errl2);
  errGs = sqrt (errGs);

end
