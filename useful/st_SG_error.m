% ST_SG_ERROR_TP: Evaluate the error in 'Schrödinger Graph' (SG)-norm.
%
%   [errGraph,errl2,errGs] = st_SG_error (space,msh,u,uex,rhs)
%
% INPUT:
%
%    spaceST:  structure representing the space of discrete functions (see sp_scalar/sp_evaluate_col)
%    mshST:    structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%    u:        vector of dof weights
%    uex:      function handle to evaluate the exact solution
%    rhs:      function handle to evaluate the source term, meaning: 
%              'source f' = 'Schrödinger op.' of the exact solution uex
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
 
function [errGraph,errl2,errGs,errGraph_e,errl2_e,errGs_e] = st_SG_error(sp,msh,u,uex,rhs)

  grad_valu = sp_eval_msh (u, sp, msh, 'gradient', true);
  grad_valu = reshape (grad_valu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);
  dt_valu   = grad_valu (:, end, :, :);
  
  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  
  valAuex  = reshape (feval (rhs, x{:}), sp.ncomp, 1, msh.nqn, msh.nel);
  
  hessian_valu = sp_eval_msh (u, sp, msh, 'hessian');
  hessian_valu = reshape (hessian_valu, sp.ncomp, msh.rdim, msh.rdim, msh.nqn, msh.nel);
  laplace_valu = zeros(1, 1, msh.nqn, msh.nel); % laplaciano somma delle derivate seconde pure (solo in spazio no tempo)
  
  for r = 1:msh.rdim-1
      laplace_valu = laplace_valu + reshape(hessian_valu(1,r,r,:,:), ...
                          1, 1, msh.nqn, msh.nel);
  end
  
  AuAuh = valAuex -1i*dt_valu + laplace_valu;
  
  valu  = sp_eval_msh (u, sp, msh);
  valu  = reshape (valu, sp.ncomp, msh.nqn, msh.nel);

  valex = reshape (feval (uex, x{:}), sp.ncomp, msh.nqn, msh.nel);
  
  uuh   = valu - valex ; 

  w = msh.quad_weights .* msh.jacdet;

  errl2_e = sum (reshape (( uuh ).*conj( uuh ), [msh.nqn, msh.nel]) .* w); % il modulo in campo complesso
  errGs_e = sum (reshape ((AuAuh).*conj(AuAuh), [msh.nqn, msh.nel]) .* w); % il modulo in campo complesso
  errl2 = sqrt (sum (errl2_e));
  errGs = sqrt (sum (errGs_e));

  errGraph  = sqrt (errl2^2 + errGs^2);

  errGraph_e = sqrt (errl2_e.^2 + errGs_e);
  errGs_e = sqrt (errGs_e);
  
end
