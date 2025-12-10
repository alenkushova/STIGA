% OP_F_V_ST_TP: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), exploiting the tensor product structure.
%
%   rhs = op_f_v_st_tp (space, spaceT, msh, mshT, coeff);
%
% INPUT:
%     
%   spaceS:  object representing the function space (see sp_vector)
%   spaceT: object representing the function space (see sp_vector)
%   mshS:    object defining the domain partition and the quadrature rule (see msh_cartesian)
%   mshT:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   coeff: function handle to compute the source function
%
% OUTPUT:
%
%   rhs: assembled right-hand side
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

function rhs = op_f_v_st_tp(spaceS,spaceT,mshS,mshT,coeff)

%  for icomp = 1:spaceS.ncomp_param
%    for idim = 1:mshS.ndim
%      size1 = size(spaceS.scalar_spaces{icomp}.sp_univ(idim).connectivity);
%      if (size1(2) ~= mshS.nel_dir(idim))
%        error ('The discrete space is not associated to the mesh')
%      end
%    end
%  end

  rhs = zeros (spaceS.ndof,spaceT.ndof);

  for iel = 1:mshS.nel_dir(1)
    msh_col = msh_evaluate_col (mshS, iel);
    sp_col  = sp_evaluate_col (spaceS, msh_col, 'value', true, 'gradient', false);
    
    msh_colT = msh_precompute(mshT);
    sp_colT = sp_precompute(spaceT,msh_colT);

    for idim = 1:mshS.rdim
      x{idim} = repmat(reshape(msh_col.geo_map(idim,:,:),msh_col.nqn,msh_col.nel),[1,1,msh_colT.nqn,msh_colT.nel]);
    end
    x{idim+1} = repmat(reshape(msh_colT.map(msh_colT.qn),1,1,msh_colT.nqn,msh_colT.nel),[msh_col.nqn,msh_col.nel,1,1]);

    rhs = rhs + op_f_v_st(sp_col,sp_colT,msh_col,msh_colT,coeff(x{:}));
  end

  rhs = rhs(:);

end