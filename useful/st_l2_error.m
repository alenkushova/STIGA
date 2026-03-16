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

function errl2 = st_l2_error (space, spaceT, msh, mshT, u, coeff)
 u=reshape(u,space.ndof,spaceT.ndof);
 eut = []; valu = [];
 for sdof = 1:space.ndof
  [eu, ~] = sp_eval_msh (u(sdof,:), spaceT, mshT);
  eut = cat(1,eut,eu(:)');
 end
 for tpnt = 1: size(eut,2)
  [eus, ~] = sp_eval_msh(eut(:,tpnt), space, msh);
  valu = cat(3,valu,eus);
 end
 valu = reshape(valu, space.ncomp, msh.nqn, msh.nel, mshT.nqn, mshT.nel);
 coeff = reshape(coeff, space.ncomp, msh.nqn, msh.nel, mshT.nqn, mshT.nel);
 ws = msh.quad_weights .* msh.jacdet;
 errl2_space = squeeze( sum ( ws.*reshape(sum((valu-coeff).^2,1),...
                                  msh.nqn, msh.nel, mshT.nqn, mshT.nel),...
                              [1 2]));
 wt = mshT.quad_weights .* mshT.jacdet;
 errl2 = sqrt(sum(wt.*errl2_space,[1 2]));
end