% OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i).
%
%   rhs = op_f_v (spv, msh, coeff);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points
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


function rhs = op_f_v_st(spv,spt,msh,mshT,coeff)

coeff = reshape (permute(coeff,[1 2 4 3 5]),spv.ncomp,msh.nqn*mshT.nqn,msh.nel,mshT.nel);

rhs   = zeros (spv.ndof,spt.ndof);
shpv  = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
shpt  = reshape (spt.shape_functions, spt.ncomp, mshT.nqn, spt.nsh_max, mshT.nel);

for ielT = 1:mshT.nel
    indicesT = find (spt.connectivity(:,ielT));
    conn_ielT = spt.connectivity(indicesT,ielT);
for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
        jacdet_weightsS = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel), 1, msh.nqn);
        jacdet_weightsT = reshape (mshT.jacdet(:, ielT) .* mshT.quad_weights(:, ielT), 1, mshT.nqn);
        jacdet_weights = kron(jacdet_weightsT,jacdet_weightsS);
        coeff_times_jw = bsxfun (@times, jacdet_weights, coeff(:,:,iel,ielT));

        shpv_iel = repmat(reshape(shpv(:,:,:,iel),spv.ncomp,msh.nqn,spv.nsh_max,1,1),[1,1,1,mshT.nqn,spt.nsh_max]);
        shpt_iel = repmat(reshape(shpt(:,:,:,ielT),spt.ncomp,1,1,mshT.nqn,spt.nsh_max),[spv.ncomp,msh.nqn,spv.nsh_max,1,1]);
        
        shpv_iel = reshape(permute((shpv_iel.*shpt_iel),[1 2 4 3 5]),spv.ncomp,msh.nqn*mshT.nqn,spv.nsh_max,spt.nsh_max);
       
        aux_val = bsxfun (@times, coeff_times_jw, shpv_iel);
        rhs_loc = reshape(sum (sum (aux_val, 1), 2),[],numel(conn_ielT));

        indices = find (spv.connectivity(:,iel));
        conn_iel = spv.connectivity(indices,iel);
        rhs(conn_iel,conn_ielT) = rhs(conn_iel,conn_ielT) + rhs_loc(indices,indicesT); 
    else
        warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_v: singular map in element number %d', iel)
    end
end
end

end
