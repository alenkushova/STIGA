% ARROW_HEAT_APPLICATION: <description>
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


function u_arrow_ = arrow_heat_application(D_arrow_Ds,g,s,Ux,Uy,Uz,Ut,rhs)
nsx = size(Ux,1); nsy = size(Uy,1); nsz = size(Uz,1);
nt  = size(Ut,1);
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nsy,nsz,nt), {Ux' Uy' Uz' Ut'}, 1:4);  
tilde_rhs = reshape (permute(tilde_rhs,[4 1 2 3]), nt, nsx*nsy*nsz);
% STEP 2 
% tilde_u = ARROW(D_arrow_Ds,g,s,tilde_rhs);
tilde_v = cat(1, tilde_rhs(1:end-1,:),...
           (tilde_rhs(end,:) + pagemtimes(g', tilde_rhs(1:end-1,:)./D_arrow_Ds(1:end-1,:)))./ s);
tilde_u = cat(1,(tilde_v(1:end-1,:)- tilde_v(end,:).*g)./D_arrow_Ds(1:end-1,:), tilde_v(end,:)); 
% STEP 3
new_tilde_u = reshape(permute(tilde_u,[2 3 4 1]), nsx, nsy, nsz, nt);
u_arrow_ = tmprod(new_tilde_u, {Ux Uy Uz Ut}, 1:4);  
u_arrow_ = u_arrow_(:);
end
