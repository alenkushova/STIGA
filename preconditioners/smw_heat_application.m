% SMW_HEAT_APPLICATION: <description>
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
 
function u_smw = smw_heat_application(D,U,V,R, Ux, Uy, Uz, Ut, rhs)
nsx = size(Ux,1); nsy = size(Uy,1); nsz = size(Uz,1);
nt  = size(Ut,1);
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nsy,nsz,nt), {Ux' Uy' Uz' Ut'}, 1:4);  
tilde_rhs = reshape (permute(tilde_rhs,[4 1 2 3]), nt, 1, nsx*nsy*nsz);
% STEP 2
% tilde_u   = SMW(D, U, R, V, tilde_rhs);
tilde_v1 = tilde_rhs./D;                                       % paralelizable elemntwise multiplication
tilde_v2 = pagemtimes(U,pagemtimes(R,pagemtimes(V,tilde_v1))); % paralelizable matrix vector multiplicaiton
tilde_u  = tilde_v1 - (tilde_v2./D);                           % paralelizable elemntwise multiplication

% STEP 3
new_tilde_u = reshape(permute(tilde_u,[2 3 4 1]), nsx, nsy, nsz, nt);
u_smw = tmprod(new_tilde_u, {Ux Uy Uz Ut}, 1:4);  
u_smw = u_smw(:);
end
