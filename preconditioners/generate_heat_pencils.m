% GENERATE_HEAT_PENCILS: <description>
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

function out = generate_heat_pencils(dim, Pmsh, Pspace)

out = {};

% matrix assembly:
Msx = op_u_v_tp (Pspace{1,2}, Pspace{2,2}, Pmsh{1,2});           % mass in space x
Mt  = op_u_v_tp (Pspace{1,end}, Pspace{2,end}, Pmsh{1,end});     % mass in time
Asx = op_gradu_gradv_tp (Pspace{1,2}, Pspace{2,2}, Pmsh{1,2});   % stif in space x
At  = op_gradu_v_tp (Pspace{1,end}, Pspace{2,end}, Pmsh{1,end}); % first derivative (Wt matrix nor sym. nor antisym.) 
At  = At'; % correct definition;

% set internal dofs and enforce symmetric matrices 
Mt  = (Mt(2:end,2:end) + Mt(2:end,2:end)')/2;
Msx = (Msx(2:end-1,2:end-1) + Msx(2:end-1,2:end-1)')/2;
Asx = (Asx(2:end-1,2:end-1) + Asx(2:end-1,2:end-1)')/2;
At = At(2:end,2:end); % not sym.

out{1,end+1} ='At';  out{2,end} = At; 
out{1,end+1} ='Mt';  out{2,end} = Mt; 
out{1,end+1} ='Asx'; out{2,end} = Asx; 
out{1,end+1} ='Msx'; out{2,end} = Msx; 

% if the dimension in space is at least two
if dim > 1 
Msy = op_u_v_tp (Pspace{1,3}, Pspace{2,3}, Pmsh{1,3});           % mass in space y
Asy = op_gradu_gradv_tp (Pspace{1,3}, Pspace{2,3}, Pmsh{1,3});   % stif in space y
% set internal dofs and enforce symmetric matrices 
Msy = (Msy(2:end-1,2:end-1) + Msy(2:end-1,2:end-1)')/2;
Asy = (Asy(2:end-1,2:end-1) + Asy(2:end-1,2:end-1)')/2;
out{1,end+1} ='Asy'; out{2,end} = Asy; 
out{1,end+1} ='Msy'; out{2,end} = Msy; 
end

if dim == 3 % if it is 3d in spacce
Msz = op_u_v_tp (Pspace{1,4}, Pspace{2,4}, Pmsh{1,4});           % mass in space z
Asz = op_gradu_gradv_tp (Pspace{1,4}, Pspace{2,4}, Pmsh{1,4});   % stif in space z
% set internal dofs and enforce symmetric matrices 
Msz = (Msz(2:end-1,2:end-1) + Msz(2:end-1,2:end-1)')/2;
Asz = (Asz(2:end-1,2:end-1) + Asz(2:end-1,2:end-1)')/2;
out{1,end+1} ='Asz'; out{2,end} = Asz; 
out{1,end+1} ='Msz'; out{2,end} = Msz; 
end

end
