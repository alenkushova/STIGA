% LU_STOKES_pressure_APPLICATION: <description>
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

function u = lu_stokes_pressure_application(dx1,dx2,dx3,dxt,Scaling,v)

nsx = dx1.MatrixSize(1); nsy = dx2.MatrixSize(1); nsz = dx3.MatrixSize(1);
nt  = dxt.MatrixSize(1);

% applichiamo l'inversa in tempo della massa.
V = reshape(v,nsx,nsy,nsz,nt);
v = permute(v,[4 1 2 3]);
v = dxt\reshape(v,nt,[]);

% scaling diagonale in spazio
v = Scaling'.*v;

% invertiamo in ordine le masse in spazio
v = reshape(v,nt,nsx,nsy,nsz);

% apply inverse of Mz
v = permute(v,[4 1 2 3]);
v = dx3\reshape(v,nsz,[]);
v = reshape(v,nsz,nt,nsx,nsy);

% apply inverse of My
v = permute(v,[4 1 2 3]);
v = dx2\reshape(v,nsy,[]);
v = reshape(v,nsy,nsz,nt,nsx);

% apply inverse of Mx
v = permute(v,[4 1 2 3]);
v = dx1\reshape(v,nsx,[]);

% scaling diagonale in spazio
u = Scaling.*reshape(v,nsx*nsy*nsz,nt);
u = u(:);

end