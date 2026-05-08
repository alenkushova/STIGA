% LU_STOKES_SETUP: <description>
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

function [precv, precp] = lu_stokes_setup(inputcell)
 % store the elements into variables with cell's lables
 names = cell2struct(inputcell(2,:), inputcell(1,:), 2); 
 data_names = fieldnames (names);
 for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= names.(data_names{iopt});']);
 end
 % Defoult 
 Dx = 0; Dy = 0; Dz = 0; 
 Dxp = 1; Dyp = 1; Dzp = 1; 
 dx2 = decomposition(1);  dx3 = decomposition(1); 
 Ux = 1; Uy = 1; Uz = 1; 
 % Fast diagonalization in univariate space direction (SETUP) 
 [Ux, Dx] = eig(full(Asx),full(Msx),'vector');  
 Dxp = full(diag(Msxp));
 dx1 = decomposition(Msxp,'banded','CheckCondition',false);
 dim = 1;
 if exist('Asy')
  [Uy, Dy] = eig(full(Asy),full(Msy),'vector');
  Dyp = full(diag(Msyp));
  dx2 = decomposition(Msyp,'banded','CheckCondition',false);
  dim = 2;
  if exist('Asz')
  [Uz, Dz] = eig(full(Asz),full(Msz),'vector');
  Dzp = full(diag(Msxp));
  dx3 = decomposition(Mszp,'banded','CheckCondition',false);
  dim = 3;
  end
 end
 Ds = reshape(Dx+Dy'+reshape(Dz,1,1,[]),[],1);
 Dsp = reshape(Dxp.*Dyp'.*reshape(Dzp,1,1,[]),[],1);
 Scaling = sqrt(Dsp./DMspF);
 dxt = decomposition(Mt , 'banded','CheckCondition',false);
 Ut = speye(size(At,1)*dim);
 % Define block diagonal matrix (At x Id + Mt x Ds) 
 B  = kron(speye(numel(Ds)),kron(At,eye(dim))) + kron(speye(numel(Ds)).*Ds,kron(Mt,eye(dim)));
 % N.B. now time is first direction hence we use 'permute' in the
 % application of the preconditioner
 dB = decomposition(B,'banded','CheckCondition',false);
 precv = @(v) lu_stokes_application(dB, Ux, Uy, Uz, Ut, v);
 precp = @(v) lu_stokes_pressure_application(dx1,dx2,dx3,dxt,Scaling,v);
end

