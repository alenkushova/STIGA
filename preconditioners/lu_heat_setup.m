% LU_HEAT_SETUP: <description>
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

function prec = lu_heat_setup(inputcell)
 % store the elements into variables with cell's lables
 names = cell2struct(inputcell(2,:), inputcell(1,:), 2); 
 data_names = fieldnames (names);
 for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= names.(data_names{iopt});']);
 end
 % Defoult 
 Dx = 0; Dy = 0; Dz = 0; 
 Ux = 1; Uy = 1; Uz = 1; Ut = speye(size(At,1));
 % Fast diagonalization in univariate space direction (SETUP) 
 [Ux, Dx] = eig(full(Asx),full(Msx),'vector');  
 if exist('Asy')
  [Uy, Dy] = eig(full(Asy),full(Msy),'vector');
  if exist('Asz')
  [Uz, Dz] = eig(full(Asz),full(Msz),'vector');
  end
 end
 Ds = reshape(Dz+Dy'+reshape(Dx,1,1,[]),[],1);
 % Define block diagonal matrix (At x Id + Mt x Ds) 
 B  = kron(speye(numel(Ds)), At) + kron(speye(numel(Ds)).*Ds, Mt);
 % N.B. now time is first direction hence we use 'permute' in the
 % application of the preconditioner
 dB = decomposition(B,'banded','CheckCondition',false);
 prec = @(v) lu_heat_application(dB, Ux, Uy, Uz, Ut, v);
end

