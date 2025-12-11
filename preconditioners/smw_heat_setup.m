% SMW_HEAT_SETUP: <description>
% 
%   CALL: 
%
%  INPUT: Ds row vector! 
%         Dt column vector!
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
 

function prec = smw_heat_setup(inputcell)
 % store the elements into variables with cell's lables
 names = cell2struct(inputcell(2,:), inputcell(1,:), 2); 
 data_names = fieldnames (names);
 for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= names.(data_names{iopt});']);
 end
 % Defoult 
 Dx = 0; Dy = 0; Dz = 0; nt = size(At,1);
 Ux = 1; Uy = 1; Uz = 1; Ut = speye(nt);
 % Fast diagonalization in univariate space direction (SETUP) 
 [Ux, Dx] = eig(full(Asx),full(Msx),'vector');  
 if exist('Asy')
  [Uy, Dy] = eig(full(Asy),full(Msy),'vector');
  if exist('Asz')
  [Uz, Dz] = eig(full(Asz),full(Msz),'vector');
  end
 end
 Ds = reshape(Dz+Dy'+reshape(Dx,1,1,[]),[],1);

% RANK = 2; % by defoult for our strategy...

 A = zeros(nt); %think of bettere definition, like keep it sparse...
 A(1:end-1,1:end-1) = At(1:end-1,1:end-1);
 R = At-A;
%________ diagonalizzazione in tempo usando il complemento di Schur ______ (by Mattia Tani)
        L = chol(Mt);
        [Ut, Dt] = schur(full(((L')\A)/L),'complex');
        Dt = sqrt(-1)*imag(diag(Dt));
        Ut = L\Ut;        
%__________________________________________________________________________

 % there are multiple choices now but we use the one from 
 %    https://arxiv.org/pdf/2403.07875   (Section 4.3),
 % where we have G and F defined as
 G = [zeros(1,nt-1) 1 ; R(end,1:nt-1) 0 ];
 F = R(:,end); F(1:nt-1,2) = zeros(nt-1,1); F(nt,2) = 1;

% Change of basis for low rank and diagonal form... 
U = Ut'*F; V = G*Ut;
D = reshape(Dt + Ds', size(Dt,1), 1, size(Ds,1)); 
% Compute the intermediate matrix in Low rank SMW formula
inv_r = pageinv(speye(2) + pagemtimes(V,(U./D)));

prec = @(rhs) smw_heat_application( D,U,V,inv_r, Ux,Uy,Uz,Ut, rhs);
end