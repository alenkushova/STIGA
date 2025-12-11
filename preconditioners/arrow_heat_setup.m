% ARROW_HEAT_SETUP: <description>
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

function prec = arrow_heat_setup(inputcell)
 % store the elements into variables with cell's lables
 names = cell2struct(inputcell(2,:), inputcell(1,:), 2); 
 data_names = fieldnames (names);
 for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= names.(data_names{iopt});']);
 end
 % Defoult 
 Dx = 0; Dy = 0; Dz = 0; 
 Ux = 1; Uy = 1; Uz = 1; 
 % Fast diagonalization in univariate space direction (SETUP) 
 [Ux, Dx] = eig(full(Asx),full(Msx),'vector');  
 if exist('Asy')
  [Uy, Dy] = eig(full(Asy),full(Msy),'vector');
  if exist('Asz')
  [Uz, Dz] = eig(full(Asz),full(Msz),'vector');
  end
 end
 Ds = reshape(Dz+Dy'+reshape(Dx,1,1,[]),[],1);
 % Time factorization
 Wt0= At(1:end-1,1:end-1); %dati interni
 Mt0= Mt(1:end-1,1:end-1); %dati interni
        
 w  = At(1:end-1,end); %ultima riga di Wt
 m  = Mt(1:end-1,end); %ultima riga di Mt

%________ diagonalizzazione in tempo usando il complemento di Schur ______ (by Mattia Tani)
        L = chol(Mt0);
        [Ut0, Dt0] = schur(full(((L')\Wt0)/L),'complex');
        Dt0 = sqrt(-1)*imag(diag(Dt0));
        Ut0 = L\Ut0;        
%__________________________________________________________________________

 % find v
 v   = Mt0 \ (-m);
 v_1 = [v;1];   
 % find r and rho
 r_rho = v_1/sqrt(v_1' * Mt * v_1);
 r     = r_rho(1:end-1);
 rho   = r_rho(end);

 % find g and sigma
 g     = Ut0' * [Wt0 w] *r_rho; 
 sigma = r_rho' * At * r_rho;

 % here is the arrow like decomposition
 Ut = full([Ut0 r;0*r' rho]); %eigenvectors
 D  = [Dt0; sigma]; %'eigenvalues' 
 D_arrow_Ds = reshape(D + Ds', numel(D), []);
 s = D_arrow_Ds(end,:) + pagemtimes(g', g./D_arrow_Ds(1:end-1,:)); 

 prec = @(rhs) arrow_heat_application(D_arrow_Ds, g, s, Ux, Uy, Uz, Ut, rhs);
end
