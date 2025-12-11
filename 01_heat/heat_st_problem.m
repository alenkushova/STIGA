% HEAT_ST_PROBLEM: <description>
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

function [Afun, rhs, u, int_dofs, Aint] = heat_st_problem(msh, space, problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Apply Dirichlet boundary conditions (strongly)  
u = zeros (space.tsp_trial.ndof*space.xsp_trial.ndof, 1);

% Only space-boundary degrees of freedom
[~, x_drchlt_dofs] = sp_drchlt_l2_proj (space.xsp_trial, msh.xmsh, ifun, drchlt_sides);

% Space-time boundary degrees of freedom ----------------------------------
[u_drchlt, drchlt_dofs, u_iniz] = heat_st_boundary_data(...
  space.xsp_trial, space.tsp_trial, msh.xmsh, msh.tmsh, dfun, drchlt_sides, 'no');
% Here 'yes' is an optional variable (default 'no'). It referst to the 
% question: is the inital condition applied weakly? If yes, then the first 
% face of the ST-cylinder is to be considered part Dirichlet boundary and 
% part Internal dof. If no, the whole first face is of Dirichlet boundary.
%--------------------------------------------------------------------------
% Applying boundary conditions in strong form
u(drchlt_dofs) = u_drchlt;
% Applying initial conditions in storng form
u(1:space.xsp_trial.ndof) = u_iniz;

%  Space-time internal degrees of freedom
int_dofs = setdiff (1:space.tsp_trial.ndof*space.xsp_trial.ndof, cat(1,(1:space.xsp_trial.ndof)',drchlt_dofs)); 
% Only space-internal degrees of freedom 
x_int_dofs = setdiff (1:space.xsp_trial.ndof, x_drchlt_dofs); 

% dimension of the vectors in space or time only 
intnx = numel(x_int_dofs); 
drchnx = numel(x_drchlt_dofs); 
intnt = space.tsp_trial.ndof-1;

switch method
  case 'Galerkin'
    % N.B. time direction is the last direction in the physic domain, so it 
    % is the first factor in the kronecker product definition and we have: 
    %     A := Wt x Ms + Mt x Ls 
    % the representative matrix for this method.

    % Assembling the matrices:
    fprintf('Building A... \n\n')
    Mt = op_u_v_tp ( space.tsp_trial, space.tsp_test, msh.tmsh);           % mass in time
    Ms = op_u_v_tp ( space.xsp_trial, space.xsp_test, msh.xmsh);           % mass in space
    Ls = op_gradu_gradv_tp ( space.xsp_trial, space.xsp_test, msh.xmsh);   % stif in space
    Wt = op_gradu_v_tp (space.tsp_trial, space.tsp_test, msh.tmsh);        % first derivative (Wt matrix) 
    % Since it is (trial rows) x (test columns), we transpose to get the
    % correct sizes and definitions
    Wt = Wt'; 
    % Enforce symmetry of symmetric matrices
     Mt = (Mt+Mt')/2; Ms = (Ms+Ms')/2; Ls = (Ls+Ls')/2;

    Afun = @(x) reshape(Ms(x_int_dofs,x_int_dofs) *reshape(x,intnx,intnt)*(Wt(2:end,2:end)') + ...
                        Ls(x_int_dofs,x_int_dofs) *reshape(x,intnx,intnt)*(Mt(2:end,2:end)'),  ...
                        intnx*intnt,1);

    Adrch2int = @(x) reshape(Ms(x_int_dofs,x_drchlt_dofs)*reshape(x,drchnx,intnt)*(Wt(2:end,2:end)') + ...
                             Ls(x_int_dofs,x_drchlt_dofs)*reshape(x,drchnx,intnt)*(Mt(2:end,2:end)'), ...
                             intnx*intnt, 1);

    Ainit2int = @(x) reshape(Ms(x_int_dofs,:)*reshape(x,[],1)*(Wt(2:end,1)') + ...
                             Ls(x_int_dofs,:)*reshape(x,[],1)*(Mt(2:end,1)'), ...
                             intnx*intnt, 1);

    % Assembling the rhs: 
    fprintf('Building rhs... \n\n')  
    if is_separable_u
      % here is the projection of f when it is separable
      fht1 = op_f_v_tp (space.tsp_test, msh.tmsh, f1);
      ghx1 = op_f_v_tp (space.xsp_test, msh.xmsh, g1);
      fht2 = op_f_v_tp (space.tsp_test, msh.tmsh, f2);
      ghx2 = op_f_v_tp (space.xsp_test, msh.xmsh, g2);
      F = kron(fht1, ghx1) + kron(fht2, ghx2);
    else
      % here is the projection of f when it is NOT separable, 
      F = op_f_v_st_tp (space.xsp_test, space.tsp_test, msh.xmsh, msh.tmsh, f); 
    end 
    fprintf('Lifting... \n\n')
    rhs  = F(int_dofs) - Adrch2int(u_drchlt) - Ainit2int(u_iniz);

    if nargout == 5
      ms= Ms(x_int_dofs,x_int_dofs);
      mt= Mt(2:end,2:end);
      ls= Ls(x_int_dofs,x_int_dofs);
      wt= Wt(2:end,2:end);
      % Enforce symmetry of symmetric matrices
      mt = (mt+mt')/2; ms = (ms+ms')/2; ls = (ls+ls')/2;
      Aint =  kron(wt,ms) + kron(mt,ls) ;
    end
  case 'Least-Squares'
    % N.B. time direction is the last direction in the physic domain, so it 
    % is the first factor in the kronecker product definition and we have: 
    %     A := Lt x Ms + Mt x Bs + Wt x Ls' + Wt' x Ls  
    % the representative matrix for this method.

    % Assembly the matrices.
    fprintf('Assembling A... \n\n')
    Mt = op_u_v_tp ( space.tsp_trial, space.tsp_test, msh.tmsh);           % mass in time
    Ms = op_u_v_tp ( space.xsp_trial, space.xsp_test, msh.xmsh);           % mass in space
    Lt = op_gradu_gradv_tp( space.tsp_trial, space.tsp_test, msh.tmsh);    % stif in time
    Ls = op_laplaceu_v_tp ( space.xsp_trial, space.xsp_test, msh.xmsh);    % laplace in space (not symmetric with b.c.)
    Bs = op_laplaceu_laplacev_tp ( space.xsp_trial, space.xsp_test, msh.xmsh); % Bilap in space
    Wt = op_gradu_v_tp (space.tsp_trial, space.tsp_test, msh.tmsh);        % first derivative (Wt matrix nor sym. nor antisym.) 
    % Since it is (trial rows) x (test columns), we transpose to get the
    % correct sizes and definitions
    Wt = Wt'; 
    % Enforce symmetry of symmetric matrices
    Mt = (Mt+Mt')/2; Ms = (Ms+Ms')/2;
    Lt = (Lt+Lt')/2; Bs = (Bs+Bs')/2;

    Afun = @(x) reshape(Ms(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(Lt(2:end,2:end)') + ...
                        Bs(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(Mt(2:end,2:end)') - ...
                        Ls(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(Wt(2:end,2:end) ) - ...
                       (Ls(x_int_dofs,x_int_dofs)')*reshape(x,intnx,intnt) *(Wt(2:end,2:end)'),  ...
                        intnx*intnt, 1);

    Adrch2int = @(x) reshape(Ms(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(Lt(2:end,2:end)') + ...
                             Bs(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(Mt(2:end,2:end)') - ...
                             Ls(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(Wt(2:end,2:end) ) - ...
                            (Ls(x_drchlt_dofs,x_int_dofs)')*reshape(x,drchnx,intnt) *(Wt(2:end,2:end)'),  ...
                             intnx*intnt, 1);

    % Assembling the rhs:  
    fprintf('Building rhs... \n\n')
    if is_separable_u
      % here is the projection of f when it is separable
      fht1m = op_f_v_tp (space.tsp_test, msh.tmsh, f1);
      ghx1m = op_f_v_tp (space.xsp_test, msh.xmsh, g1);
      fht2m = op_f_v_tp (space.tsp_test, msh.tmsh, f2);
      ghx2m = op_f_v_tp (space.xsp_test, msh.xmsh, g2);

      % N.B. required OP_F_GRADV_TP and OP_F_LAPV_TP
      fht1dt = op_f_gradv_tp (space.tsp_test, msh.tmsh, f1);
      ghx1lap = op_f_lapv_tp  (space.xsp_test, msh.xmsh, g1);
      fht2dt = op_f_gradv_tp (space.tsp_test, msh.tmsh, f2);
      ghx2lap = op_f_lapv_tp  (space.xsp_test, msh.xmsh, g2);

      F = kron(fht1dt,ghx1m) + kron(fht1m,ghx1lap) +...
          kron(fht2dt,ghx2m) + kron(fht2m,ghx2lap) ;
    else
      % here is the projection of f when it is NOT separable, 
      % it requires xt_structures: geo/ msh/ space. (TO BE IMPLEMENTED YET)
      % N.B. requires also op_f_heatv_tp operator. 
      F = op_f_heatv_tp (space.xtsp_test, msh.xtmsh, f);  
    end

    fprintf('Lifting... \n\n')
    rhs  = F(int_dofs) - Adrch2int(u_drchlt(space.xsp_trial.ndof+1:end));

    if nargout == 5
      ms= Ms(x_int_dofs,x_int_dofs);
      mt= Mt(2:end,2:end);
      ls= Ls(x_int_dofs,x_int_dofs);
      lt= Lt(2:end,2:end);
      bs= Bs(x_int_dofs,x_int_dofs);      
      wt= Wt(2:end,2:end);
      % Enforce symmetry of symmetric matrices
      mt = (mt+mt')/2; ms = (ms+ms')/2; 
      lt = (lt+lt')/2; ls = (ls+ls')/2; 
      bs = (bs+bs')/2;
      Aint =  kron(lt,ms) + kron(mt,bs) -kron(wt,ls') -kron(wt',ls);
    end
end
end

