% SCHRODINGER_ST_PROBLEM: <description>
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

function [Aint, rhs, u, int_dofs, Afun, matrices] = schrodinger_st_problem(msh, space, problem_data, method_data)
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

% Assembly the matrices.
fprintf('Assembling A... \n\n')
Mt = op_u_v_tp               ( space.tsp_trial, space.tsp_test, msh.tmsh); % mass in time
Lt = op_gradu_gradv_tp       ( space.tsp_trial, space.tsp_test, msh.tmsh); % stif in time
Wt = op_gradu_v_tp           ( space.tsp_trial, space.tsp_test, msh.tmsh); % Wt rows = trial, cols = test, so we have to transpose it. 
Wt = Wt'; 

% when we have a geometry
Ms = op_u_v_tp               ( space.xsp_trial, space.xsp_test, msh.xmsh); % mass in space
Bs = op_laplaceu_laplacev_tp ( space.xsp_trial, space.xsp_test, msh.xmsh); % Blap in space
Ls = op_laplaceu_v_tp        ( space.xsp_trial, space.xsp_test, msh.xmsh); % Lap(u) vs w 
% Ls = op_gradu_gradv_tp     ( space.xsp_trial, space.xsp_test, msh.xmsh); % stiff in space 
% Ls = -Ls ;
% symmetrize 
Mt = (Mt+Mt')/2; Ms = (Ms+Ms')/2;
Lt = (Lt+Lt')/2; Bs = (Bs+Bs')/2;

matrices.Mt = Mt(2:end,2:end);
matrices.Lt = Lt(2:end,2:end);
matrices.Wt = Wt(2:end,2:end);
matrices.Ms = Ms(2:end-1,2:end-1);
matrices.Ls =-1*(Ls(2:end-1,2:end-1)' + Ls(2:end-1,2:end-1) )/2;
matrices.Bs = Bs(2:end-1,2:end-1);
% %___ Per calcolare condizionamento di B*LsMsLs
% mat = matrices;
% DDs  = trial_degree(1); DDt  = trial_degree(end); 
% NNs = nsub(1); NNt = nsub(end);
% filename = 'test_condizionamento_LML_B';
% filename2 = ['eigs_nel_' num2str(NNs) '_p_' num2str(DDs)];
% load(filename)
% eigs_A = eigs(mat.Bs,(mat.Ls*(mat.Ms\mat.Ls)),size(mat.Bs,1));
% cond_LML_B(DDs-1,log2(NNs)-2) = max(eigs_A)/min(eigs_A);
% save(filename,"cond_LML_B")
% save(filename2,"eigs_A","mat")
% %___

fprintf('Building rhs... \n\n')
F = op_f_sv_st_tp (space.xsp_test, space.tsp_test, msh.xmsh, msh.tmsh, f); % here is the projection of f in S(V)

if isequal(solver,'MB') 
  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
  Aint = A(int_dofs,int_dofs); 
  Afun = @(x) Aint*x;

  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt - A(int_dofs, 1:space.xsp_test.ndof)*u_iniz ; 
  % La riga sopra va completata! 
elseif isequal(preconditioner ,'ichol') || isequal(preconditioner ,'ilu') || isequal(solver,'CG')
  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
  Aint = A(int_dofs,int_dofs); 
  Afun = @(x) reshape(Ms(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Lt(2:end,2:end)') + ...
                      Bs(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Mt(2:end,2:end)') + ...
                      Ls(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(1i* Wt(2:end,2:end))  - ...
                     (Ls(x_int_dofs,x_int_dofs)')*reshape(x,intnx,intnt) *(1i*(Wt(2:end,2:end)')),  ...
                      intnx*intnt, 1);
  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt;
else
  Aint = 0;
%  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
  Afun = @(x) reshape(Ms(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Lt(2:end,2:end)') + ...
                      Bs(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Mt(2:end,2:end)') + ...
                      Ls(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(1i* Wt(2:end,2:end))  - ...
                     (Ls(x_int_dofs,x_int_dofs)')*reshape(x,intnx,intnt) *(1i*(Wt(2:end,2:end)')),  ...
                      intnx*intnt, 1);

  Adrch2int = @(x) reshape(Ms(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(    Lt(2:end,2:end)') + ...
                           Bs(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(    Mt(2:end,2:end)') + ...
                           Ls(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(1i* Wt(2:end,2:end) ) - ...
                          (Ls(x_drchlt_dofs,x_int_dofs)')*reshape(x,drchnx,intnt) *(1i*(Wt(2:end,2:end)')),  ...
                           intnx*intnt, 1);
  Ainit2int = @(x) reshape(Ms(x_int_dofs,:)  *reshape(x,[],1) *(    Lt(2:end,1)') + ...
                           Bs(x_int_dofs,:)  *reshape(x,[],1) *(    Mt(2:end,1)') + ...
                           Ls(x_int_dofs,:)  *reshape(x,[],1) *(1i* Wt(1,2:end) ) - ...
                          (Ls(:,x_int_dofs)')*reshape(x,[],1) *(1i*(Wt(2:end,1)')),  ...
                           intnx*intnt, 1);

  % fprintf('Lifting... \n\n')
  % rhs  = F(int_dofs) - Adrch2int(u_drchlt(space.xsp_trial.ndof+1:end));
  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - Adrch2int(u_drchlt) - Ainit2int(u_iniz); 
end

% we should also give in output the univariate matrices.. 
end
