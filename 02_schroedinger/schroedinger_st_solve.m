% SCHROEDINGER_ST_SOLVE: <description>
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
%

function [geo, msh, space, u, report] = schroedinger_st_solve (problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% define Degree, Num elem in space and in time:
DDs  = trial_degree(1); DDt  = trial_degree(end); NNs = nsub(1); NNt = nsub(end);

% Building the discrete spaces:
[geo, msh, space] = discretize(problem_data,method_data);

% Assembling the linear system and the right hand side:
if solver == 'MB'
  [Aint, rhs, u, int_dofs, Afun, mat] = schrodinger_st_problem(msh, space, problem_data, method_data);
else 
  [Aint, rhs, u, int_dofs, Afun, mat] = schrodinger_st_problem(msh, space, problem_data, method_data);
end  
%force Matlab to use one single thread:
N = maxNumCompThreads; 
maxNumCompThreads(1);

% Solving: 
switch solver
  case 'MB'
    % Solution with backslash
    fprintf('Solving with MATLAB backslash... \n\n')
    tic;
    u(int_dofs)  = Aint\rhs;
    time = toc;
    report.flag = 'Solved with MATLAB backslash';
    report.time = time;
  case 'CG'
    tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
    % Solution with CG
    fprintf('Solving with CG (without prec.)... \n\n')
    tic
    [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit);
    time = toc;
    u(int_dofs) = u_inner;
    % Output informations of CG solver iterations
    report.flag    = flag;
    report.rel_res = rel_res;
    report.iter    = iter;
    report.res_vec = res_vec;
    report.time    = time;
  case 'PCG'
    switch preconditioner 
        case 'LUFD' % solving with PCG and block LU with FD in space
          fprintf('Assembling LU-block preconditioner.. \n\n ')
          tol = 10^(-10); maxit = min(numel(int_dofs),200);   
          file_preconditioner = ['data/DATA_PREC_p_' num2str(DDs) '_Nt_' num2str(NNt) '_Ns_' num2str(NNs^3) '.mat'];
          switch space_dimension
            case '1D'
%              load(file_preconditioner,"Lt","Mt","Wt","Lsx","Msx");
              tic
%              P = LU_SETUP_x_SCHRODINGER_1D(Lt,Mt,Wt,Lsx,Msx);
              P = LU_SETUP_x_SCHRODINGER_1D(mat.Lt,mat.Mt,mat.Wt,mat.Ls,mat.Ms);
              fprintf('Solving with PCG and LU-block preconditioner... \n\n')
              [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit,@(x) P(x));
              time = toc;
            case '2D'
              load(file_preconditioner,"Lt","Mt","Wt","Lsx","Msx","Lsy","Msy");
              tic
              P = LU_SETUP_x_SCHRODINGER_2D(Lt,Mt,Wt,Lsx,Msx,Lsy,Msy);
              fprintf('Solving with PCG and LU-block preconditioner... \n\n')
              [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit,@(x) P(x));
              time = toc;
            case '3D'
              load(file_preconditioner,"Lt","Mt","Wt","Lsx","Msx","Lsy","Msy","Lsz","Msz");
              tic
              P = LU_SETUP_x_SCHRODINGER(Lt,Mt,Wt,Lsx,Msx,Lsy,Msy,Lsz,Msz);
              fprintf('Solving with PCG and LU-block preconditioner... \n\n')
              [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit,@(x) P(x));
              time = toc;
          end
          u(int_dofs) = u_inner;
          % Output informations of GMRES solver iterations
          report.flag    = flag;
          report.rel_res = rel_res;
          report.iter    = iter;
          report.res_vec = res_vec;
          report.time    = time;
      case 'ilu'
        tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
        fprintf('Solving with CG and ILU preconditioner... \n\n')
        tic
        [L,U] = ilu(Aint); 
        [u_inner, flag, rel_res, iter, res_vec] = pcg(Aint,rhs,tol,maxit,L,U);
        time = toc;
        u(int_dofs) = u_inner;
        % Output informations of CG solver iterations
        report.flag    = flag;
        report.rel_res = rel_res;
        report.iter    = iter;
        report.res_vec = res_vec;
        report.time    = time;
      case 'ichol'
        tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
        fprintf('Solving with CG and ICHOL preconditioner... \n\n')
        tic
        L = ichol(Aint); 
        [u_inner, flag, rel_res, iter, res_vec] = pcg(Aint,rhs,tol,maxit,L,L');
        time = toc;
        u(int_dofs) = u_inner;
        % Output informations of CG solver iterations
        report.flag    = flag;
        report.rel_res = rel_res;
        report.iter    = iter;
        report.res_vec = res_vec;
        report.time    = time;
    end
end

maxNumCompThreads(N);

end

function [Aint, rhs, u, int_dofs, Afun, matrices] = old_schrodinger_st_problem(msh, space, problem_data, method_data)
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
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt;
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
  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
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

  % fprintf('Lifting... \n\n')
  % rhs  = F(int_dofs) - Adrch2int(u_drchlt(space.xsp_trial.ndof+1:end));
  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt;
end

% we should also give in output the univariate matrices.. 
end

%% NOTA BENE: 
%  abbiamo il seguente sistema da invertire Au = f con 
%
%   A = Ms x Kt + \nu^2 Bs x Mt -i\nu Ls' x Wt + i\nu Ls x Wt' 
%  
%  e lo risolviamo iterativamente con PCG dove il precondizionatore é 
%
%   P = Ms x Kt + \nu^2 KsMs^(-1)Ks x Mt + i\nu Ks x Wt - i\nu Ks x Wt' 
%
%  questo perchè Ks è la matrice di stiffness ovvero -Ls (salvo errori)
%
%  CONTROLLA IL CASO 2D E IL CASO GENERICO, QUESTIONE DI SEGNI....

function prec = LU_SETUP_x_SCHRODINGER_1D(Lt,Mt,Wt,Asx,Msx)
nsx = size(Asx,1); nt = size(Lt,1);
%Fast diagonalization in univariate space direction (SETUP) 
[Us, Ds] = eig(full(Asx),full(Msx),'vector'); 
B = kron(speye(numel(Ds)), Lt) + ...
    kron(spdiags(reshape(Ds.^2,numel(Ds),1), 0, numel(Ds),numel(Ds)),Mt) - ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt') + ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt);
dB = decomposition(B,'banded','CheckCondition',false);
prec = @(v) LU_SOLVER_x_SCHRODINGER_1D(dB,Us,v,nsx,nt);
end

function prec = LU_SETUP_x_SCHRODINGER_2D(Lt,Mt,Wt,Lsx,Msx,Lsy,Msy)
nsx = size(Lsx,1); nsy = size(Lsy,1); nt = size(Lt,1);
%Fast diagonalization in univariate space direction (SETUP) 
[Ux, Dx] = eig(full(Lsx),full(Msx),'vector'); 
[Uy, Dy] = eig(full(Lsy),full(Msy),'vector');
% assembling and decomposing:
Ds = reshape(Dy + Dx',[],1);
B = kron(speye(numel(Ds)), Lt) + ...
    kron(spdiags(reshape(Ds.^2,numel(Ds),1), 0, numel(Ds),numel(Ds)),Mt) - ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt') + ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt);
dB = decomposition(B,'banded','CheckCondition',false);
prec = @(v) LU_SOLVER_x_SCHRODINGER_2D(dB,Ux,Uy,v,nsx,nsy,nt);
end

function u_lu = LU_SOLVER_x_SCHRODINGER_1D(lu_B,Ux,rhs,nsx,nt)
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nt), {Ux'}, 1);  
tilde_rhs = reshape (permute(tilde_rhs,[2 1]), nt*nsx,1);

% STEP 2
tilde_u = lu_B\tilde_rhs;
% STEP 3
new_tilde_u = reshape(permute(reshape(tilde_u,nt,nsx),[2 1]), nsx, nt);
u_lu = tmprod(new_tilde_u, {Ux }, 1);  
u_lu = u_lu(:);
end

function u_lu = LU_SOLVER_x_SCHRODINGER_2D(lu_B,Ux,Uy,rhs,nsx,nsy,nt)
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nsy,nt), {Ux' Uy'}, 1:2);  
tilde_rhs = reshape (permute(tilde_rhs,[3 1 2]), nt*nsx*nsy,1);
% STEP 2
tilde_u = lu_B\tilde_rhs;
% STEP 3
new_tilde_u = reshape(permute(reshape(tilde_u,nt,nsx,nsy),[2 3 1]), nsx, nsy, nt);
u_lu = tmprod(new_tilde_u, {Ux Uy}, 1:2);  
u_lu = u_lu(:);
end


