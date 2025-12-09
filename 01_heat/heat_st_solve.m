% HEAT_ST_SOLVE: <description>
% 
%   CALL: 
%
%  INPUT:
%
% OUTPUT:
%
%
% ProjectName - STIGA
% Copyright (C) 2025 Alen Kushova
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

function [geo, msh, space, u, report] = heat_st_solve(problem_data,method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Building the discrete spaces:
[geo, msh, space] = discretize(problem_data,method_data);

% Assembling the linear system and the right hand side:
if solver == 'M'
  [Afun, rhs, u, int_dofs, Aint] = heat_st_problem(msh, space, problem_data, method_data);
else 
  [Afun, rhs, u, int_dofs] = heat_st_problem(msh, space, problem_data, method_data);
end    

% Force computation in one singular thread
N = maxNumCompThreads; 
maxNumCompThreads(1);

file_preconditioner = '';

% Solving: 
switch solver
  case 'M'
    % Solution with backslash
    fprintf('Solving with MATLAB backslash... \n\n')
    tic;
    u(int_dofs)  = Aint\rhs;
    time = toc;
    report.flag = 'Solved with MATLAB backslash';
    report.time = time;
  case 'GMRES'
    tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
    % Solution with GMRES 
    fprintf('Solving with GMRES (without prec.)... \n\n')
    tic
    [u_inner, flag, rel_res, iter, res_vec] = gmres(Afun,rhs,[],tol,maxit);
    time = toc;
    u(int_dofs) = u_inner;
    % Output informations of GMRES solver iterations
    report.flag    = flag;
    report.rel_res = rel_res;
    report.iter    = iter;
    report.res_vec = res_vec;
    report.time    = time;
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
  case 'LU'
    fprintf('Assembling LU-block preconditioner.. \n\n ')
    tol   = 10^(-8); maxit = min(numel(int_dofs),200);        
    % Load the matrices of the correct size... 
    load(file_preconditioner,"At","Mt","Asx","Msx","Asy","Msy","Asz","Msz");
    tic
    P = LU_SETUP(At,Mt,Asx,Msx,Asy,Msy,Asz,Msz);
    % Solution with GMRES and Arrow preconditioner
    fprintf('Solving with GMRES and LU-block preconditioner... \n\n')
    [u_inner, flag, rel_res, iter, res_vec] = gmres(Afun,rhs,[],tol,maxit,@(x) P(x));
    time = toc;
    if max(abs(imag(u_inner)))<=10^(-15)
        u_inner = real(u_inner); % the solution is real, but may be affected by complex decomposition numbers.
    end
    u(int_dofs) = u_inner;
    % Output informations of GMRES solver iterations
    report.flag    = flag;
    report.rel_res = rel_res;
    report.iter    = iter;
    report.res_vec = res_vec;
    report.time    = time;
  case 'AR' 
    % Building the preconditioner
    fprintf('Assembling Arrow (AR) preconditioner... \n\n')
    tol   = 10^(-8); maxit = min(numel(int_dofs),200);        
    % Load the matrices of the correct size... 
    load(file_preconditioner,"At","Mt","Asx","Msx","Asy","Msy","Asz","Msz");
    tic
    P = ARROW_SETUP(At,Mt,Asx,Msx,Asy,Msy,Asz,Msz);
    % Solution with GMRES and Arrow preconditioner
    fprintf('Solving with GMRES and Arrow (AR) preconditioner... \n\n')
    [u_inner, flag, rel_res, iter, res_vec] = gmres(Afun,rhs,[],tol,maxit,@(x) P(x));
    time = toc;
    if max(abs(imag(u_inner)))<=10^(-15)
        u_inner = real(u_inner); % the solution is real, but may be affected by complex decomposition numbers.
    end
    u(int_dofs) = u_inner;
    % Output informations of GMRES solver iterations
    report.flag    = flag;
    report.rel_res = rel_res;
    report.iter    = iter;
    report.res_vec = res_vec;
    report.time    = time;
  case 'LR' 
    % Building the preconditioner
    tol   = 10^(-8); maxit = min(numel(int_dofs),200);        
    % Load the matrices of the correct size... 
    load(file_preconditioner,"At","Mt","Asx","Msx","Asy","Msy","Asz","Msz");
    tic
    P = SMW_SETUP(At,Mt,Asx,Msx,Asy,Msy,Asz,Msz);
    % Solution with GMRES and Arrow preconditioner
    fprintf('Solving with GMRES and Low-Rank preconditioner... \n\n')
    [u_inner, flag, rel_res, iter, res_vec] = gmres(Afun,rhs,[],tol,maxit,@(x) P(x));
    time = toc;
    if max(abs(imag(u_inner)))<=10^(-15)
        u_inner = real(u_inner); % the solution is real, but may be affected by complex decomposition numbers.
    end
    u(int_dofs) = u_inner;
    % Output informations of GMRES solver iterations
    report.flag    = flag;
    report.rel_res = rel_res;
    report.iter    = iter;
    report.res_vec = res_vec;
    report.time    = time;
end

% reset Matlab to full computational threads
report.NumThreads = maxNumCompThreads; % = 1 (single thread computations)
maxNumCompThreads(N);

end