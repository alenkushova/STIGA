% PROBLEM_NAME: <description>
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
clear; close all; clc;

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
T = 1; problem_data.T = T ;    

% problem_data.xt_geo_name = 'geo_rectangle.txt'; % for T = 2
problem_data.xt_geo_name = 'geo_square.txt'; % for T = 1
problem_data.x_geo_name = nrbline ([0 0], [1 0]); % univariate geo in space
problem_data.t_geo_name = nrbline ([0 0], [T 0]); % univariate geo in time

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic


% Parameters:
Mm  = 1.5;
T_0 = 1.5;
beta= 2.5;
R = @(t) sqrt(T_0^2 - 1i*beta*t);

% Exact solution:
problem_data.uex     = @(x, t) (Mm*T_0./R(t)).*exp(-(x.^2)./R(t));
problem_data.graduex = @(x, t) (cat (1, ...
                reshape (-2*x./R(t).*problem_data.uex(x,t) , [1, size(x)]), ...
                reshape (1i*beta*(R(t)-x.^2)./(2*R(t).^3).*problem_data.uex(x,t) , [1, size(x)])));

% Source term
problem_data.f = @(x, t) ...
            (beta*x.^2 - (beta+8*x.^2).*R(t) + 4*R(t).^2)./(2*R(t).^3).*problem_data.uex(x,t);

% Dirichlet data 
problem_data.dfun = @(x, t, iside) problem_data.uex(x, t);
problem_data.ifun = @(x, iside) problem_data.uex(x, 0);
problem_data.space_dimension = '1D'; 

% coefficients
problem_data.gmm = 1i;
problem_data.eta = 1;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
n = 4; % elements per univariate direction
d = 2;  % polynomial degrees 
method_data.trial_degree     = [d d-1];                                    % Degree of the trial splines (last is time dir)
method_data.trial_regularity = method_data.trial_degree-1;                 % Regularity of trial the splines
method_data.test_degree      = [d d-1];                                    % Degree of the test splines (last is time dir)
method_data.test_regularity  = method_data.test_degree-1;                  % Regularity of the test splines
method_data.nsub       = [n n*T];                                          % Number of subdivisions
method_data.nquad      = method_data.trial_degree+1;                       % Points for the Gaussian quadrature rule

% LIST OF SOLVERS:
% 'MB'  = Matlab Backslash 
% 'CG'  = Conjugate Gradients
% 'PCG' = Preconditioned Conjugate Gradients
method_data.solver     = 'PCG'; 

% LIST OF PRECONDITIONERS:
% 'LUFD'  = block LU with FD in space
% 'ilu'   = incomplete LU factorization
% 'ichol' = incomplete Cholesky factorization
if isequal(method_data.solver,'MB') || isequal(method_data.solver,'CG')
  method_data.preconditioner = 'none'; % no preconditioner
else    
  method_data.preconditioner = 'LUFD'; % or choose another preconditioner 
end

% 3) CALL TO THE SOLVER 
[geometry, msh, space, u, report] = schroedinger_st_solve (problem_data, method_data);
report

%% 4) POST PROCESSING 
vtk_pts = {linspace(0, 1, 100), linspace(0, T, 100)};
[eu, F] = sp_eval (u, space.xtsp_trial, geometry.xtgeo, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure ('Units', 'pixels', 'Position', [150 200 1000 350])
subplot (1,2,1)
h1 = pcolor (X, Y, real(eu));
colorbar
colormap jet
h1.EdgeColor = 'none';
h1.FaceColor = 'interp';
title ('Numerical solution: $$\Re(u_h)$$','Interpreter','latex'), axis tight
xlabel('Space','Interpreter','latex')
ylabel('Time','Interpreter','latex')
subplot (1,2,2)
h2 = pcolor (X, Y, real(problem_data.uex (X,Y)));
colorbar
colormap jet
h2.EdgeColor = 'none';
h2.FaceColor = 'interp';
title ('Exact solution: $$\Re(u)$$','Interpreter','latex'), axis tight
xlabel('Space','Interpreter','latex')
ylabel('Time','Interpreter','latex')

%% 5) COMPUTE THE ERRORS
Uex = @(x, t) (problem_data.uex(x, t));
rhs = @(x, t) (problem_data.f(x, t));
[error_Graph, error_l2, error_sGraph] = schroedinger_graph_error(space.xtsp_trial, msh.xtmsh, u, Uex, rhs); % ABSOLUTE errors of the approximation
[U_Graph, U_l2, U_sGraph] = schroedinger_graph_error(space.xtsp_trial, msh.xtmsh, 0*u, Uex, rhs); % Norm of solution u
REL_ERR_Graph = error_Graph/U_Graph % RELATIVE error in GRAPH-norm
REL_ERR_l2 = error_l2/U_l2          % RELATIVE error in L2-norm

%% 6) SAVE NUMERICAL SOLUTION
filename = ['smooth_ps_pt/ST_SCHRODINGER_SMOOTH_1D_' method_data.preconditioner '_' method_data.solver '_ps_' num2str(d) '_pt_' num2str(d-1) '_Nt_' num2str(n*T) '_final_time_' num2str(T) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename); 
