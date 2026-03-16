close all
clear 
clc
M = 1000;
T = 2;
%% PROBLEM DATA
problem_data.T = T ;  % Final time.
problem_data.xt_geo_name = 'geo_rectangle.txt'; % for T = 2
problem_data.x_geo_name = nrbline ([0 0], [1 0]); % univariate geo in space
problem_data.t_geo_name = nrbline ([0 0], [T 0]); % univariate geo in time

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 ];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 ];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Solution
solutions = ['solutions' num2str(M) '.mat'];
load(solutions); % or you can load previously built solutions and rhs.
problem_data.uex     = @(x, t) uex(x,t);
problem_data.graduex = @(x, t) grad_uex(x,t);

% Source term
problem_data.f = @(x, t) f(x,t);

% Dirichlet boundary conditions
problem_data.h = @(x, t, ind) problem_data.uex(x,t);
problem_data.hx = @(x, ind) zeros (size (x)); % auxiliary function. 
problem_data.space_dimension = '1D'; 

% coefficients
problem_data.gmm = 1i;
problem_data.eta = 1;

%% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
n = 8; % elements per univariate direction
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
method_data.solver     = 'MB'; %'PCG'; 

% LIST OF PRECONDITIONERS:
% 'LUFD'  = block LU with FD in space
% 'ilu'   = incomplete LU factorization
% 'ichol' = incomplete Cholesky factorization
if isequal(method_data.solver,'MB') || isequal(method_data.solver,'CG')
  method_data.preconditioner = 'none'; % no preconditioner
else    
  method_data.preconditioner = 'LUFD'; % or choose another preconditioner 
end
%% 3) CALL TO THE SOLVER 
[geometry, msh, space, u, report] = solve_schrodinger_st (problem_data, method_data);
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
filename = ['furier_ps_pt/ST_SCHRODINGER_FURIER_1D_' method_data.preconditioner '_' method_data.solver '_ps_' num2str(d) '_pt_' num2str(d-1) '_Nt_' num2str(n*T) '_final_time_' num2str(T) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);
