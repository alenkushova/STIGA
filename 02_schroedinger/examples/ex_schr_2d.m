close all
clear 
clc
T = 1;
%% PROBLEM DATA
problem_data.T = T ;  % Final time.
problem_data.xt_geo_name = 'geo_cube.txt'; % for T = 2
problem_data.x_geo_name  = 'geo_square.txt'; % univariate geo in space
problem_data.t_geo_name  = nrbline ([0 0], [T 0]); % univariate geo in time
problem_data.x1_geo_name = nrbline ([0 0], [1 0]); % univariate geo in time
problem_data.x2_geo_name = nrbline ([0 0], [1 0]); % univariate geo in time
problem_data.x3_geo_name = nrbline ([0 0], [1 0]); % univariate geo in time

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 4 5];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 3 4];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Parameters:
omg = 0.2;
a0  = (2/(omg^2))^(1/4);
problem_data.a0  = a0;
problem_data.omg = omg;

% Exact solution:
problem_data.uex     = @(x, y, t) (a0*exp(-1i*(x.^2 + y.^2 + t.^2)/(omg^2)));
problem_data.graduex = @(x, y, t) (cat (1, ...
                reshape (-2i*x/(omg^2).*problem_data.uex(x, y, t) , [1, size(x)]), ...
                reshape (-2i*y/(omg^2).*problem_data.uex(x, y, t) , [1, size(y)]), ...
                reshape (-2i*t/(omg^2).*problem_data.uex(x, y, t) , [1, size(t)])));
esp = @(x) (exp(-1i*(x.^2)/(omg^2)));
problem_data.uex1    = @(x) esp(x);
problem_data.uex2    = @(y) esp(y);
problem_data.uex3    = @(z) esp(z);
problem_data.uext    = @(t) esp(t);

% Source term
problem_data.f   = @(x, y, t) (4i*omg^2 + 4*x.^2 + 4*y.^2 + 2*omg^2*t)/omg^4.*problem_data.uex(x, y, t);
problem_data.fx1 = {@(x) (4i*a0/omg^2)*esp(x), @(x) (4*a0*x.^2/omg^4).*esp(x), @(x) esp(x), @(x) esp(x)};
problem_data.fx2 = {@(y) esp(y), @(y) esp(y), @(y) (4*a0*y.^2/omg^4).*esp(y), @(y) esp(y)};
problem_data.ft  = {@(t) esp(t), @(t) esp(t), @(t) esp(t), @(t) (2*a0*t/omg^2).*esp(t)};
problem_data.h = @(x, y, t, ind) problem_data.uex(x, y, t);
problem_data.hx = @(x, y, ind) zeros (size (x)); % auxiliary function. 
problem_data.gmm = 1i;
problem_data.eta = 1;
problem_data.space_dimension = '2D'; 

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
n = 8;
d = 2;
method_data.trial_degree     = [d d d-1]; % Degree of the splines (last is time dir)
method_data.trial_regularity = method_data.trial_degree-1; % Regularity of the splines
method_data.test_degree      = [d d d-1]; % Degree of the splines (last is time dir)
method_data.test_regularity  = method_data.test_degree-1; % Regularity of the splines
method_data.nsub       = [n n n]; % Number of subdivisions
method_data.nquad      = method_data.trial_degree+1; % Points for the Gaussian quadrature rule

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
  method_data.preconditioner = 'none'; % no preconditioner! DO NOT CHANGE.
else    
  method_data.preconditioner = 'LUFD'; % or choose another preconditioner 
end

%% 3) CALL TO THE SOLVER 
%[geometry, msh, space, u, report] = solve_schrodinger_st (problem_data, method_data);
[geometry, msh, space, u, report] = solve_schrodinger_st_on_cartesian_domains (problem_data, method_data);
report

%% 4) POST-PROCESSING
% output_file = 'smooth_2D/Schrodinger_st_2D';
% 
% vtk_pts = {linspace(0, 1, 200), linspace(0, 1, 200), linspace(0, 1, 20)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, space.xtsp_trial, geometry.xtgeo, vtk_pts, output_file, 'u')

% Plot of the figure:
% figure
% sp_plot_solution (real(u), space.xtsp_trial, geometry.xtgeo, [40 40 40], [20 20 20]);

%% 5) DISPLAY ERRORS 
% compute the error for the real part:
% Uex = @(x, y, t) (problem_data.uex(x, y, t));
% rhs = @(x, y, t) (problem_data.f(x, y,t));
% [error_Graph, error_l2, error_sGraph] = schroedinger_graph_error(space.xtsp_trial, msh.xtmsh, u, Uex, rhs) % ABSOLUTE errors of the approximation
% [U_Graph, U_l2, U_sGraph] = schroedinger_graph_error(space.xtsp_trial, msh.xtmsh, 0*u, Uex, rhs) % Norm of solution u
% REL_ERR_Graph = error_Graph/U_Graph; % RELATIVE error in GRAPH-norm
% REL_ERR_l2 = error_l2/U_l2;          % RELATIVE error in L2-norm
% 

%% 6) SAVE NUMERICAL SOLUTION
%filename = ['smooth_2D/ST_SCHRODINGER_SMOOTH_2D_' method_data.preconditioner '_' method_data.solver '_degree_' num2str(d) '_Nt_' num2str(n*T) '_final_time_' num2str(T) '.mat'];
filename = ['smooth_2D_ps_pt/ST_SCHRODINGER_SMOOTH_2D_' method_data.preconditioner '_' method_data.solver '_ps_' num2str(d) '_pt_' num2str(d-1) '_Nt_' num2str(n*T) '_final_time_' num2str(T) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);


