T = 1; % final time 
n = 5; % number of subdivisions in space
d = 1; % polynomial degree of pressure space

viscosity = 1; % funzione costante

% definiamo i dati esatti del problema (soluzioni e sorgenti) la soluzione
% è presa dall'articolo di Stainbach - 
syms x y t;
sym_velexx = (exp(-t)-1)*x^2*(x-1)^2*(4*y^3-6*y^2+2*y);
sym_velexy = (1-exp(-t))*y^2*(y-1)^2*(4*x^3-6*x^2+2*x);
sym_presex = (1+x-exp(-x*y*t))*t^2;
problem_data = computeProblemData2d(sym_velexx,sym_velexy,sym_presex,viscosity);

square = nrbextrude( nrbline ([0 0], [1 0]), [0,1]); % square NURBS surface
problem_data.geo_time = nrbline ([0 0], [T 0]); % NURBS line in time 
problem_data.geo_space = 'geo_square.txt'; % square NURBS surface as .txt
problem_data.geo_space_time = nrbextrude(square, [0 0 1]); % NURBS volume

% Dirichlet/Neumann sides only in space
problem_data.drchlt_sides = [1 2 3 4]; 
problem_data.nmnn_sides = []; 

% viscosity of the fluid
problem_data.viscosity =  @(x, y) ones (size (x)); 
 
% Dirichlet data 
problem_data.dfun = @(x,y,t,iside) problem_data.velex(x,y,t);

method_data.trial_degree     = [d d d];  % degree of the trial pressure space 
method_data.trial_regularity = method_data.trial_degree-1; % regularity of the trial pressure space
method_data.test_degree     = [d d d];  % degree of the trial pressure space 
method_data.test_regularity = method_data.test_degree-1; % regularity of the trial pressure space
method_data.nsub       = [n n n];  % number of subdivisions 
method_data.nquad      = method_data.trial_degree+2; % number of quadrature points (+2 cuz vel \in degree+1)

% CALL TO THE SOLVER
[geo, msh, space, vel, pres, report] = stokes_st_solve(problem_data, method_data);

% Absolute error computation
% pres_errl2 = st_l2_error_pressures_tp(space.spp,space.spt_pres,msh.xmsh,msh.tmsh,pres,problem_data.presex);
% [vel_errl2, vel_errh1s, vel_errh1t] = st_h1_error_tp(space.spv, space.spt_vel, msh.xmsh, msh.tmsh, ...
%     vel, problem_data.velex, problem_data.dt_velex, problem_data.grad_velex);
% 
% % Norm of exact solutions
% presex_l2norm = st_l2_error_pressures_tp(space.spp,space.spt_pres,msh.xmsh,msh.tmsh,0*pres,problem_data.presex);
% [velex_l2norm, velex_sh1norm, velex_th1norm] = st_h1_error_tp(space.spv, space.spt_vel, msh.xmsh, msh.tmsh, ...
%     0*vel, problem_data.velex, problem_data.dt_velex, problem_data.grad_velex);
% 
% pres_errl2_rel = pres_errl2/presex_l2norm;
% vel_errl2_rel  = vel_errl2/velex_l2norm;
% vel_errh1s_rel = vel_errh1s/velex_sh1norm;
% vel_errh1t_rel = vel_errh1t/velex_th1norm;
% 
% %filname = ['Explicit_solution_TH_results_d=' num2str(d) '_n=' num2str(n) '_T=1.mat'];
% filname = ['Explicit_solution_TH_gmres_results_d=' num2str(d) '_n=' num2str(n) '_T=1.mat'];
% save(filname)

%% post-processing
nframes = 11;
plot_vel_pres(vel, pres, space, geo, nframes, 'explicit_solution');
plot_exact_vel_pres(problem_data, geo, nframes, 'exact_explicit_solution');

