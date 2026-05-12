T = 1; % final time 
ns = 5; % number of subdivisions in space
nt = 5; % number of subdivisions in time
d = 3; % polynomial degree of pressure space

problem_data.geo_time = nrbline ([0 0], [T 0]); % NURBS line in time 
problem_data.geo_space = 'QoR.txt'; % square NURBS surface as .txt
problem_data.geo_space_time = nrbextrude(geo_load(problem_data.geo_space).nurbs, [0 0 1]); % NURBS volume

% Dirichlet/Neumann sides only in space
problem_data.drchlt_sides = [1 2 3 4]; 
problem_data.nmnn_sides = []; 

% viscosity of the fluid
problem_data.viscosity =  @(x, y) ones (size (x)); 

% Initial term to be projected 
v0x = @(x, y) zeros([1, size(x)]);
v0y = @(x, y) zeros([1, size(x)]);
problem_data.ifun= @(x, y) cat(1, reshape (v0x (x, y), [1, size(x)]), ...
                     reshape (v0y (x, y), [1, size(x)]));

% Dirichlet data 
problem_data.dfun= @lid_driven_dfun;

method_data.trial_degree     = [d d d];  % degree of the trial pressure space 
method_data.trial_regularity = method_data.trial_degree-1; % regularity of the trial pressure space
method_data.test_degree     = [d d d];  % degree of the trial pressure space 
method_data.test_regularity = method_data.test_degree-1; % regularity of the trial pressure space
method_data.nsub       = [ns ns nt];  % number of subdivisions 
method_data.nquad      = method_data.trial_degree+2; % number of quadrature points (+2 cuz vel \in degree+1)


% CALL TO THE SOLVER
[geo, msh, space, vel, pres, report] = stokes_st_solve(problem_data, method_data);
report

% post-processing
nframes = 11;
plot_vel_pres(vel, pres, space, geo, 11, 'QoR_solution');


function h  = lid_driven_dfun (x, y, t, iside) 
 switch (iside)
  case {1, 3, 4}
   h = zeros ([2, size(x)]);
  case 2
   h1 = ones ([1, size(x)]);
   h2 = zeros ([1, size(x)]);
   h  = cat (1, h1, h2);
 end
end
