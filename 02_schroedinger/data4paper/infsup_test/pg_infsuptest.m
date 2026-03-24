% INF_SUP_TEST for Schroedinger space-time discrete weak formulaiton
%  p     =  degree of B-splines
%  nel   =  number of subdivisions in space direction!
%  'Trial_N'= trial space norm: string 'l2' for L^2-norm, 'G' for graph norm. 
%  'Test_N' = trial space norm: string 'l2' for L^2-norm, 'G' for graph norm.  
% 
% This tests if the infsup discrete condition is achieved.
% 
% we check inf_v sup_w (A*v,w)/ |v||w| 
% A is the matrix associated to the bilinear form.
% Mv is the matrix representing the norm on trial space.
% Mw is the matrix representing the norm on test space.
% 
% It also can compute the eigenfunction associated to the smallestabs
% eigenvalue, in order to visualize its behaviour. 
function [mu, D, geo, msh, space, eigvect] = pg_infsuptest (p, nel)
% Set the discretization
input_args = method_data (p, nel, 1);
dataset = set_discretization(input_args);
geo      = dataset.xtgeo;
msh      = dataset.xtmsh;
mshx     = dataset.xmsh;
msht     = dataset.tmsh;
space    = dataset.xtspace_trial;
xspace   = dataset.xspace_trial;
tspace   = dataset.tspace_trial;
sp_test  = dataset.xtspace_test;
xsp_test = dataset.xspace_test;
tsp_test = dataset.tspace_test;

% Assembly the matrices.
Ms = op_u_v_tp               ( xspace, xspace, mshx); % mass in space
Bs = op_laplaceu_laplacev_tp ( xspace, xspace, mshx); % Blap in space
Ls = op_laplaceu_v_tp        ( xspace, xspace, mshx); % Lap(u) vs w 
Mt = op_u_v_tp               ( tspace, tspace, msht); % mass in time
Lt = op_gradu_gradv_tp       ( tspace, tspace, msht); % stif in time
Wt = op_gradu_v_tp           ( tspace, tspace, msht); % Wt rows = trial, cols = test, so we have to transpose it. 
Wt = Wt'; 

A   = 1i*kron(Wt,Ms) - kron(Mt,Ls);
Mvs = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
L2  = kron(Mt,Ms);

Mv = Mvs + L2;
Mw = L2; % L^2-norm on test space A*(splines) --> so graph semi-norm.

% Always symmetrize the mass matrices
Mv = (Mv+Mv')/2;
Mw = (Mw+Mw')/2;

% Extract internal dofs
int_dof_t = (1:tspace.ndof)>1;
int_dof_x = ((1:xspace.ndof)>1).*((1:xspace.ndof)<xspace.ndof) ;
int_dofs = (1:space.ndof).*kron(int_dof_t,int_dof_x);
int_dofs = nonzeros(int_dofs);

Atilde  = A  (int_dofs, int_dofs);  % schroedinger op matrix int-dofs
Mwtilde = Mw (int_dofs, int_dofs);  % norm in test space
Mvtilde = Mv (int_dofs, int_dofs);  % norm in trial space

% Compute the inf-sup constants:
F = @(u) Atilde' \ (Mvtilde * (Atilde\u)); 
[w, D] = eigs ( F, size(Atilde,1), Mwtilde, 1, 'smallestabs', 'IsFunctionSymmetric', true);
mu = sqrt(real(D)); mu = sqrt(real(mu));
eigvect = zeros(space.ndof,1);
eigvect(int_dofs,1) = w;
end

function output_args = method_data (p, nel, T)
% Geometries.
output_args.xt_geo_name = 'geo_square.txt';
output_args.x_geo_name  = nrbline ([0 0], [1 0]);
output_args.t_geo_name  = nrbline ([0 0], [T 0]);

% Choice of discretization parameters. 
output_args.trial_degree     = [p p ];                                     % Degree of the splines (last is time dir)
output_args.trial_regularity = output_args.trial_degree-1;                 % Regularity of the splines ( " )
output_args.test_degree      = [p p ];                                     % Degree of the splines (last is time dir)
output_args.test_regularity  = output_args.test_degree-1;                  % Regularity of the splines ( " )
output_args.nsub  = [nel T*nel];                                           % Number of subdivisions ( " )
output_args.nquad = max(output_args.trial_degree, output_args.test_degree) +1;                             % Points for the Gaussian quadrature rule ( " )
end

function dataset = set_discretization(input_args)

% read the fields from structures
data_names = fieldnames (input_args);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= input_args.(data_names{iopt});']);
end

% define geometies
dataset.xtgeo = geo_load(xt_geo_name);
dataset.xgeo  = geo_load(x_geo_name);
dataset.tgeo  = geo_load(t_geo_name);

% define mesh structures :
rdim = numel(trial_degree);
% space time knots for trial functions
[knots, zeta]    = kntrefine(dataset.xtgeo.nurbs.knots, nsub-1, trial_degree, trial_regularity);
knots   = kntunclamp(knots, trial_degree, trial_regularity, []);

% space knots for trial functions in space
[x_knots, x_zeta]= kntrefine(dataset.xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, trial_degree(1:rdim-1), trial_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, trial_degree(1:rdim-1), trial_regularity(1:rdim-1), []);

%time knots for trial functions in time
[t_knots, t_zeta]= kntrefine(dataset.tgeo.nurbs.knots,...
                    nsub(rdim)-1, trial_degree(rdim), trial_regularity(rdim));
t_knots = kntunclamp(t_knots, trial_degree(rdim), trial_regularity(rdim), []);

% define quadrature rules 
rule         = msh_gauss_nodes (nquad);
[qn, qw]     = msh_set_quad_nodes (zeta, rule);
dataset.xtmsh= msh_cartesian (zeta, qn, qw, dataset.xtgeo);

[xqn, xqw]   = msh_set_quad_nodes (x_zeta, rule(1:end-1));
dataset.xmsh = msh_cartesian (x_zeta, xqn, xqw, dataset.xgeo);

[tqn, tqw]   = msh_set_quad_nodes (t_zeta, rule(end));
dataset.tmsh = msh_cartesian (t_zeta, tqn, tqw, dataset.tgeo);

% define space structures for trial functions
dataset.xtspace_trial = sp_bspline (knots, trial_degree, dataset.xtmsh);
dataset.xspace_trial  = sp_bspline (x_knots, trial_degree(1:rdim-1), dataset.xmsh);
dataset.tspace_trial  = sp_bspline (t_knots, trial_degree(rdim), dataset.tmsh);

% space time knots for trial functions
[knots, ~]    = kntrefine(dataset.xtgeo.nurbs.knots, nsub-1, test_degree, test_regularity);
knots   = kntunclamp(knots, test_degree, test_regularity, []);

% space knots for trial functions in space
[x_knots, ~]= kntrefine(dataset.xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, test_degree(1:rdim-1), test_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, test_degree(1:rdim-1), test_regularity(1:rdim-1), []);

%time knots for trial functions in time
[t_knots, ~]= kntrefine(dataset.tgeo.nurbs.knots,...
                    nsub(rdim)-1, test_degree(rdim), test_regularity(rdim));
t_knots = kntunclamp(t_knots, test_degree(rdim), test_regularity(rdim), []);

% define space structures for trial functions
dataset.xtspace_test = sp_bspline (knots, test_degree, dataset.xtmsh);
dataset.xspace_test  = sp_bspline (x_knots, test_degree(1:rdim-1), dataset.xmsh);
dataset.tspace_test  = sp_bspline (t_knots, test_degree(rdim), dataset.tmsh);

end