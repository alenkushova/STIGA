function [lam, cond_num] = conditioning(p, nel)
% Set the discretization
data    = method_data (p, nel, 1);
dataset = set_discretization (data);
msh    = dataset.xmsh;
space  = dataset.xspace_trial;

% univariate matrices
Ms = op_u_v_tp               ( space, space, msh); % mass in space
Bs = op_laplaceu_laplacev_tp ( space, space, msh); % Blap in space
Ls = op_laplaceu_v_tp        ( space, space, msh); % Lap(u) vs w 

% symmetrize and take internal data
Ms = (Ms+Ms')/2; Ms = Ms(2:end-1,2:end-1);
Bs = (Bs+Bs')/2; Bs = Bs(2:end-1,2:end-1);
Ls =-1*(Ls(2:end-1,2:end-1)' + Ls(2:end-1,2:end-1) )/2; 

% univariate eigenvalues and conditioning number of [(L(M)-1L)-1]*B
lam.dim1 = eigs(Bs,(Ls*(Ms\Ls)),size(Bs,1));
cond_num = max(lam.dim1)/min(lam.dim1);

% for the case of 32 number of elements we compute also the 2D eigenvalues.
if nel == 32
  R = 2*kron(Ls,Ls);
  B = kron(Bs,Ms) + kron(Ms,Bs);
  LML = kron(Ls*(Ms\Ls),Ms)+kron(Ms,Ls*(Ms\Ls));
  lam.dim2 = eigs(B+R, LML+R, size(B,1));
else
  lam.dim2 = "empty computation";
end
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