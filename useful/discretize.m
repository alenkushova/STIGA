% DISCRETIZE: generate a discretization in space-time keeping space and
%             time as separated components. This function does not build
%             any space-time cylinder or whole space-time msh/space of
%             splines. The reason is that geopdes does not handle 4d
%             geometry files. Plus, the best approach whoul be to keep
%             TENSOR-PRODUCT structure between SPACE and TIME so that the
%             Galerkin matrices storage cost is comparable to Galerkin +
%             time-stepping schemes. (could be better if we have cartesian
%             domains. For that case we will design an improoved function.
%
%   Copyright (C) Alen Kushova, Freiburg, 4/12/2025
%

function [geo, msh, space] = discretize(problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% load geometies
% xtgeo= geo_load(xt_geo_name);                                            % if neeeded 
xgeo = geo_load(x_geo_name);
tgeo = geo_load(t_geo_name);
% define a stucture "geo.__" with the above fields
geo = struct('xgeo',xgeo,'tgeo',tgeo); 
%geo = struct('xtgeo',xtgeo,'xgeo',xgeo,'tgeo',tgeo);                      % if neeeded 

% define mesh structures :
rdim = numel(trial_degree); % how many univariate dimensions we have.

% space time knots for trial functions                                     % if neeeded 
% [knots, zeta]    = kntrefine(xtgeo.nurbs.knots, nsub-1, trial_degree, trial_regularity);
% knots   = kntunclamp(knots, trial_degree, trial_regularity, []);

% space knots for trial functions in space
[x_knots, x_zeta]= kntrefine(xgeo.nurbs.knots, nsub(1:rdim-1)-1, trial_degree(1:rdim-1), trial_regularity(1:rdim-1));
 x_knots         = kntunclamp(x_knots, trial_degree(1:rdim-1), trial_regularity(1:rdim-1), []);

%time knots for trial functions in time
[t_knots, t_zeta]= kntrefine(tgeo.nurbs.knots, nsub(rdim)-1, trial_degree(rdim), trial_regularity(rdim));
 t_knots         = kntunclamp(t_knots, trial_degree(rdim), trial_regularity(rdim), []);

% define quadrature rules 
rule         = msh_gauss_nodes (nquad);

% [qn, qw]  = msh_set_quad_nodes (zeta, rule);                             % if neeeded 
%   xtmsh   = msh_cartesian (zeta, qn, qw, xtgeo);                         

[xqn, xqw] = msh_set_quad_nodes (x_zeta, rule(1:end-1));
   xmsh    = msh_cartesian (x_zeta, xqn, xqw, xgeo);

[tqn, tqw] = msh_set_quad_nodes (t_zeta, rule(end)); 
   tmsh    = msh_cartesian (t_zeta, tqn, tqw, tgeo);

% define a stucture "msh.__" with the above fields
msh = struct('xmsh',xmsh,'tmsh',tmsh);
% msh = struct('xtmsh',xtmsh,'xmsh',xmsh,'tmsh',tmsh);                     % if neeeded 

% define space structures for trial functions
% xtsp_trial = sp_bspline (knots,   trial_degree,          xtmsh);         % if neeeded 
xsp_trial  = sp_bspline (x_knots, trial_degree(1:rdim-1), xmsh);
tsp_trial  = sp_bspline (t_knots, trial_degree(rdim),     tmsh);

% space time knots for test functions                                      % if neeeded 
% [knots, ~]    = kntrefine(xtgeo.nurbs.knots, nsub-1, test_degree, test_regularity); 
% knots   = kntunclamp(knots, test_degree, test_regularity, []);          

% space knots for test functions in space
[x_knots, ~]= kntrefine(xgeo.nurbs.knots, nsub(1:rdim-1)-1, test_degree(1:rdim-1), test_regularity(1:rdim-1));
 x_knots    = kntunclamp(x_knots, test_degree(1:rdim-1), test_regularity(1:rdim-1), []);

%time knots for test functions in time
[t_knots, ~]= kntrefine(tgeo.nurbs.knots, nsub(rdim)-1, test_degree(rdim), test_regularity(rdim));
 t_knots    = kntunclamp(t_knots, test_degree(rdim), test_regularity(rdim), []);

% define space structures for test functions
% xtsp_test = sp_bspline (knots,   test_degree,          xtmsh);           % if neeeded 
xsp_test  = sp_bspline (x_knots, test_degree(1:rdim-1), xmsh);
tsp_test  = sp_bspline (t_knots, test_degree(rdim),     tmsh);

% define a stucture "space.__" with the above fields
space = struct('xsp_trial',xsp_trial,'tsp_trial',tsp_trial,'xsp_test', xsp_test, 'tsp_test', tsp_test);
end
