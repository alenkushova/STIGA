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

% Regularity parameter alpha > 0 for u = |x - 1/2|^alpha * exp(-t):
% - alpha>= 2: Classical solution (C^2). The second derivative is 
%              continuous.
% - 1<alpha<2: Strong solution (W^{2,p}). The first derivative is 
%              continuous, but the Laplacian (second derivative) blows up 
%              at the singularity x = 1/2.
% - alpha = 1: Weak solution (Lipschitz). The second derivative behaves 
%              like a Dirac Delta distribution at x = 1/2.
% - 0<alpha<1: Singular/Hölder solution. The gradient (slope) is infinite 
%              at x = 1/2.
alpha = 1.5; problem_data.alpha = alpha;


% NO NEED OF SPACE-TIME DOMAIN.
problem_data.xt_geo_is_needed = false;
problem_data.x_geo_name = nrbline ([0 0], [1 0]); % space domain (1D) 
problem_data.t_geo_name = nrbline ([0 0], [T 0]); % time domain (1D)

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Exact solution: 
problem_data.uex      = @(x, t)          (abs(x-1/2).^alpha  ).*exp(-t);   % separable solution

problem_data.grad_uex = @(x, t) reshape (alpha*(abs(x-1/2).^(alpha-1)).*sign(x-1/2).*exp(-t),...
                                         [1, size(x)]);
problem_data.dt_uex   = @(x, t) reshape (-(abs(x-1/2).^alpha  ).*exp(-t), ...
                                         [1, size(x)]);

% -------------------------------------------------------------------------
% Let  u(x,y,z,t) = g(x,y,z)*f(t) be separable, then:
%
%    Heat( u(x,y,z,t) ) = g(x,y,z)*(\dt f(t)) + (-\Lap g(x,y,z))*f(t) = rhs
%
% In this case the source term is separable aswell and of the kind:
%                     
%    rhs =  g1(x,y,z)*f1(t) + g2(x,y,z)*f2(t)
%
% Then we insert a flag
%
problem_data.is_separable_u = true; 
%
% that controls the separability.
%
if problem_data.is_separable_u 
  % if it is true you define 
  problem_data.ux  = @(x) (abs(x-1/2)).^alpha; % space component of solution
  problem_data.ut  = @(t) exp(-t);   % time component of solution
  % and for the right hand side
  problem_data.f1=@(t) -1*exp(-t);
  problem_data.g1=@(x) problem_data.ux(x);
  problem_data.f2=@(t) problem_data.ut(t);
  problem_data.g2=@(x) -alpha*(alpha-1)*(abs(x-1/2)).^(alpha-2);
else
  % write instead of 0 your non-separable right hand side:
  problem_data.f = @(x, t) -((abs(x-1/2)).^(alpha) + alpha*(alpha-1)*(abs(x-1/2)).^(alpha-2)).*exp(-t); 
end
% N.B. 
% if you want to use non-separable rhs format for separables rhs then use
%
% problem_data.f = @(x,y,z,t) problem_data.g1(x,y,z).*problem_data.f1(t) + ...
%                             problem_data.g2(x,y,z).*problem_data.f2(t); 
%
% -------------------------------------------------------------------------

% Dirichlet data 
problem_data.dfun = @(x, t, iside) problem_data.uex(x, t);
problem_data.ifun = @(x, iside) problem_data.uex(x, 0);

problem_data.gmm = 1; % parameter
problem_data.eta = 1; % parameter

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 

p = 6; % polynomial degree of spline spaces
i = 6; % number of dyadic refinements
Nt = 2^i; nel_i = Nt-p+2; nel_t = Nt-p+1;

method_data.trial_degree     = [p p];                         % Degree of the trial splines (last is time dir)
method_data.trial_regularity = method_data.trial_degree-1;    % Regularity of trial the splines
method_data.test_degree      = [p p];                         % Degree of the test splines (last is time dir)
method_data.test_regularity  = method_data.test_degree-1;     % Regularity of the test splines
method_data.nsub  = [nel_i nel_t];                            % Number of subdivisions
method_data.nquad = method_data.trial_degree+1;               % Points for the Gaussian quadrature rule

% list of methods: 
% 'Galerkin'
% 'Least-Squares' (NOT READY    )
%
method_data.method     = 'Galerkin';  

% list of solvers: 
% 'GMRES' = Generalized Minimal RESiduals without preconditioners.
%    'CG' = Conjugate Gradients without preconditioners.(4 symmetric cases)
%    'LU' = (p)GMRES solver with FD in space + LU decomposition in time.
%    'AR' = (p)GMRES solver with FD in space + ARROW decomposition in time.
%    'LR' = (p)GMRES solver with FD in space + Low-Rank decomposition in
%           time. It uses Sherman-Morrison-Woodbury formula for the
%           inversion of the time block factors.
% 
method_data.solver = 'LU'; 

% 3.0) SOLVE FIRST THE PROBLEM WITH CRANK-NICOLSON METHOD IN TIME _________
method_data.n_time_intervals = 9;
problem_data.f = @(x, t)   -((abs(x-1/2)).^(alpha) + alpha*(alpha-1)*(abs(x-1/2)).^(alpha-2)).*exp(-t);  
[geo, msh, space, u] = heat_crank_nicolson_solve (problem_data, method_data);
% plot of solution with crank nicolson
xp={linspace(0,1,50)}; % punti di valutazione in spazio
t={linspace(0,1,method_data.n_time_intervals+1)}; % punti di valutazione in tempo
lables_t = geo.tgeo.map(t);
esolx = []; euex  = [];
for tpnt = 1: numel(t{:})
    tt = lables_t(tpnt);
    euex = cat(3,euex, problem_data.uex(xp{:}, tt));
    [evs, F] = sp_eval (u(:,tpnt), space.xsp_trial, geo.xgeo,xp);
    esolx = cat(3,esolx,evs);
end
min_eu = min(esolx(:));
max_eu = max(esolx(:));
min_eexu = min(euex(:));
max_eexu = max(euex(:));
min_esol = min(min_eu,min_eexu);
max_esol = max(max_eu,max_eexu);

fig=figure();
% fig.WindowState = 'maximized';
for it=1:numel(t{:})
subplot(2,(method_data.n_time_intervals+1)/2,it) % per vedere in un unico frame.
plot(squeeze(F(1,:,1)),esolx(:,:,it),'LineWidth',1.5);
hold on 
plot(squeeze(F(1,:,1)),euex(:,:,it),':','LineWidth',1.5);
%axis tight; axis square;
title(sprintf('Heat solution at t=%.2f', lables_t(it)), 'FontSize', 10);
hold off 
ylim([min_esol,max_esol]);
legend('u_h','u')
%   pause(2)
end 
%__________________________________________________________________________

% 3) CALL TO THE SOLVER
[geo, msh, space, u, report] = heat_st_solve (problem_data, method_data);

report

% 4) VISUALIZE THE SOLUTION
 heat_st_plot(u, problem_data.uex, space, geo, 10)

% 5) COMPUTE THE ERRORS 
[errl2, errh1s, errh1t] = ... Asbsolute errors
    st_h1_error_tp(space.xsp_trial,  space.tsp_trial,...
                        msh.xmsh,  msh.tmsh,  u,  problem_data.uex, ...
                             problem_data.dt_uex,   problem_data.grad_uex);

[l2_uex, h1s_uex, h1t_uex] = ... L2-norm and h1-(semi)-norm of solution
    st_h1_error_tp(space.xsp_trial,  space.tsp_trial,...
                        msh.xmsh,  msh.tmsh, 0*u,  problem_data.uex, ...
                             problem_data.dt_uex,   problem_data.grad_uex);

report.rel_errl2  = errl2/l2_uex;
report.rel_errh1s = errh1s/h1s_uex;
report.rel_errh1t = errh1t/h1t_uex;

% format output fancy
labels = { ...
    '‖ u - u_{ex} ‖_{L²(Ωₛ × Ωₜ)}', ...
    '‖ ∇(u - u_{ex}) ‖_{L²(Ωₛ × Ωₜ)}', ...
    '‖ ∂ₜ(u - u_{ex}) ‖_{L²(Ωₛ × Ωₜ)}' ...
};

values = [report.rel_errl2, report.rel_errh1s, report.rel_errh1t]; 
fprintf('%-35s | %12s\n', 'Error type', 'Value'); 
fprintf('%s\n', repmat('-',1,50)); 

for k = 1:length(values) 
    fprintf('%-35s | %12s\n', labels{k}, values(k)); 
end
