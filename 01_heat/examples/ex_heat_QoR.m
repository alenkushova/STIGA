% PROBLEM_NAME: <description>
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

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
T = 1; problem_data.T = T ;                                              

% NO NEED OF SPACE-TIME DOMAIN.
% 2D domain
ring = geo_load('geo_ring.txt');
% 3D domain  
problem_data.x_geo_name = nrbrevolve(ring.nurbs, [-1 -1 0], [ 1 0 0], pi/2);  
problem_data.t_geo_name = nrbline ([0 0], [T 0]); % time domain (1D)

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 4 5 6];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 3 4 5 6];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Exact solution: 
problem_data.uex = @(x, y, z, t) -(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*sin(t).*sin(z); % separable

problem_data.graduex = @(x, y, z, t) (cat (1, ...
                reshape ( -2*x.*(x.^2+y.^2-4).*x.*y.^2.*sin(t).*sin(z) ...
                          -2*(x.^2+y.^2-1).*x.^2.*y.^2.*sin(t).*sin(z) ...
                          -(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*sin(t).*sin(z), [1, size(x)]), ...
                reshape ( -2*y.*(x.^2+y.^2-4).*x.*y.^2.*sin(t).*sin(z) ...
                          -2*(x.^2+y.^2-1).*x.*y.^3.*sin(t).*sin(z) ...
                          -2*(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.*sin(t).*sin(z), [1, size(y)]), ...
                reshape ( -(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*sin(t).*cos(z), [1, size(z)]), ...
                reshape ( -(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*cos(t).*sin(z), [1, size(t)])));

% -------------------------------------------------------------------------
% Let  u(x,y,z,t) = g(x,y,z)*f(t) be separable, then:
%
%     Heat u(x,y,z,t) = g(x,y,z)*(\dt f(t)) + (-\Lap g(x,y,z))*f(t) = rhs
%
% In this case the source term is separable aswell and of the kind:
%                     
%     rhs =  g1(x,y,z)*f1(t) - g2(x,y,z)*f2(t)
%
% Then we insert a flag
%
problem_data.is_separable_u = true; 
%
% that controls the separability.
%
if problem_data.is_separable_u 
  % if it is true you define 
  problem_data.ux  = @(x, y, z) -x.*y.^2.*(x.^2+y.^2-1).*(x.^2+y.^2-4).*sin(z); % space component of solution
  problem_data.ut  = @(t) sin(t);                                               % time component of solution
  % and for the right hand side
  problem_data.f1=@(t) cos(t);
  problem_data.g1=@(x, y, z) problem_data.ux(x, y, z);
  problem_data.f2=@(t) problem_data.ut(t);
  problem_data.g2=@(x,y,z) x.*sin(z).*(40*y.^2.*(x.^2+y.^2-2.0) + (2-y.^2).*(x.^2+y.^2-1).*(x.^2+y.^2-4));
else
  % write instead of 0 your non-separable right hand side:
  problem_data.f = @(x,y,z,t)     0  ; 
end
% N.B. 
% if you want to use non-separable rhs format for separables rhs then use
%
% problem_data.f = @(x,y,z,t) problem_data.g1(x,y,z).*problem_data.f1(t) + ...
%                             problem_data.g2(x,y,z).*problem_data.f2(t); 
%
% -------------------------------------------------------------------------
% Dirichlet data 
problem_data.dfun = @(x, y, z, t, iside) problem_data.uex(x, y, z, t);
problem_data.ifun = @(x, y, z, iside) problem_data.uex(z, y, z, 0);
problem_data.gmm = 1;
problem_data.eta = 1; 

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 

p = 2; % polynomial degree of spline spaces
i = 2; % number of dyadic refinements
Nt = 2^i; nel_i = Nt-p+2; nel_t = Nt-p+1;

method_data.trial_degree     = [p p p p];                                  % Degree of the trial splines (last is time dir)
method_data.trial_regularity = method_data.trial_degree-1;                 % Regularity of trial the splines
method_data.test_degree      = [p p p p];                                  % Degree of the test splines (last is time dir)
method_data.test_regularity  = method_data.test_degree-1;                  % Regularity of the test splines
method_data.nsub  = [nel_i nel_i nel_i nel_t];                        % Number of subdivisions
method_data.nquad = method_data.trial_degree+1;                       % Points for the Gaussian quadrature rule

% list of methods: 
method_data.method     = 'Galerkin';  

% list of solvers: 
method_data.solver = 'GMRES'; 

% 3) CALL TO THE SOLVER
[geo, msh, space, u, report] = heat_st_solve (problem_data, method_data);

report