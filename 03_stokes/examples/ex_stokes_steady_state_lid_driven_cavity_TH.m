% EX_STOKES_SQUARE_TH: solve the Stokes problem in the unit square with generalized Taylor-Hood elements.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:4;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
fx = @(x, y) 0*x;
fy = @(x, y) 0*x;
problem_data.f  = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Boundary terms
problem_data.h = @(x, y, iside) lid_driven_dfun(x, y, iside);

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'th';   % Element type for discretization
method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
method_data.nsub         = [10 10];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes (problem_data, method_data);

% 4) POST-PROCESSING
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};

[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure()
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
hold on 
contourf(X, Y, sqrt(squeeze(eu(1,:,:)).^2 + squeeze(eu(2,:,:)).^2));

title('Computed solution')


function h  = lid_driven_dfun (x, y, iside) 
 switch (iside)
  case {1, 2, 3}
   h = zeros ([2, size(x)]);
  case 4
   h1 = ones ([1, size(x)]);
   h2 = zeros ([1, size(x)]);
   h  = cat (1, h1, h2);
 end
end
