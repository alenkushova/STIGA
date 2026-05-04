% STOKES_ST_SOLVE: <description>
% 
%   CALL: 
%
%  INPUT:
%
% OUTPUT:
%
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

function [geo, msh, space, vel, pres, report] = stokes_st_solve(problem_data,method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

[geo, msh, space] = discretize_stokes(problem_data,method_data);
[vel, pres, report] = stokes_st_problem(msh, space, problem_data, method_data);

end

 

function  [geo, msh, space] = discretize_stokes(problem_data,method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

tgeo  = geo_load(geo_time); 
xgeo  = geo_load(geo_space);   
xtgeo = geo_load(geo_space_time);   

geo.tgeo = tgeo; geo.xgeo = xgeo; geo.xtgeo = xtgeo;

[~, zeta] = kntrefine (xtgeo.nurbs.knots, nsub-1, trial_degree, trial_regularity);
rule      = msh_gauss_nodes (nquad);
[qn, qw]  = msh_set_quad_nodes (zeta, rule);

tmsh      = msh_cartesian (zeta(end), qn(end), qw(end), tgeo);
xmsh      = msh_cartesian (zeta(1:end-1), qn(1:end-1), qw(1:end-1), xgeo);
xtmsh     = msh_cartesian (zeta, qn, qw, xtgeo);

msh.tmsh = tmsh; msh.xmsh = xmsh; msh.xtmsh = xtmsh;

% Let us split the degree regularity and number of subdivision variables
% into space and time. 
degree_ps = trial_degree(1:end-1);          degree_pt = trial_degree(end);
regularity_ps = trial_regularity(1:end-1);  regularity_pt = trial_regularity(end);
nsub_ps = nsub(1:end-1);              nsub_pt = nsub(end);

% And the spline spaces in time direction of degree pt (both for pres/vel)
knots_pt = kntrefine (tgeo.nurbs.knots, nsub_pt-1, degree_pt, regularity_pt);
spt_pt = sp_bspline (knots_pt, degree_pt, tmsh);

% Then we may also need a spline space in time with degree pt+1               
% ( we use ptp1 for varibale names )
knots_ptp1 = kntrefine (tgeo.nurbs.knots, nsub_pt-1, degree_pt+1, regularity_pt+1);
spt_ptp1  = sp_bspline (knots_ptp1, degree_pt + 1, tmsh);

% This is a black box construction already implemented in GeoPDEs
[sps_v, sps_ps] = sp_bspline_fluid ('TH', xgeo.nurbs.knots, nsub_ps, degree_ps, regularity_ps, msh.xmsh);

% save also the space structure
space.spt_vel = spt_pt; 
space.spt_pres = spt_pt;
space.spp = sps_ps; 
space.spv = sps_v; 

end

function  [vel, pres, report] = stokes_st_problem(msh, space, problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

dim = space.spv.ncomp;

sizev = space.spv.ndof*space.spt_vel.ndof;
sizep = space.spp.ndof*space.spt_vel.ndof;

vel  = zeros(sizev,1);  

fprintf('Projecting boundary conditions... \n\n') 
% we impose boundary conditions and here we deal initial conditions weakly
[vel_drchlt, full_drchlt_dofs, vel_iniz] = st_stokes_boundary_data(space.spv, space.spt_vel, msh.xmsh, msh.tmsh, dfun, drchlt_sides,'yes');
vel(full_drchlt_dofs) = vel_drchlt;
% notice 'yes'  is an optional variable (default 'no'). It referst to the
% question: Is the inital condition applied weakly? Usually no, here yes!
% If yes, then the first face of the ST-cuboid is to be considered part
% Dirichlet boundary and part Internal dof. If no, the whole first face
% is of Dirichlet boundary.

% Now we compute only the space-boundary degrees of freedom
if dim == 2
[~, x_drchlt_dofs] = sp_drchlt_l2_proj (space.spv, msh.xmsh, @(x, y, iside) velex(x,y,0), drchlt_sides);
elseif dim == 3
[~, x_drchlt_dofs] = sp_drchlt_l2_proj (space.spv, msh.xmsh, @(x, y, z, iside) velex(x,y,z,0), drchlt_sides);
end
x_int_dofs = setdiff (1:space.spv.ndof, x_drchlt_dofs); 
int_dofs = setdiff(1:sizev,full_drchlt_dofs);

% Dim. of the vectors in space (or time) only (internal or Dirichlet dofs)
intnx = numel(x_int_dofs); 
drchnx = numel(x_drchlt_dofs); 
intnt = space.spt_vel.ndof;
% for both space and time here is the total number of internal dofs
nintdofs = numel(int_dofs);
presnx = space.spp.ndof;

fprintf('Assembling space factor matrices... \n\n') 
Kx = op_gradu_gradv_tp (space.spv, space.spv, msh.xmsh, viscosity); 
Mx = op_u_v_tp (space.spv, space.spv, msh.xmsh); 
Dx = op_div_v_q_tp (space.spv, space.spp, msh.xmsh);
Mxp= op_u_v_tp (space.spp, space.spp, msh.xmsh); 
Kx = (Kx+Kx')/2;  Mx = (Mx+Mx')/2;  Mxp = (Mxp+Mxp')/2;

fprintf('Assembling time factor matrices... \n\n')
Mt = op_u_v_tp (space.spt_pres, space.spt_pres, msh.tmsh);
Wt = op_vel_dot_gradu_v_tp (space.spt_pres, space.spt_pres, msh.tmsh, @(x) ones(size(x))); % (dt u, v)
Mt = (Mt+Mt')/2;
% Since the initial data are imposed weakly I should add kron(At,Mx) with 
%      |1, 0, ..., 0|   
% At = |0, ...   , 0| square matrix of size(Mt) with only 1 in first entry, 
%      |            |
%      |0, ...   , 0| 
% but this can be done as kron(Wt+At, Mx), and so directly as: 
Wt(1,1) = Wt(1,1) + 1; % since all the rest is equal to zero. 

% The unique solvability is guaranteed despite constant shifts for the 
% pressures. We therefore want to impose 0 mean value to the pressures.
% For this reason we consider the following:
fun_one = @(varargin) ones(size(varargin{1})); % function of all ones
means = op_f_v_tp (space.spp, msh.xmsh, fun_one).'; % \int {pres.basis} \d \Omega
% Notice that it is stored as a horizontal vector!

fprintf('Assembling right hand side f... \n\n')
if exist('f','var')
  rhs_vel = op_f_v_st_tp (space.spv, space.spt_vel, msh.xmsh, msh.tmsh, f);
else
  rhs_vel = zeros (sizev, 1);
end        

% The matrix A that we want to invert is: 
%        _________________________
%       |            |       |    |
%       |            |       |    |
%       |    H       |  BT   |    |
%       |            |       |    |
% A =   | ___________|_______|____|
%       |            |       |    |
%       |   B        |       | CT |
%       |____________|_______|____|
%       |            |   C   |    |
%       |____________|_______|____|
%
% It applies to vectors:  x = ( vel. , pres. , mean_pres. )'
% So we need the following matrix vector multiplication: 
%
%   H , B : they must be applied to velocity fields
%   H , B : they must have an application to dirichlet data of velocities
%   BT, C : they must be applied to pressure fields
%   CT    : it must be applied to the mean values of the pressures 
%


Hdrchl2int = @(x) reshape(Mx(x_int_dofs,x_drchlt_dofs) *reshape(x,drchnx,intnt)*(Wt') + ...
                          Kx(x_int_dofs,x_drchlt_dofs) *reshape(x,drchnx,intnt)*(Mt'),  ...
                          nintdofs,1);
Bdrchl2int = @(x) -reshape(Dx(:,x_drchlt_dofs) *reshape(x,drchnx,intnt)*(Mt'), ...
                          presnx*intnt,1);

Hfun = @(x) reshape(Mx(x_int_dofs,x_int_dofs) *reshape(x,intnx,intnt)*(Wt') + ...
                    Kx(x_int_dofs,x_int_dofs) *reshape(x,intnx,intnt)*(Mt'),  ...
                    nintdofs,1);

Bfun = @(x) -reshape(Dx(:,x_int_dofs) *reshape(x,intnx,intnt)*(Mt'), ...
                     presnx*intnt,1);

BTfun = @(x) -reshape((Dx(:,x_int_dofs))'*reshape(x,presnx,intnt)*Mt, ...
                     nintdofs,1);

Cfun  = @(x) reshape( means * reshape(x,presnx,intnt) * Mt', intnt, 1 );

CTfun = @(x) reshape(means'*reshape(x,1,intnt)*Mt, presnx*intnt, 1);

% Finally we can define Afun: 
Afun = @(x) cat(1, Hfun(x(1:nintdofs)) + BTfun(x(nintdofs+1:nintdofs+sizep)),...
                   Bfun(x(1:nintdofs)) + CTfun(x(nintdofs+sizep+1:end)),...
                   Cfun(x(nintdofs+1:nintdofs+sizep)));

% eliminiamo i prodotti kronecker.... Quindi non abbiamo A e B. 
vel_iniz = kron([1;zeros(space.spt_vel.ndof-1,1)], vel_iniz);
rhs_vel  = rhs_vel(int_dofs) + vel_iniz(int_dofs) - Hdrchl2int(vel(full_drchlt_dofs)); 
rhs_pres = -Bdrchl2int(vel(full_drchlt_dofs));
rhs = [rhs_vel ; rhs_pres ; zeros(intnt,1)];

fprintf('Assembling Arrow (AR) preconditioner for velocity block... \n\n')
parametric_problem_data = problem_data; 
square = nrbextrude( nrbline ([0 0], [1 0]), [0,1]); % square NURBS surface
parametric_problem_data.geo_space = 'geo_square.txt'; % square NURBS surface as .txt
parametric_problem_data.geo_space_time = nrbextrude(square, [0 0 1]); % NURBS volume

% define parametric spaces
[Pgeo, Pmsh, Pspace] = discretize_stokes(problem_data,method_data);
% Pspace contains paremtric spaces velocity and pressure (and time space)
% Pspace.spv = spazio spline per le velocità sul dominio parametrico.
% Pspace.spp = spazio spline per le pressioni sul dominio parametrico.
% Pspace.spt_vel  = spazio spline in tempo per le velocità
% ""      ""_pres = spazio spline in tempo per le pressioni ma è lo stesso
%                   spazio in tempo delle velocità.
cell_geo = cell(1, dim + 2); 
cell_mesh = cell(1, dim + 2);
cell_space_1 = cell(2, dim + 2);
cell_space_2 = cell(2, dim + 2);

for i = 2:dim+1
  j = i-1;
  cell_geo{i}   = geo_load(nrbline ([0 0], [1 0]));
  rule     = msh_gauss_nodes (nquad(j));
  [qn, qw] = msh_set_quad_nodes (Pmsh.xmsh.breaks(j), rule);
  cell_mesh{i}   = msh_cartesian (Pmsh.xmsh.breaks(j), qn, qw, cell_geo{i});
  cell_space_1{1,i} = Pspace.spv.scalar_spaces{1}.sp_univ(j);
  cell_space_1{2,i} = cell_space_1{1,i}; 
  cell_space_2{1,i} = Pspace.spv.scalar_spaces{2}.sp_univ(j);
  cell_space_2{2,i} = cell_space_2{1,i}; 
end

cell_geo{end} = Pgeo.tgeo;
cell_mesh{end} = Pmsh.tmsh;
cell_space_1{1,end} = Pspace.spt_vel;
cell_space_1{2,end} = cell_space_1{1,end};
cell_space_2{1,end} = Pspace.spt_vel;
cell_space_2{2,end} = cell_space_2{1,end};

% arrow preconditioner for first compontent of velocity field
varout1 = generate_heat_pencils(dim, cell_mesh, cell_space_1);
P1 = arrow_heat_setup(varout1); 
% arrow preconditioner for second compontent of velocity field 
varout2 = generate_heat_pencils(dim, cell_mesh, cell_space_2);
P2 = arrow_heat_setup(varout2); 

% MASS = 10;%

Pfun = @(x) cat(1, P1(x(1:nintdofs/2)), ...
                   P2(x(nintdofs/2+1,nintdofs))); %,...
                  % MASS(x(nintdofs+1:nintdofs+sizep)),...
                  % x(nintdofs+sizep+1:end));

[sol, flag, rel_res, iter, res_vec] = gmres(Hfun,rhs(1:nintdofs),[],1e-8, 200, @(x) Pfun(x));


fprintf('Solving the linear system with GMRES... \n\n')
[sol, flag, rel_res, iter, res_vec] = gmres(Afun, rhs, [], 1e-8, 1000);
report.flag    = flag;
report.rel_res = rel_res;
report.iter    = iter;
report.res_vec = res_vec;
fprintf('Done. \n\n')

vel(int_dofs) = sol(1:nintdofs);
pres = sol(nintdofs+1:nintdofs+sizep);

% The preconditioner is block diagonal given by 
%        _________________________
%       |            |       |    |
%       |            |       |    |
%       |    H       |       |    |
%       |            |       |    |
% P =   | ___________|_______|____|
%       |            |       |    |
%       |            | MASS  |    |
%       |____________|_______|____|
%       |            |       |1   |
%       |____________|_______|___1|
%
% we define its forward application as follows 
% P = @(x) cat(1, inv_Hfun (x(1:sizev)) ,...
%                 inv_Mfun (x(sizev+1:sizev+sizep)),...
%                 IDfun (x(sizev+sizep+1:end)));

end
