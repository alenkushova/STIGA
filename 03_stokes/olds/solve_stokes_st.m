function [geo, msh, space, vel, pres, report] = solve_stokes_st(problem_data, method_data)
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

[~, zeta] = kntrefine (xtgeo.nurbs.knots, nsub-1, degree, regularity);
rule      = msh_gauss_nodes (nquad);
[qn, qw]  = msh_set_quad_nodes (zeta, rule);

tmsh      = msh_cartesian (zeta(end), qn(end), qw(end), tgeo);
xmsh      = msh_cartesian (zeta(1:end-1), qn(1:end-1), qw(1:end-1), xgeo);
xtmsh     = msh_cartesian (zeta, qn, qw, xtgeo);

msh.tmsh = tmsh; msh.xmsh = xmsh; msh.xtmsh = xtmsh;

knots_pt = kntrefine (tgeo.nurbs.knots, nsub(end)-1, degree(end), regularity(end));
knots_ptp1 = kntrefine (tgeo.nurbs.knots, nsub(end)-1, degree(end)+1, regularity(end)+1);
spt_p = sp_bspline (knots_pt, degree(end), tmsh);
spt_pp1  = sp_bspline (knots_ptp1, degree(end)+1, tmsh);

degree_p = degree(1:end-1);
regularity_p = regularity(1:end-1);
nsub_p = nsub(1:end-1);
knots_p = kntrefine (xgeo.nurbs.knots, nsub_p-1, degree_p, regularity_p);
spp = sp_bspline (knots_p, degree(1:end-1), xmsh);
xtspp = sp_bspline ([knots_p(:)' {knots_pt}], degree, xtmsh);

degree_v = degree_p + 1; 
regularity_v = regularity_p;
nsub_v = nsub_p;
knots_v = kntrefine (xgeo.nurbs.knots, nsub_v-1, degree_v, regularity_v);

scalar_space = sp_bspline (knots_v, degree_v, xmsh);
for idim = 1:xmsh.ndim
    scalar_spaces{idim} = scalar_space;
end
spv = sp_vector (scalar_spaces, xmsh);

space.spt_vel = spt_p; 
space.spt_pres = spt_p;
space.spp = spp; 
space.spv = spv; 

fprintf('Assembling space factor matrices... \n\n') 
tic;
Kx = op_gradu_gradv_tp (spv, spv, xmsh, viscosity); 
Mx = op_u_v_tp (spv, spv, xmsh); 
Dx = op_div_v_q_tp (spv, spp, xmsh);
Mxp= op_u_v_tp (spp, spp, xmsh); 
report.space_factors = toc;
fprintf('Assembling time factor matrices... \n\n')
tic; 
Mt = op_u_v_tp (spt_p, spt_p, tmsh);
Wt = op_vel_dot_gradu_v_tp (spt_p, spt_p, tmsh, @(x) ones(size(x)));% (dt u, v)
report.time_factors = toc;

% a  = [1; zeros(size(Mt,1)-1,1)]; 
% At = spdiags(a,0,numel(a),numel(a));
% e ora vorrei aggiungere kron(At,Mx), quindi posso fare direttamente 
% kron(Wt+At,Mx) ma allora posso prima prendere questo: 
Wt(1,1) = Wt(1,1)+1; % perchè tutto il resto fa zero!
fprintf('Assembling A... \n\n')
tic;
A  = kron(Mt,Kx) + kron(Wt,Mx); % matrice totale assemblata. 
report.matrix_A = toc;
fprintf('Assembling B... \n\n')
tic;
B  = kron(Mt,Dx); B = -B; % pressures = 0 at t=0.
report.matrix_B = toc;

vel  = zeros(spv.ndof*spt_p.ndof,1); 
full_drchlt_dofs = [];
% project initial data and dirichlet data 
if exist('ifun','var')
  vel = reshape(vel,spv.ndof,spt_p.ndof);
  [vel_drchlt, drchlt_dofs, vel_iniz] = stokes_boundary_data(spv, xmsh, dfun, drchlt_sides, ifun);
  for idof = 1:spt_p.ndof
        vel(drchlt_dofs, idof) = vel_drchlt; % imponiamo le condizioni di Dirichlet 
        full_drchlt_dofs = cat(1, full_drchlt_dofs, drchlt_dofs + (idof-1)*(spv.ndof));
  end
  vel = vel(:);
else 
  tic;  
  [vel_drchlt, full_drchlt_dofs, vel_iniz] = st_stokes_boundary_data(spv, spt_p, xmsh, tmsh, dfun, drchlt_sides,'yes');
  report.boundary_data = toc;
  % notice 'yes'  is an optional variable (default 'no'). It referst to the
  % question: Is the inital condition applied weakly? Usually no, here yes!
  % If yes, then the first face of the ST-cuboid is to be considered part
  % Dirichlet boundary and part Internal dof. If no, the whole first face
  % is of Dirichlet boundary.
  vel(full_drchlt_dofs) = vel_drchlt;
end

int_dofs = setdiff(1:spv.ndof*spt_p.ndof,full_drchlt_dofs);
nintdofs = numel(int_dofs);

if exist('f','var')
    fprintf('Assembling right hand side f... \n\n')
    tic;
    rhs_vel = op_f_v_st_tp(spv,spt_p,xmsh,tmsh,f);
    report.vector_f = toc;
else
    rhs_vel = zeros(spv.ndof*spt_p.ndof,1);
end        

vel_iniz = kron([1;zeros(spt_p.ndof-1,1)], vel_iniz);

rhs_vel  = rhs_vel(int_dofs) + vel_iniz(int_dofs) -A(int_dofs, full_drchlt_dofs)*vel(full_drchlt_dofs); 
rhs_pres = -B(:, full_drchlt_dofs)*vel(full_drchlt_dofs);
rhs = [rhs_vel ; rhs_pres];

fun_one = @(varargin) ones(size(varargin{1}));
E = op_f_v_tp (spp, xmsh, fun_one).';       

fprintf('Assembling the whole space time matrix + constraint on pressure average... \n\n')
tic;
mat = [A(int_dofs, int_dofs)       B(:,int_dofs)'               sparse(nintdofs,size(Mt,1)) ; ...
       B(:,int_dofs)               sparse(size(B,1),size(B,1))  (kron(Mt,E))'; ...
       sparse(size(Mt,1),nintdofs) kron(Mt,E)                   sparse(size(Mt,1),size(Mt,1)) ];
report.whole_matrix = toc;

P = blkdiag(A(int_dofs,int_dofs),kron(Mt,Mxp),speye(size(Mt)));

%fprintf('Solving the linear system with matrix backslash... \n\n')
tic;
%sol = mat\[rhs;zeros(spt_p.ndof,1)];
fprintf('Solving the linear system with preconditioned GMRES... \n\n')
[sol, flag, rel_res, iter, res_vec] = gmres(mat,[rhs;zeros(spt_p.ndof,1)],[],1e-8,200,P);
report.solving_time = toc;
report.solution_details.flag = flag;
report.solution_details.rel_res = rel_res;
report.solution_details.iter = iter;
report.solution_details.res_vec = res_vec;

fprintf('Done. \n\n')
vel(int_dofs) = sol(1:nintdofs);
pres = sol(1+nintdofs:end-spt_p.ndof);  

% questo solo per dubug rhs space-time.
% Mp =  op_u_v_tp (spp, spp, xmsh);
% rhs_pres = op_f_v_tp(xtspp,xtmsh,presex);
% 
% rhs_pres_st = op_f_v_st_tp(spp,spt_p,xmsh,tmsh,presex);
% pres = kron(Mt,Mp)\ rhs_pres;
% 
end

