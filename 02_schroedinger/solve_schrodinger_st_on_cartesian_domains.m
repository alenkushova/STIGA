function [geo, msh, space, u, report] = solve_schrodinger_st_on_cartesian_domains (problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Univariate discrete spaces for Cartesian domains:
[geo, msh, space] = univariate_spaces(problem_data,method_data);



% Assembling the linear system and the right hand side over Cartesian
% domains and using full properties of tensorization among each direction
[Aint, rhs, u, int_dofs, Afun, mat] = schrodinger_st_problem_on_cartesian_domains(msh, space, problem_data, method_data);
%end

%force Matlab to use one single thread:
N = maxNumCompThreads; 
maxNumCompThreads(1);

% Solving: 
switch solver
  case 'MB'
    % Solution with backslash
    fprintf('Solving with MATLAB backslash... \n\n')
    tic;
    u(int_dofs)  = Aint\rhs;
    time = toc;
    report.flag = 'Solved with MATLAB backslash';
    report.time = time;
  case 'CG'
    tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
    % Solution with CG
    fprintf('Solving with CG (without prec.)... \n\n')
    tic
    [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit);
    time = toc;
    u(int_dofs) = u_inner;
    % Output informations of CG solver iterations
    report.flag    = flag;
    report.rel_res = rel_res;
    report.iter    = iter;
    report.res_vec = res_vec;
    report.time    = time;
  case 'PCG'
    switch preconditioner 
        case 'LUFD' % solving with PCG and block LU with FD in space
          fprintf('Assembling LU-block preconditioner.. \n\n ')
          tol = 10^(-8); maxit = min(numel(int_dofs),200);   
          switch space_dimension
            case '1D'
              tic
              P = LU_SETUP_x_SCHRODINGER_1D(mat.Lt,mat.Mt,mat.Wt,mat.L1,mat.M1);
              fprintf('Solving with PCG and LU-block preconditioner... \n\n')
              [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit,@(x) P(x));
              time = toc;
            case '2D'
              tic
              P = LU_SETUP_x_SCHRODINGER_2D(mat.Lt,mat.Mt,mat.Wt,mat.L1,mat.M1,mat.L2,mat.M2);
              fprintf('Solving with PCG and LU-block preconditioner... \n\n')
              [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit,@(x) P(x));
              time = toc;
            case '3D'
              tic
              P = LU_SETUP_x_SCHRODINGER_3D(mat.Lt,mat.Mt,mat.Wt,mat.L1,mat.M1,mat.L2,mat.M2,mat.L3,mat.M3);
              fprintf('Solving with PCG and LU-block preconditioner... \n\n')
              [u_inner, flag, rel_res, iter, res_vec] = pcg(Afun,rhs,tol,maxit,@(x) P(x));
              time = toc;
          end
          u(int_dofs) = u_inner;
          % Output informations of GMRES solver iterations
          report.flag    = flag;
          report.rel_res = rel_res;
          report.iter    = iter;
          report.res_vec = res_vec;
          report.time    = time;
      case 'ilu'
        tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
        fprintf('Solving with CG and ILU preconditioner... \n\n')
        tic
        [L,U] = ilu(Aint); 
        [u_inner, flag, rel_res, iter, res_vec] = pcg(Aint,rhs,tol,maxit,L,U);
        time = toc;
        u(int_dofs) = u_inner;
        % Output informations of CG solver iterations
        report.flag    = flag;
        report.rel_res = rel_res;
        report.iter    = iter;
        report.res_vec = res_vec;
        report.time    = time;
      case 'ichol'
        tol   = 10^(-8);  maxit = min(numel(int_dofs),200);  
        fprintf('Solving with CG and ICHOL preconditioner... \n\n')
        tic
        L = ichol(Aint); 
        [u_inner, flag, rel_res, iter, res_vec] = pcg(Aint,rhs,tol,maxit,L,L');
        time = toc;
        u(int_dofs) = u_inner;
        % Output informations of CG solver iterations
        report.flag    = flag;
        report.rel_res = rel_res;
        report.iter    = iter;
        report.res_vec = res_vec;
        report.time    = time;
    end
end

maxNumCompThreads(N);

end

function [geo, msh, space] = discretize(problem_data, method_data)
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% load geometies
xtgeo = geo_load(xt_geo_name);
xgeo  = geo_load(x_geo_name);
tgeo  = geo_load(t_geo_name);

geo = struct('xgeo',xgeo,'tgeo',tgeo,'xtgeo',xtgeo);

% define mesh structures :
rdim = numel(trial_degree);
% space time knots for trial functions
[knots, zeta]    = kntrefine(xtgeo.nurbs.knots, nsub-1, trial_degree, trial_regularity);
knots   = kntunclamp(knots, trial_degree, trial_regularity, []);

% space knots for trial functions in space
[x_knots, x_zeta]= kntrefine(xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, trial_degree(1:rdim-1), trial_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, trial_degree(1:rdim-1), trial_regularity(1:rdim-1), []);

%time knots for trial functions in time
[t_knots, t_zeta]= kntrefine(tgeo.nurbs.knots,...
                    nsub(rdim)-1, trial_degree(rdim), trial_regularity(rdim));
t_knots = kntunclamp(t_knots, trial_degree(rdim), trial_regularity(rdim), []);

% define quadrature rules 
rule         = msh_gauss_nodes (nquad);
[qn, qw]     = msh_set_quad_nodes (zeta, rule);
xtmsh  = msh_cartesian (zeta, qn, qw, xtgeo);

[xqn, xqw]   = msh_set_quad_nodes (x_zeta, rule(1:end-1));
xmsh = msh_cartesian (x_zeta, xqn, xqw, xgeo);

[tqn, tqw]   = msh_set_quad_nodes (t_zeta, rule(end));
tmsh = msh_cartesian (t_zeta, tqn, tqw, tgeo);

msh = struct('xmsh',xmsh,'tmsh',tmsh,'xtmsh',xtmsh);

% define space structures for trial functions
xtsp_trial = sp_bspline (knots,   trial_degree,          xtmsh);
xsp_trial  = sp_bspline (x_knots, trial_degree(1:rdim-1), xmsh);
tsp_trial  = sp_bspline (t_knots, trial_degree(rdim),     tmsh);

% space time knots for trial functions
[knots, ~]    = kntrefine(xtgeo.nurbs.knots, nsub-1, test_degree, test_regularity);
knots   = kntunclamp(knots, test_degree, test_regularity, []);

% space knots for test functions in space
[x_knots, ~]= kntrefine(xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, test_degree(1:rdim-1), test_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, test_degree(1:rdim-1), test_regularity(1:rdim-1), []);

%time knots for test functions in time
[t_knots, ~]= kntrefine(tgeo.nurbs.knots,...
                    nsub(rdim)-1, test_degree(rdim), test_regularity(rdim));
t_knots = kntunclamp(t_knots, test_degree(rdim), test_regularity(rdim), []);

% define space structures for test functions
xtsp_test = sp_bspline (knots,   test_degree,          xtmsh);
xsp_test  = sp_bspline (x_knots, test_degree(1:rdim-1), xmsh);
tsp_test  = sp_bspline (t_knots, test_degree(rdim),     tmsh);

space = struct('xsp_trial',xsp_trial,'tsp_trial',tsp_trial,'xtsp_trial',xtsp_trial,...
               'xsp_test' ,xsp_test ,'tsp_test' ,tsp_test ,'xtsp_test' ,xtsp_test);
end


function [geo,msh,space] = univariate_spaces(problem_data,method_data)
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% univariate directions: 
rdim = numel(trial_degree);

tgeo  = geo_load(t_geo_name);  %  t_geo_name = nrbline ([0 0], [T 0]); 
x1geo = geo_load(x1_geo_name); % x1_geo_name = nrbline ([0 0], [1 0]);
%xgeo  = geo_load(x_geo_name);

% knots:
[t_knots, t_zeta] = kntrefine(tgeo.nurbs.knots,...
                    nsub(rdim)-1, trial_degree(rdim), trial_regularity(rdim));
t_knots = kntunclamp(t_knots, trial_degree(rdim), trial_regularity(rdim), []);

[x1_knots, x1_zeta] = kntrefine(x1geo.nurbs.knots,...
                    nsub(1)-1, trial_degree(1), trial_regularity(1));
x1_knots = kntunclamp(x1_knots, trial_degree(1), trial_regularity(1), []);

% [x_knots, x_zeta] = kntrefine(xgeo.nurbs.knots,...
%                     nsub(1:rdim-1)-1, trial_degree(1:rdim-1), trial_regularity(1:rdim-1));
% x_knots = kntunclamp(x_knots, trial_degree(1:rdim-1), trial_regularity(1:rdim-1), []);

% Define quadrature rules
rule         = msh_gauss_nodes (nquad);

[tqn, tqw]   = msh_set_quad_nodes (t_zeta, rule(end));
tmsh = msh_cartesian (t_zeta, tqn, tqw, tgeo);

[x1qn, x1qw]   = msh_set_quad_nodes (x1_zeta, rule(1));
x1msh = msh_cartesian (x1_zeta, x1qn, x1qw, x1geo);

% [xqn, xqw]   = msh_set_quad_nodes (x_zeta, rule(1:end-1));
% xmsh = msh_cartesian (x_zeta, xqn, xqw, xgeo);

% Define trial spaces
tsp_trial   = sp_bspline (t_knots,  trial_degree(rdim),   tmsh);
x1sp_trial  = sp_bspline (x1_knots, trial_degree(1),     x1msh);
%xsp_trial   = sp_bspline (x_knots,  trial_degree(1:rdim-1),xmsh);

if (rdim == 2)
    % geo = struct('x1geo',x1geo,'tgeo',tgeo,'xgeo',xgeo);
    % msh = struct('x1msh',x1msh,'tmsh',tmsh,'xmsh',xmsh);
    % space = struct('x1sp_trial',x1sp_trial,'tsp_trial',tsp_trial,'xsp_trial',xsp_trial);
    geo = struct('x1geo',x1geo,'tgeo',tgeo);
    msh = struct('x1msh',x1msh,'tmsh',tmsh);
    space = struct('x1sp_trial',x1sp_trial,'tsp_trial',tsp_trial);
elseif (rdim == 3)
    x2geo = geo_load(x2_geo_name); % x2_geo_name = nrbline ([0 0], [1 0]);
    geo = struct('x1geo',x1geo,'x2geo',x2geo,'tgeo',tgeo);
    
    % knots:
    [x2_knots, x2_zeta] = kntrefine(x2geo.nurbs.knots,...
                        nsub(2)-1, trial_degree(2), trial_regularity(2));
    x2_knots = kntunclamp(x2_knots, trial_degree(2), trial_regularity(2), []);
            
    % Define quadrature rules
    rule         = msh_gauss_nodes (nquad);
    
    [x2qn, x2qw]   = msh_set_quad_nodes (x2_zeta, rule(2));
    x2msh = msh_cartesian (x2_zeta, x2qn, x2qw, x2geo);
    
    msh = struct('x1msh',x1msh,'x2msh',x2msh,'tmsh',tmsh);

    % Define trial spaces
    x2sp_trial  = sp_bspline (x2_knots, trial_degree(2),  x2msh);
    
    space = struct('x1sp_trial',x1sp_trial,'x2sp_trial',x2sp_trial,'tsp_trial',tsp_trial);


elseif (rdim == 4)
    x2geo = geo_load(x2_geo_name); % x2_geo_name = nrbline ([0 0], [1 0]);
    x3geo = geo_load(x3_geo_name); % x3_geo_name = nrbline ([0 0], [1 0]);
    geo = struct('x1geo',x1geo,'x2geo',x2geo,'x3geo',x3geo,'tgeo',tgeo);
    
    % knots:
    [x2_knots, x2_zeta] = kntrefine(x2geo.nurbs.knots,...
                        nsub(2)-1, trial_degree(2), trial_regularity(2));
    x2_knots = kntunclamp(x2_knots, trial_degree(2), trial_regularity(2), []);
        
    [x3_knots, x3_zeta] = kntrefine(x3geo.nurbs.knots,...
                        nsub(3)-1, trial_degree(3), trial_regularity(3));
    x3_knots = kntunclamp(x3_knots, trial_degree(3), trial_regularity(3), []);
    
    % Define quadrature rules
    rule         = msh_gauss_nodes (nquad);
    
    [x2qn, x2qw]   = msh_set_quad_nodes (x2_zeta, rule(2));
    x2msh = msh_cartesian (x2_zeta, x2qn, x2qw, x2geo);
    
    [x3qn, x3qw]   = msh_set_quad_nodes (x3_zeta, rule(3));
    x3msh = msh_cartesian (x3_zeta, x3qn, x3qw, x3geo);
    
    msh = struct('x1msh',x1msh,'x2msh',x2msh,'x3msh',x3msh,'tmsh',tmsh);

    % Define trial spaces
    x2sp_trial  = sp_bspline (x2_knots, trial_degree(2),  x2msh);
    x3sp_trial  = sp_bspline (x3_knots, trial_degree(3),  x3msh);
    
    space = struct('x1sp_trial',x1sp_trial,'x2sp_trial',x2sp_trial,...
                   'x3sp_trial',x3sp_trial,'tsp_trial',tsp_trial);

end


end

function [Aint, rhs, u, int_dofs, Afun, matrices] = schrodinger_st_problem_on_cartesian_domains(msh, space, problem_data, method_data)
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% univariate directions: 
rdim = numel(trial_degree);

% Assembly the matrices.
fprintf('Assembling univariate factors of A... \n\n')
Mt = op_u_v_tp               ( space.tsp_trial, space.tsp_trial, msh.tmsh); % mass in time
Lt = op_gradu_gradv_tp       ( space.tsp_trial, space.tsp_trial, msh.tmsh); % stif in time
Wt = op_gradu_v_tp           ( space.tsp_trial, space.tsp_trial, msh.tmsh); % Wt rows = trial, cols = test, so we have to transpose it. 
Wt = Wt'; 
% and symmetrize mass and stiffnes in time, and take dimension of internal
% dofs in time direction ---> all dofs are nt+1.
Mt = (Mt+Mt')/2; Lt = (Lt+Lt')/2; 
Mtint = Mt(2:end,2:end); Ltint = Lt(2:end,2:end); Wtint = Wt(2:end,2:end); 
nt = size(Mtint,1);
matrices.Mt = Mt(2:end,2:end);
matrices.Lt = Lt(2:end,2:end);
matrices.Wt = Wt(2:end,2:end);

switch rdim 
    case 2 % 1D space problems
        M1 = op_u_v_tp               ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % mass in space x1
        B1 = op_laplaceu_laplacev_tp ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % Blap in space x1
        L1 = op_laplaceu_v_tp        ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % Lap(u) vs w in space x1
        M1 = (M1+M1')/2; B1 = (B1+B1')/2; % mass and bi-Laplace symmetrized  
        M1int = M1(2:end-1,2:end-1); M1int = (M1int+M1int')/2; % internal mass symmetrized
        B1int = B1(2:end-1,2:end-1); B1int = (B1int+B1int')/2; % internal biLap symmetrized
        L1int = L1(2:end-1,2:end-1); L1int = (L1int+L1int')/2; % internal Lap symmetrized
        nx1 = size(M1int,1);
        matrices.M1 = M1int;
        matrices.L1 =-1*L1int;

    case 3 % 2D space problems
        M1 = op_u_v_tp               ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % mass in space x1
        B1 = op_laplaceu_laplacev_tp ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % Blap in space x1
        L1 = op_laplaceu_v_tp        ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % Lap(u) vs w in space x1
        M1 = (M1+M1')/2; B1 = (B1+B1')/2; % mass and bi-Laplace symmetrized 
        M1int = M1(2:end-1,2:end-1); M1int = (M1int+M1int')/2; % internal mass symmetrized
        B1int = B1(2:end-1,2:end-1); B1int = (B1int+B1int')/2; % internal biLap symmetrized
        L1int = L1(2:end-1,2:end-1); L1int = (L1int+L1int')/2; % internal Lap symmetrized
        nx1 = size(M1int,1);
        matrices.M1 = M1int;
        matrices.L1 =-1*L1int;

        M2 = op_u_v_tp               ( space.x2sp_trial, space.x2sp_trial, msh.x2msh); % mass in space x2
        B2 = op_laplaceu_laplacev_tp ( space.x2sp_trial, space.x2sp_trial, msh.x2msh); % Blap in space x2
        L2 = op_laplaceu_v_tp        ( space.x2sp_trial, space.x2sp_trial, msh.x2msh); % Lap(u) vs w in space x2
        M2 = (M2+M2')/2; B2 = (B2+B2')/2; % mass and bi-Laplace symmetrized 
        M2int = M2(2:end-1,2:end-1); M2int = (M2int+M2int')/2; % internal mass symmetrized
        B2int = B2(2:end-1,2:end-1); B2int = (B2int+B2int')/2; % internal biLap symmetrized
        L2int = L2(2:end-1,2:end-1); L2int = (L2int+L2int')/2; % internal Lap symmetrized
        nx2 = size(M2int,1);
        matrices.M2 = M2int;
        matrices.L2 =-1*L2int;

    case 4 % 3D space problems
        M1 = op_u_v_tp               ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % mass in space x1
        B1 = op_laplaceu_laplacev_tp ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % Blap in space x1
        L1 = op_laplaceu_v_tp        ( space.x1sp_trial, space.x1sp_trial, msh.x1msh); % Lap(u) vs w in space x1
        M1 = (M1+M1')/2; B1 = (B1+B1')/2; % mass and bi-Laplace symmetrized 
        M1int = M1(2:end-1,2:end-1); M1int = (M1int+M1int')/2; % internal mass symmetrized
        B1int = B1(2:end-1,2:end-1); B1int = (B1int+B1int')/2; % internal biLap symmetrized
        L1int = L1(2:end-1,2:end-1); L1int = (L1int+L1int')/2; % internal Lap symmetrized
        nx1 = size(M1int,1);
        matrices.M1 = M1int;
        matrices.L1 =-1*L1int;

        M2 = op_u_v_tp               ( space.x2sp_trial, space.x2sp_trial, msh.x2msh); % mass in space x2
        B2 = op_laplaceu_laplacev_tp ( space.x2sp_trial, space.x2sp_trial, msh.x2msh); % Blap in space x2
        L2 = op_laplaceu_v_tp        ( space.x2sp_trial, space.x2sp_trial, msh.x2msh); % Lap(u) vs w in space x2
        M2 = (M2+M2')/2; B2 = (B2+B2')/2; % mass and bi-Laplace symmetrized 
        M2int = M2(2:end-1,2:end-1); M2int = (M2int+M2int')/2; % internal mass symmetrized
        B2int = B2(2:end-1,2:end-1); B2int = (B2int+B2int')/2; % internal biLap symmetrized
        L2int = L2(2:end-1,2:end-1); L2int = (L2int+L2int')/2; % internal Lap symmetrized
        nx2 = size(M2int,1);
        matrices.M2 = M2int;
        matrices.L2 =-1*L2int;

        M3 = op_u_v_tp               ( space.x3sp_trial, space.x3sp_trial, msh.x3msh); % mass in space x3
        B3 = op_laplaceu_laplacev_tp ( space.x3sp_trial, space.x3sp_trial, msh.x3msh); % Blap in space x3
        L3 = op_laplaceu_v_tp        ( space.x3sp_trial, space.x3sp_trial, msh.x3msh); % Lap(u) vs w in space x3
        M3 = (M3+M3')/2; B3 = (B3+B3')/2; % mass and bi-Laplace symmetrized 
        M3int = M3(2:end-1,2:end-1); M3int = (M3int+M3int')/2; % internal mass symmetrized
        B3int = B3(2:end-1,2:end-1); B3int = (B3int+B3int')/2; % internal biLap symmetrized
        L3int = L3(2:end-1,2:end-1); L3int = (L3int+L3int')/2; % internal Lap symmetrized
        nx3 = size(M3int,1);
        matrices.M3 = M3int;
        matrices.L3 =-1*L3int;
end

% Assembly the rhs.
fprintf('Assembling the rhs... \n\n')
switch rdim 
    case 2
        % project boundary conditions and initial data:
        dof = reshape(1:(nx1+2)*(nt+1),nx1+2,nt+1);
        u   = zeros(size(dof));
        app1 = u; app2 = u; app3 = u;
        ux1 = op_f_v_tp(space.x1sp_trial, msh.x1msh, uex1);
        ut  = op_f_v_tp(space.tsp_trial,  msh.tmsh,  uext);

        UF0 = reshape(a0*ux1*uext(0),(nx1+2),1);
        UF1 = reshape(a0*ut*uex1(0), 1,(nt +1));
        UF2 = reshape(a0*ut*uex1(1), 1,(nt +1));

        app1(dof(:,1)) = tmprod(UF0,{inv(M1)}, 1); % projecting the initial data
        app2(dof(1,:)) = tmprod(UF1,{inv(Mt)}, 2); % projecting the Dirichlet data on F1
        app3(dof(end,:)) = tmprod(UF2,{inv(Mt)}, 2); % projecting the Dirichlet data on F1
        u = (app1+app2+app3)./((app1 ~= 0)+(app2 ~= 0)+(app3 ~= 0)+((app1+app2+app3)==0));
        u = u(:);

        % find space-time internal degrees of freedom
        drchlt_dof = find(u);
        int_dofs = setdiff (1:numel(dof), drchlt_dof); % internal d.o.f.

        % buil the right hand side by tensor product construction
        F = zeros(size(dof)); 

        for index = 1:rdim+1
            FM1 = op_f_v_tp  (space.x1sp_trial,msh.x1msh, fx1{index});
            FMt = op_f_v_tp  (space.tsp_trial ,msh.tmsh , ft{index});
            FL1 = op_f_lapv_tp(space.x1sp_trial,msh.x1msh,fx1{index});
            FKt = op_f_gradv_tp(space.tsp_trial,msh.tmsh, ft{index});
            F = F(:) -gmm*kron(FKt,FM1) -eta*kron(FMt,FL1);
        end
        
        if isequal(solver,'MB') ||  isequal(solver,'CG') || isequal(preconditioner ,'ichol') || isequal(preconditioner ,'ilu')
            A = kron(Lt,M1) + eta^2*kron(Mt,B1) -gmm*eta*kron(Wt,L1') + gmm*eta*kron(Wt',L1);
            Aint = A(int_dofs,int_dofs); 
            Afun = @(x) reshape(tmprod(reshape(x,nx1,nt), {M1int Ltint}, 1:2) +... 
                         eta^2* tmprod(reshape(x,nx1,nt), {B1int Mtint}, 1:2) +...
                       gmm*eta* tmprod(reshape(x,nx1,nt), {L1int Wtint'}, 1:2)-...
                       gmm*eta* tmprod(reshape(x,nx1,nt), {L1int' Wtint}, 1:2) ...
                               ,nx1*nt,1);
    
            [~,~,u_drchlt] = find(u); % find Dirichlet data.

            fprintf('Lifting... \n\n')
            rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt; % lifting.
        else
            Aint = 0;

            Afun = @(x) reshape(tmprod(reshape(x,nx1,nt), {M1int Ltint}, 1:2) +... 
                         eta^2* tmprod(reshape(x,nx1,nt), {B1int Mtint}, 1:2) +...
                       gmm*eta* tmprod(reshape(x,nx1,nt), {L1int Wtint'}, 1:2)-...
                       gmm*eta* tmprod(reshape(x,nx1,nt), {L1int' Wtint}, 1:2) ...
                               ,nx1*nt,1);
    
            Ad2i = @(x) reshape(tmprod(reshape(x,nx1+2,nt+1), {M1(2:end-1,:) Lt(2:end,:)}, 1:2) +... 
                         eta^2* tmprod(reshape(x,nx1+2,nt+1), {B1(2:end-1,:) Mt(2:end,:)}, 1:2) +...
                       gmm*eta* tmprod(reshape(x,nx1+2,nt+1), {L1(2:end-1,:) Wt(:,2:end)'}, 1:2)-...
                       gmm*eta* tmprod(reshape(x,nx1+2,nt+1), {L1(:,2:end-1)' Wt(2:end,:)}, 1:2) ...
                                ,nx1*nt,1);
            
            % Dirichlet conditions.
            fprintf('Lifting... \n\n')
            rhs  = F(int_dofs) - Ad2i(u);
        end                
    case 3
        dof = reshape(1:(nx1+2)*(nx2+2)*(nt+1),nx1+2,nx2+2,nt+1);
        u   = zeros(size(dof));
        app1 = u; app2 = u; app3 = u; app4 = u; app5 = u;

        ux1 = op_f_v_tp(space.x1sp_trial, msh.x1msh, uex1);
        ux2 = op_f_v_tp(space.x2sp_trial, msh.x2msh, uex2);
        ut  = op_f_v_tp(space.tsp_trial,  msh.tmsh,  uext);

        UF0 = reshape(a0*kron(ux2,ux1)*uext(0),(nx1+2),(nx2+2),1);
        UF1 = reshape(a0*kron(ut, ux2)*uex1(0),1,(nx2+2),(nt +1));
        UF2 = reshape(a0*kron(ut, ux2)*uex1(1),1,(nx2+2),(nt +1));
        UF3 = reshape(a0*kron(ut, ux1)*uex2(0),(nx1+2),1,(nt +1));
        UF4 = reshape(a0*kron(ut, ux1)*uex2(1),(nx1+2),1,(nt +1));

        app1(dof(:,:,1))   = tmprod(UF0,{inv(M1) inv(M2)}, 1:2); % projecting the initial data
        app2(dof(1,:,:))   = tmprod(UF1,{inv(M2) inv(Mt)}, 2:3); % projecting the Dirichlet data on F1
        app3(dof(end,:,:)) = tmprod(UF2,{inv(M2) inv(Mt)}, 2:3); % projecting the Dirichlet data on F2
        app4(dof(:,1,:))   = tmprod(UF3,{inv(M1) inv(Mt)}, [1 3]); % projecting the Dirichlet data on F3
        app5(dof(:,end,:)) = tmprod(UF4,{inv(M1) inv(Mt)}, [1 3]); % projecting the Dirichlet data on F4

        u = (app1+app2+app3+app4+app5)./((app1 ~= 0)+(app2 ~= 0)+(app3 ~= 0)+(app4 ~= 0)+(app5 ~= 0)+((app1+app2+app3+app4+app5)==0));
        u = u(:);

        % find space-time internal degrees of freedom
        drchlt_dof = find(u);
        int_dofs = setdiff (1:numel(dof), drchlt_dof); % internal d.o.f.

        % buil the right hand side by tensor product construction
        F = zeros(size(dof)); 
        
        for index = 1:rdim+1
        FM1 = op_f_v_tp  (space.x1sp_trial,msh.x1msh,fx1{index});
        FM2 = op_f_v_tp  (space.x2sp_trial,msh.x2msh,fx2{index});
        FMt = op_f_v_tp  (space.tsp_trial ,msh.tmsh ,ft{index});
        FL1 = op_f_lapv_tp  (space.x1sp_trial,msh.x1msh,fx1{index});
        FL2 = op_f_lapv_tp  (space.x2sp_trial,msh.x2msh,fx2{index});
        FKt = op_f_gradv_tp (space.tsp_trial ,msh.tmsh ,ft{index});

        F = F(:) -gmm* kron(FKt,kron(FM2,FM1)) -eta*(kron(FMt, kron(FL2,FM1) + kron(FM2,FL1)));

        end

        if isequal(solver,'MB') ||  isequal(solver,'CG') || isequal(preconditioner ,'ichol') || isequal(preconditioner ,'ilu')
            A = kron(Lt,kron(M2,M1)) + eta^2*kron(Mt,kron(B2,M1)+kron(M2,B1)+2*kron(L2,L1)) ...
                -gmm*eta*kron(Wt ,kron(M2,L1')+kron(L2',M1)) ...
                +gmm*eta*kron(Wt',kron(M2,L1 )+kron(L2 ,M1));

            Aint = A(int_dofs,int_dofs); 
            Afun = @(x) reshape(tmprod(reshape(x,nx1,nx2,nt), {M1int M2int Ltint}, 1:3) +... 
                         eta^2*(tmprod(reshape(x,nx1,nx2,nt), {B1int M2int Mtint}, 1:3) +...
                                tmprod(reshape(x,nx1,nx2,nt), {M1int B2int Mtint}, 1:3) +...
                             2* tmprod(reshape(x,nx1,nx2,nt), {L1int L2int Mtint}, 1:3))+...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nt), {L1int M2int Wtint'}, 1:3)+...
                                tmprod(reshape(x,nx1,nx2,nt), {M1int L2int Wtint'}, 1:3))-...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nt), {L1int' M2int Wtint}, 1:3)+...
                                tmprod(reshape(x,nx1,nx2,nt), {M1int L2int' Wtint}, 1:3))...
                               ,nx1*nx2*nt,1);
    
            [~,~,u_drchlt] = find(u); % find Dirichlet data.

            fprintf('Lifting... \n\n')
            rhs = F(int_dofs) - A(int_dofs, drchlt_dof)*u_drchlt; % lifting.
        else
            Aint = 0;
            Afun = @(x) reshape(tmprod(reshape(x,nx1,nx2,nt), {M1int M2int Ltint}, 1:3) +... 
                         eta^2*(tmprod(reshape(x,nx1,nx2,nt), {B1int M2int Mtint}, 1:3) +...
                                tmprod(reshape(x,nx1,nx2,nt), {M1int B2int Mtint}, 1:3) +...
                             2* tmprod(reshape(x,nx1,nx2,nt), {L1int L2int Mtint}, 1:3))+...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nt), {L1int M2int Wtint'}, 1:3)+...
                                tmprod(reshape(x,nx1,nx2,nt), {M1int L2int Wtint'}, 1:3))-...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nt), {L1int' M2int Wtint}, 1:3)+...
                                tmprod(reshape(x,nx1,nx2,nt), {M1int L2int' Wtint}, 1:3))...
                               ,nx1*nx2*nt,1);
    
            Ad2i = @(x) reshape(tmprod(reshape(x,nx1+2,nx2+2,nt+1), {M1(2:end-1,:) M2(2:end-1,:) Lt(2:end,:)}, 1:3) +... 
                         eta^2*(tmprod(reshape(x,nx1+2,nx2+2,nt+1), {B1(2:end-1,:) M2(2:end-1,:) Mt(2:end,:)}, 1:3) +...
                                tmprod(reshape(x,nx1+2,nx2+2,nt+1), {M1(2:end-1,:) B2(2:end-1,:) Mt(2:end,:)}, 1:3) +...
                             2* tmprod(reshape(x,nx1+2,nx2+2,nt+1), {L1(2:end-1,:) L2(2:end-1,:) Mt(2:end,:)}, 1:3))+...
                       gmm*eta*(tmprod(reshape(x,nx1+2,nx2+2,nt+1), {L1(2:end-1,:) M2(2:end-1,:) Wt(:,2:end)'}, 1:3)+...
                                tmprod(reshape(x,nx1+2,nx2+2,nt+1), {M1(2:end-1,:) L2(2:end-1,:) Wt(:,2:end)'}, 1:3))-...
                       gmm*eta*(tmprod(reshape(x,nx1+2,nx2+2,nt+1), {L1(:,2:end-1)' M2(2:end-1,:) Wt(2:end,:)}, 1:3)+...
                                tmprod(reshape(x,nx1+2,nx2+2,nt+1), {M1(2:end-1,:) L2(:,2:end-1)' Wt(2:end,:)}, 1:3))...
                                ,nx1*nx2*nt,1);
            % Dirichlet conditions.
            fprintf('Lifting... \n\n')
            rhs  = F(int_dofs) - Ad2i(u);
        end
    case 4
        dof = reshape(1:(nx1+2)*(nx2+2)*(nx3+2)*(nt+1),nx1+2,nx2+2,nx3+2,nt+1);
        u   = zeros(size(dof));
        app1 = u; app2 = u; app3 = u; app4 = u; app5 = u; app6 = u; app7 = u;

        ux1 = op_f_v_tp(space.x1sp_trial, msh.x1msh, uex1);
        ux2 = op_f_v_tp(space.x2sp_trial, msh.x2msh, uex2);
        ux3 = op_f_v_tp(space.x3sp_trial, msh.x3msh, uex3);
        ut  = op_f_v_tp(space.tsp_trial,  msh.tmsh,  uext);

        UF0 = reshape(a0*kron(ux3,kron(ux2,ux1))*uext(0),(nx1+2),(nx2+2),(nx3+2),1);
        UF1 = reshape(a0*kron(ut, kron(ux3,ux2))*uex1(0),1,(nx2+2),(nx3+2),(nt +1));
        UF2 = reshape(a0*kron(ut, kron(ux3,ux2))*uex1(1),1,(nx2+2),(nx3+2),(nt +1));
        UF3 = reshape(a0*kron(ut, kron(ux3,ux1))*uex2(0),(nx1+2),1,(nx3+2),(nt +1));
        UF4 = reshape(a0*kron(ut, kron(ux3,ux1))*uex2(1),(nx1+2),1,(nx3+2),(nt +1));
        UF5 = reshape(a0*kron(ut, kron(ux2,ux1))*uex3(0),(nx1+2),(nx2+2),1,(nt +1));
        UF6 = reshape(a0*kron(ut, kron(ux2,ux1))*uex3(1),(nx1+2),(nx2+2),1,(nt +1));

        app1(dof(:,:,:,1))   = tmprod(UF0,{inv(M1) inv(M2) inv(M3)}, 1:3); % projecting the initial data
        app2(dof(1,:,:,:))   = tmprod(UF1,{inv(M2) inv(M3) inv(Mt)}, 2:4); % projecting the Dirichlet data on F1
        app3(dof(end,:,:,:)) = tmprod(UF2,{inv(M2) inv(M3) inv(Mt)}, 2:4); % projecting the Dirichlet data on F2
        app4(dof(:,1,:,:))   = tmprod(UF3,{inv(M1) inv(M3) inv(Mt)}, [1 3 4]); % projecting the Dirichlet data on F3
        app5(dof(:,end,:,:)) = tmprod(UF4,{inv(M1) inv(M3) inv(Mt)}, [1 3 4]); % projecting the Dirichlet data on F4
        app6(dof(:,:,1,:))   = tmprod(UF5,{inv(M1) inv(M2) inv(Mt)}, [1 2 4]); % projecting the Dirichlet data on F3
        app7(dof(:,:,end,:)) = tmprod(UF6,{inv(M1) inv(M2) inv(Mt)}, [1 2 4]); % projecting the Dirichlet data on F4

        u = (app1+app2+app3+app4+app5+app6+app7)./...
            ((app1 ~= 0)+(app2 ~= 0)+(app3 ~= 0)+(app4 ~= 0)+(app5 ~= 0)+(app6 ~= 0)+(app7 ~= 0)+...
            ((app1+app2+app3+app4+app5+app6+app7)==0));
        u = u(:);

        % find space-time internal degrees of freedom
        drchlt_dof = find(u);
        int_dofs = setdiff (1:numel(dof), drchlt_dof); % internal d.o.f.

        % buil the right hand side by tensor product construction
        F = zeros(size(dof)); 
                
        for index = 1:rdim+1
        FM1 = op_f_v_tp  (space.x1sp_trial,msh.x1msh,fx1{index});
        FM2 = op_f_v_tp  (space.x2sp_trial,msh.x2msh,fx2{index});
        FM3 = op_f_v_tp  (space.x3sp_trial,msh.x3msh,fx3{index});
        FMt = op_f_v_tp  (space.tsp_trial ,msh.tmsh , ft{index});
        FL1 = op_f_lapv_tp  (space.x1sp_trial,msh.x1msh,fx1{index});
        FL2 = op_f_lapv_tp  (space.x2sp_trial,msh.x2msh,fx2{index});
        FL3 = op_f_lapv_tp  (space.x3sp_trial,msh.x3msh,fx3{index});
        FKt = op_f_gradv_tp (space.tsp_trial ,msh.tmsh , ft{index});

        F = F(:) -gmm* kron(FKt,  kron(kron(FM3,FM2),FM1))...
                 -eta*(kron(FMt, (kron(kron(FL3,FM2),FM1) ...
                                + kron(kron(FM3,FL2),FM1) ...
                                + kron(kron(FM3,FM2),FL1))));
        end

        if isequal(solver,'MB') ||  isequal(solver,'CG') || isequal(preconditioner ,'ichol') || isequal(preconditioner ,'ilu')
            A =       kron(Lt ,kron(M3,kron(M2,M1))) +...
                eta^2*kron(Mt ,kron(B3,kron(M2,M1)) + kron(M3,kron(B2,M1)) + kron(M3,kron(M2,B1)) ...
                          +2* (kron(M3,kron(L2,L1)) + kron(L3,kron(M2,L1)) + kron(L3,kron(L2,M1)))) ...
             -gmm*eta*kron(Wt ,kron(M3,kron(M2,L1'))+ kron(M3,kron(L2',M1))+ kron(L3',kron(M2,M1))) ...
             +gmm*eta*kron(Wt',kron(M3,kron(M2,L1 ))+ kron(M3,kron(L2 ,M1))+ kron(L3 ,kron(M2,M1)));

            Aint = A(int_dofs,int_dofs); 

            Afun = @(x) reshape(tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int M3int Ltint}, 1:4) +... 
                         eta^2*(tmprod(reshape(x,nx1,nx2,nx3,nt), {B1int M2int M3int Mtint}, 1:4) +...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int B2int M3int Mtint}, 1:4) +...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int B3int Mtint}, 1:4) +...
                             2* tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int L2int L3int Mtint}, 1:4) +...
                             2* tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int M2int L3int Mtint}, 1:4) +...
                             2* tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int L2int M3int Mtint}, 1:4))+...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int M2int M3int Wtint'}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int L2int M3int Wtint'}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int L3int Wtint'}, 1:4))-...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int' M2int M3int Wtint}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int L2int' M3int Wtint}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int L3int' Wtint}, 1:4))...
                                ,nx1*nx2*nx3*nt,1);
    
            [~,~,u_drchlt] = find(u); % find Dirichlet data.

            fprintf('Lifting... \n\n')
            rhs = F(int_dofs) - A(int_dofs, drchlt_dof)*u_drchlt; % lifting.
        else
            Aint = 0;
            Afun = @(x) reshape(tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int M3int Ltint}, 1:4) +... 
                         eta^2*(tmprod(reshape(x,nx1,nx2,nx3,nt), {B1int M2int M3int Mtint}, 1:4) +...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int B2int M3int Mtint}, 1:4) +...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int B3int Mtint}, 1:4) +...
                             2* tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int L2int L3int Mtint}, 1:4) +...
                             2* tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int M2int L3int Mtint}, 1:4) +...
                             2* tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int L2int M3int Mtint}, 1:4))+...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int M2int M3int Wtint'}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int L2int M3int Wtint'}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int L3int Wtint'}, 1:4))-...
                       gmm*eta*(tmprod(reshape(x,nx1,nx2,nx3,nt), {L1int' M2int M3int Wtint}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int L2int' M3int Wtint}, 1:4)+...
                                tmprod(reshape(x,nx1,nx2,nx3,nt), {M1int M2int L3int' Wtint}, 1:4))...
                                ,nx1*nx2*nx3*nt,1);
    
            Ad2i = @(x) reshape(tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) M2(2:end-1,:) M3(2:end-1,:) Lt(2:end,:)}, 1:4) +... 
                         eta^2*(tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {B1(2:end-1,:) M2(2:end-1,:) M3(2:end-1,:) Mt(2:end,:)}, 1:4) +...
                                tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) B2(2:end-1,:) M3(2:end-1,:) Mt(2:end,:)}, 1:4) +...
                                tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) M2(2:end-1,:) B3(2:end-1,:) Mt(2:end,:)}, 1:4) +...
                             2* tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) L2(2:end-1,:) L3(2:end-1,:) Mt(2:end,:)}, 1:4) +...
                             2* tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {L1(2:end-1,:) M2(2:end-1,:) L3(2:end-1,:) Mt(2:end,:)}, 1:4) +...
                             2* tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {L1(2:end-1,:) L2(2:end-1,:) M3(2:end-1,:) Mt(2:end,:)}, 1:4))+...
                       gmm*eta*(tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {L1(2:end-1,:) M2(2:end-1,:) M3(2:end-1,:) Wt(:,2:end)'}, 1:4)+...
                                tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) L2(2:end-1,:) M3(2:end-1,:) Wt(:,2:end)'}, 1:4)+...
                                tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) M2(2:end-1,:) L3(2:end-1,:) Wt(:,2:end)'}, 1:4))-...
                       gmm*eta*(tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {L1(:,2:end-1)' M2(2:end-1,:) M3(2:end-1,:) Wt(2:end,:)}, 1:4)+...
                                tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) L2(:,2:end-1)' M3(2:end-1,:) Wt(2:end,:)}, 1:4)+...
                                tmprod(reshape(x,nx1+2,nx2+2,nx3+2,nt+1), {M1(2:end-1,:) M2(2:end-1,:) L3(:,2:end-1)' Wt(2:end,:)}, 1:4))...
                                ,nx1*nx2*nx3*nt,1);
            
            % Dirichlet conditions.
            fprintf('Lifting... \n\n')
            rhs  = F(int_dofs) - Ad2i(u);
        end

end


end

function [Aint, rhs, u, int_dofs, Afun] = schrodinger_st_problem(msh, space, problem_data, method_data)
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Apply Dirichlet bpoundary conditions
u = zeros (space.xtsp_trial.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space.xtsp_trial, msh.xtmsh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt; % boundary d.o.f.
% find space-time internal degrees of freedom
int_dofs = setdiff (1:space.xtsp_trial.ndof, drchlt_dofs); % internal d.o.f.
% find only space-boundary degrees of freedom
[~, x_drchlt_dofs] = sp_drchlt_l2_proj (space.xsp_trial, msh.xmsh, hx, x_drchlt_sides); % DEFINISCI HX in script
% find only space-internal degrees of freedom
x_int_dofs = setdiff (1:space.xsp_trial.ndof, x_drchlt_dofs); 
% number of elements of internal and Dirichlet degrees of freedom
intnx = numel(x_int_dofs); drchnx = numel(x_drchlt_dofs); 
intnt = space.tsp_trial.ndof-1;

% Assembly the matrices.
fprintf('Assembling A... \n\n')
Mt = op_u_v_tp               ( space.tsp_trial, space.tsp_test, msh.tmsh); % mass in time
Lt = op_gradu_gradv_tp       ( space.tsp_trial, space.tsp_test, msh.tmsh); % stif in time
Wt = op_gradu_v_tp           ( space.tsp_trial, space.tsp_test, msh.tmsh); % Wt rows = trial, cols = test, so we have to transpose it. 
Wt = Wt'; 

% when we have a geometry
Ms = op_u_v_tp               ( space.xsp_trial, space.xsp_test, msh.xmsh); % mass in space
Bs = op_laplaceu_laplacev_tp ( space.xsp_trial, space.xsp_test, msh.xmsh); % Blap in space
Ls = op_laplaceu_v_tp        ( space.xsp_trial, space.xsp_test, msh.xmsh); % Lap(u) vs w 
% Ls = op_gradu_gradv_tp     ( space.xsp_trial, space.xsp_test, msh.xmsh); % stiff in space 
% Ls = -Ls ;
% symmetrize 
Mt = (Mt+Mt')/2; Ms = (Ms+Ms')/2;
Lt = (Lt+Lt')/2; Bs = (Bs+Bs')/2;


fprintf('Building rhs... \n\n')
F = op_f_sv_tp (space.xtsp_test, msh.xtmsh, f); % here is the projection of f in S(V)
% F = op_f_v_tp (space.xtsp_test, msh.xtmsh, f); % here is the projection of f in S(V)

if isequal(solver,'MB') 
  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
  Aint = A(int_dofs,int_dofs); 
  Afun = @(x) Aint*x;

  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt;
elseif isequal(preconditioner ,'ichol') || isequal(preconditioner ,'ilu') || isequal(solver,'CG')
  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
  Aint = A(int_dofs,int_dofs); 
  Afun = @(x) reshape(Ms(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Lt(2:end,2:end)') + ...
                      Bs(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Mt(2:end,2:end)') + ...
                      Ls(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(1i* Wt(2:end,2:end))  - ...
                     (Ls(x_int_dofs,x_int_dofs)')*reshape(x,intnx,intnt) *(1i*(Wt(2:end,2:end)')),  ...
                      intnx*intnt, 1);
  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt;
else
  Aint = 0;
  A = kron(Lt,Ms) + kron(Mt,Bs) -1i*kron(Wt,Ls') + 1i*kron(Wt',Ls);
  Afun = @(x) reshape(Ms(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Lt(2:end,2:end)') + ...
                      Bs(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(    Mt(2:end,2:end)') + ...
                      Ls(x_int_dofs,x_int_dofs)  *reshape(x,intnx,intnt) *(1i* Wt(2:end,2:end))  - ...
                     (Ls(x_int_dofs,x_int_dofs)')*reshape(x,intnx,intnt) *(1i*(Wt(2:end,2:end)')),  ...
                      intnx*intnt, 1);

  Adrch2int = @(x) reshape(Ms(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(    Lt(2:end,2:end)') + ...
                           Bs(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(    Mt(2:end,2:end)') + ...
                           Ls(x_int_dofs,x_drchlt_dofs)  *reshape(x,drchnx,intnt) *(1i* Wt(2:end,2:end) ) - ...
                          (Ls(x_drchlt_dofs,x_int_dofs)')*reshape(x,drchnx,intnt) *(1i*(Wt(2:end,2:end)')),  ...
                           intnx*intnt, 1);

  fprintf('Lifting... \n\n')
  rhs  = F(int_dofs) - Adrch2int(u_drchlt(space.xsp_trial.ndof+1:end));
  fprintf('Lifting... \n\n')
  rhs = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt;
end

% we should also give in output the univariate matrices.. 
end

%% NOTA BENE: 
%  abbiamo il seguente sistema da invertire Au = f con 
%
%   A = Ms x Kt + \nu^2 Bs x Mt -i\nu Ls' x Wt + i\nu Ls x Wt' 
%  
%  e lo risolviamo iterativamente con PCG dove il precondizionatore é 
%
%   P = Ms x Kt + \nu^2 KsMs^(-1)Ks x Mt + i\nu Ks x Wt - i\nu Ks x Wt' 
%
%  questo perchè Ks è la matrice di stiffness ovvero -Ls (salvo errori)
%
%  CONTROLLA IL CASO 2D E IL CASO GENERICO, QUESTIONE DI SEGNI....

function prec = LU_SETUP_x_SCHRODINGER_1D(Lt,Mt,Wt,Asx,Msx)
nsx = size(Asx,1); nt = size(Lt,1);
%Fast diagonalization in univariate space direction (SETUP) 
[Us, Ds] = eig(full(Asx),full(Msx),'vector'); 
B = kron(speye(numel(Ds)), Lt) + ...
    kron(spdiags(reshape(Ds.^2,numel(Ds),1), 0, numel(Ds),numel(Ds)),Mt) - ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt') + ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt);
dB = decomposition(B,'banded','CheckCondition',false);
prec = @(v) LU_SOLVER_x_SCHRODINGER_1D(dB,Us,v,nsx,nt);
end

function prec = LU_SETUP_x_SCHRODINGER_2D(Lt,Mt,Wt,Lsx,Msx,Lsy,Msy)
nsx = size(Lsx,1); nsy = size(Lsy,1); nt = size(Lt,1);
%Fast diagonalization in univariate space direction (SETUP) 
[Ux, Dx] = eig(full(Lsx),full(Msx),'vector'); 
[Uy, Dy] = eig(full(Lsy),full(Msy),'vector');
% assembling and decomposing:
Ds = reshape(Dy + Dx',[],1);
B = kron(speye(numel(Ds)), Lt) + ...
    kron(spdiags(reshape(Ds.^2,numel(Ds),1), 0, numel(Ds),numel(Ds)),Mt) - ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt') + ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt);
dB = decomposition(B,'banded','CheckCondition',false);
prec = @(v) LU_SOLVER_x_SCHRODINGER_2D(dB,Ux,Uy,v,nsx,nsy,nt);
end

function prec = LU_SETUP_x_SCHRODINGER_3D(Lt,Mt,Wt,Lsx,Msx,Lsy,Msy,Lsz,Msz)
nsx = size(Lsx,1); nsy = size(Lsy,1); nsz = size(Lsz,1); nt = size(Lt,1);
%Fast diagonalization in univariate space direction (SETUP) 
[Ux, Dx] = eig(full(Lsx),full(Msx),'vector'); 
[Uy, Dy] = eig(full(Lsy),full(Msy),'vector');
[Uz, Dz] = eig(full(Lsz),full(Msz),'vector');
% assembling and decomposing:
Ds = reshape(Dz + Dy' + reshape(Dx,1,1,numel(Dx)),[],1);
B = kron(speye(numel(Ds)), Lt) + ...
    kron(spdiags(reshape(Ds.^2,numel(Ds),1), 0, numel(Ds),numel(Ds)),Mt) - ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt') + ...
    kron(spdiags(reshape(Ds,numel(Ds),1), 0, numel(Ds),numel(Ds)), 1i*Wt);
dB = decomposition(B,'banded','CheckCondition',false);
prec = @(v) LU_SOLVER_x_SCHRODINGER_3D(dB,Ux,Uy,Uz,v,nsx,nsy,nsz,nt);
end


function u_lu = LU_SOLVER_x_SCHRODINGER_1D(lu_B,Ux,rhs,nsx,nt)
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nt), {Ux'}, 1);  
tilde_rhs = reshape (permute(tilde_rhs,[2 1]), nt*nsx,1);

% STEP 2
tilde_u = lu_B\tilde_rhs;
% STEP 3
new_tilde_u = reshape(permute(reshape(tilde_u,nt,nsx),[2 1]), nsx, nt);
u_lu = tmprod(new_tilde_u, {Ux }, 1);  
u_lu = u_lu(:);
end

function u_lu = LU_SOLVER_x_SCHRODINGER_2D(lu_B,Ux,Uy,rhs,nsx,nsy,nt)
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nsy,nt), {Ux' Uy'}, 1:2);  
tilde_rhs = reshape (permute(tilde_rhs,[3 1 2]), nt*nsx*nsy,1);
% STEP 2
tilde_u = lu_B\tilde_rhs;
% STEP 3
new_tilde_u = reshape(permute(reshape(tilde_u,nt,nsx,nsy),[2 3 1]), nsx, nsy, nt);
u_lu = tmprod(new_tilde_u, {Ux Uy}, 1:2);  
u_lu = u_lu(:);
end

function u_lu = LU_SOLVER_x_SCHRODINGER_3D(lu_B,Ux,Uy,Uz,rhs,nsx,nsy,nsz,nt)
% STEP 1
tilde_rhs = tmprod(reshape(rhs,nsx,nsy,nsz,nt), {Ux' Uy' Uz'}, 1:3);  
tilde_rhs = reshape (permute(tilde_rhs,[4 1 2 3]), nt*nsx*nsy*nsz,1);
% STEP 2
tilde_u = lu_B\tilde_rhs;
% STEP 3
new_tilde_u = reshape(permute(reshape(tilde_u,nt,nsx,nsy,nsz),[2 3 4 1]), nsx, nsy, nsz, nt);
u_lu = tmprod(new_tilde_u, {Ux Uy Uz}, 1:3);  
u_lu = u_lu(:);
end


