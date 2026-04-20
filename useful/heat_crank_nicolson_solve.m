function [geo, msh, space, u] = ...
                      heat_crank_nicolson_solve (problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Building the discrete spaces:
[geo, msh, space] = discretize(problem_data,method_data); 

% Assembling the linear system and the right hand side:
Ms = op_u_v_tp ( space.xsp_trial, space.xsp_test, msh.xmsh);               % mass in space
Ks = op_gradu_gradv_tp ( space.xsp_trial, space.xsp_test, msh.xmsh);       % stif in space
Ms = (Ms+Ms')/2; Ks = (Ks+Ks')/2; 

% Initialize the solution as null vector
u = zeros (space.xsp_trial.ndof, n_time_intervals+1);
rhs = zeros (space.xsp_test.ndof, n_time_intervals+1);
% We have uniform partitions of size dt
dt = T/n_time_intervals;

% Dirichlet dofs and internal dofs (just the dofs we want)
[~, drchlt_dofs] = sp_drchlt_l2_proj (space.xsp_trial, msh.xmsh,@(x, iside) dfun(x,0,iside), drchlt_sides);
int_dofs = setdiff (1:space.xsp_trial.ndof, drchlt_dofs); 

% Initial data by L2-projection
u(:,1) = Ms \ op_f_v_tp (space.xsp_trial, msh.xmsh, @(x) uex(x,0));
rhs(:,1) = op_f_v_tp(space.xsp_trial, msh.xmsh, @(x) f(x,0));

% Loop for the Crank-Nicolson method in time
for n_time = 1:n_time_intervals
 [u_drchlt, ~] = sp_drchlt_l2_proj (space.xsp_trial, msh.xmsh, @(x, iside) dfun(x,0+dt*n_time,iside), drchlt_sides);
 u(drchlt_dofs, n_time+1) = u_drchlt; 
 rhs(:,n_time+1) = op_f_v_tp(space.xsp_trial, msh.xmsh, @(x) f(x,0+dt*n_time));
 u(int_dofs,n_time+1) = (Ms(int_dofs,int_dofs) + dt/2*Ks(int_dofs,int_dofs))\... 
                   ((Ms(int_dofs,:) - dt/2*Ks(int_dofs,:))*u(:,n_time) + ...
                     dt/2*(rhs(int_dofs,n_time+1)+rhs(int_dofs,n_time)) - ...
                    (Ms(int_dofs,drchlt_dofs) + dt/2*Ks(int_dofs,drchlt_dofs))*u(drchlt_dofs,n_time+1));
end

end