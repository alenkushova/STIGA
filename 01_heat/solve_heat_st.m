function [geo, msh, space, u, report] = solve_heat_st(problem_data,method_data)
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
[Aint, rhs, u, int_dofs, Afun] = heat_st_problem(msh, space, problem_data, method_data);

end