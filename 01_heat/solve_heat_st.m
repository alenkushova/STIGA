% SOLVE_HEAT_ST
% 
% DOC HERE:
%
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