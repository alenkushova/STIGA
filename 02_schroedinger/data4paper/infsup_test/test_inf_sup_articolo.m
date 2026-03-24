%========================
% THIS FUNCTION PRODUCES THE DATA IN THE DATABASE.
% you may change this two parameters to produce the infsup constant for the
% plots in Figure 1 of the paper. It takes obvious time to reproduce all
% the data in database. 

p = 2; % polinomial degree
disc = 'gal'; % discretization: plain Galerkin 'gal' or least squares 'ls'
%========================

eigval = zeros(1,5); infsup = zeros(1,5); h = zeros(1,5);
geometries = {}; meshes = {}; spaces = {}; eigvect = {}; lamMww = {};
v = {}; Av = {};

% 'filename' of the file where to save the results
filename = [disc '_degree_' num2str(p) '_infsup_trial_graf_test_l2_norms.mat'];
%load(filename)
index = 0;
for i = index+1 : index+128
  nel = i+4; j = i;
  switch disc
    case 'gal'
      [mu, ei, geo, msh, space, w]  = pg_infsuptest (p, nel);
    case 'ls'
      [mu, ei, geo, msh, space, w]  = ls_infsuptest (p, nel);
  end
  h (j) = 1/nel;   infsup (j) = mu;   eigval (j) = ei;
  eigvect    = cat(2, eigvect, {w});
  geometries = cat(2, geometries, {geo});
  meshes     = cat(2, meshes, {msh});
  spaces     = cat(2, spaces, {space});
% Saving at any loop iteration:
  save( filename , 'h', 'infsup', 'eigval', 'eigvect', 'geometries', 'meshes', 'spaces')
  fprintf(['Iteration ' num2str(i) '. Number of elements ' num2str(nel) '. Results saved in file: ' filename '\n\n'])
end
clear 

