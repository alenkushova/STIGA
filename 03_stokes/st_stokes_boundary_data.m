function [vel_drchlt, drchlt_dofs, vel_iniz] = ... 
    st_stokes_boundary_data(spv, spt, xmsh, tmsh, h, drchlt_sides,isitweak)
  arguments
    spv
    spt
    xmsh
    tmsh
    h
    drchlt_sides
    isitweak = 'no'
  end
% proiezione dei dati iniziali.
dim = spv.ncomp;
switch dim
 case 2 % 2D in spazio
  vel0 = op_f_v_tp (spv, xmsh, @(x,y) h(x,y,0)); % valutazione di \int_\Omega u0 * v dx
 case 3 % 3D in spazio
  vel0 = op_f_v_tp (spv, xmsh, @(x,y,z) h(x,y,z,0)); % valutazione di \int_\Omega u0 * v dx
end
Mx = op_u_v_tp (spv, spv, xmsh); % L2 repr. matrix to project the data
fprintf('Projecting initial data... \n\n')
vel_iniz = Mx\vel0; % L2 projection of initial data

% DIRICHLET DATA IN SPACE-TIME
fprintf('Projecting Dirichlet boundary data... \n\n')
[vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_st (spv, spt, xmsh, tmsh, h, drchlt_sides,isitweak);

end