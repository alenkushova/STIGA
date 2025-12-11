function  errh1t = st_seminorm_h1t_error (space, spaceT, msh, mshT, u, coeff)
u=reshape(u,space.ndof,spaceT.ndof);
edtut = []; valdtu = [];
for sdof = 1:space.ndof
    [edtu, ~] = sp_eval_msh (u(sdof,:), spaceT, mshT, 'gradient');
    edtut = cat(1,edtut,edtu(:)');
end
for tpnt = 1: size(edtut,2)
    [eus, ~] = sp_eval_msh(edtut(:,tpnt), space, msh);
    valdtu = cat(3,valdtu,eus);
end
valdtu = reshape(valdtu, space.ncomp, msh.nqn, msh.nel, mshT.nqn, mshT.nel);
coeff = reshape(coeff, space.ncomp, msh.nqn, msh.nel, mshT.nqn, mshT.nel);
ws = msh.quad_weights .* msh.jacdet;
errl2_space = squeeze( sum ( ws.*reshape(sum((valdtu-coeff).^2,1),...
                                  msh.nqn, msh.nel, mshT.nqn, mshT.nel),...
                              [1 2]));
wt = mshT.quad_weights .* mshT.jacdet;
errh1t = sqrt(sum(wt.*errl2_space,[1 2]));
end