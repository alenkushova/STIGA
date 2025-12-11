function  errh1s = st_seminorm_h1s_error (space,spaceT,msh,mshT,u,coeff)
u=reshape(u,space.ndof,spaceT.ndof);
eut = []; val_gradu = [];
for sdof = 1:space.ndof
    [eu, ~] = sp_eval_msh (u(sdof,:), spaceT, mshT);
    eut = cat(1,eut,eu(:)');
end
for tpnt = 1: size(eut,2)
    [egradus, ~] = sp_eval_msh(eut(:,tpnt), space, msh, 'gradient');
    val_gradu = cat(4,val_gradu,egradus);
end
val_gradu = reshape(val_gradu, space.ncomp, msh.ndim, msh.nqn, msh.nel, mshT.nqn, mshT.nel);
coeff = reshape(coeff, space.ncomp, msh.ndim, msh.nqn, msh.nel, mshT.nqn, mshT.nel);
ws = msh.quad_weights .* msh.jacdet;
errl2_space = squeeze( sum ( ws.*reshape(sum((val_gradu-coeff).^2,[1 2]),...
                                  msh.nqn, msh.nel, mshT.nqn, mshT.nel),...
                              [1 2]));
wt = mshT.quad_weights .* mshT.jacdet;
errh1s = sqrt(sum(wt.*errl2_space,[1 2]));
end
