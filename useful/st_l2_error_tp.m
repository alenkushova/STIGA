function errl2 = st_l2_error_tp(spaceS,spaceT,mshS,mshT,p,pex)

  if (numel(p) ~= spaceS.ndof*spaceT.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  errl2 = 0;
  
  for iel = 1:mshS.nel_dir(1)
    msh_col = msh_evaluate_col (mshS, iel);
    sp_col  = sp_evaluate_col (spaceS, msh_col, 'value', true, 'gradient', false);
    
    msh_colT = msh_precompute(mshT);
    sp_colT = sp_precompute(spaceT,msh_colT);

    for idim = 1:mshS.rdim
      x{idim} = repmat(reshape(msh_col.geo_map(idim,:,:),msh_col.nqn,msh_col.nel),[1,1,msh_colT.nqn,msh_colT.nel]);
    end
    x{idim+1} = repmat(reshape(msh_colT.map(msh_colT.qn),1,1,msh_colT.nqn,msh_colT.nel),[msh_col.nqn,msh_col.nel,1,1]);

    errl2 = errl2 + (st_l2_error (sp_col,sp_colT,msh_col,msh_colT,p,pex(x{:}))).^2;
  end
  
  errl2 = sqrt (errl2);

end
