function errl2 = st_l2_error_pressures_tp(spaceS,spaceT,mshS,mshT,p,pex)
  if (numel(p) ~= spaceS.ndof*spaceT.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  errl2 = 0;
  for iel = 1:mshT.nel_dir(1)
    msh_colT = msh_evaluate_col (mshT, iel);
    sp_colT  = sp_evaluate_col (spaceT, msh_colT, 'value', true, 'gradient', false);
    
    msh_col = msh_precompute(mshS);
    sp_col = sp_precompute(spaceS,msh_col);

    for idim = 1:mshS.rdim
      x{idim} = repmat(reshape(msh_col.geo_map(idim,:,:),msh_col.nqn,msh_col.nel),[1,1,msh_colT.nqn,msh_colT.nel]);
    end
    x{idim+1} = repmat(reshape(msh_colT.geo_map,1,1,msh_colT.nqn,msh_colT.nel),[msh_col.nqn,msh_col.nel,1,1]);

    ws = msh_col.quad_weights .* msh_col.jacdet; 

    coef = pex(x{:})- sum(ws.*pex(x{:}),[1,2]);

    errl2 = errl2 + (st_l2_error (sp_col,sp_colT,msh_col,msh_colT,p,coef)).^2;
  end
  
  errl2 = sqrt (errl2);

end
