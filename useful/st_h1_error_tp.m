function [errl2, errh1s, errh1t] = st_h1_error_tp(spaceS,spaceT,mshS,mshT,v,vex,dtvex,gradvex)

  if (numel(v) ~= spaceS.ndof*spaceT.ndof)
    error ('Wrong size of the vector of degrees of freedom')
  end

  errl2 = 0;
  errh1t = 0;
  errh1s = 0;
  
  for iel = 1:mshS.nel_dir(1)
    msh_col = msh_evaluate_col (mshS, iel);
    sp_col  = sp_evaluate_col (spaceS, msh_col, 'value', true, 'gradient', true);
    
    msh_colT = msh_precompute(mshT);
    sp_colT = sp_precompute(spaceT,msh_colT, 'value', true, 'gradient', true);

    for idim = 1:mshS.rdim
      x{idim} = repmat(reshape(msh_col.geo_map(idim,:,:),msh_col.nqn,msh_col.nel),[1,1,msh_colT.nqn,msh_colT.nel]);
    end
    x{idim+1} = repmat(reshape(msh_colT.map(msh_colT.qn),1,1,msh_colT.nqn,msh_colT.nel),[msh_col.nqn,msh_col.nel,1,1]);

    errl2  = errl2 + (st_l2_error (sp_col,sp_colT,msh_col,msh_colT,v,vex(x{:}))).^2;
        
    errh1t = errh1t + (st_seminorm_h1t_error (sp_col,sp_colT,msh_col,msh_colT,v,dtvex(x{:}))).^2;

    errh1s = errh1s + (st_seminorm_h1s_error (sp_col,sp_colT,msh_col,msh_colT,v,gradvex(x{:}))).^2;

  end
  
  errh1s = sqrt(errh1s);
  errh1t = sqrt(errh1t);
  errl2  = sqrt (errl2);

end
