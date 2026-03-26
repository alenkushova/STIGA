%% This routine requires 'miga' packege 
% To compute the errors on a the finer grids, which are saved in the file 
% 'errors_on_finer_mesh.mat', we implemented this routine that requires 
% 'miga' packege. To get the 'miga' packege please contact 
%
%  Dr. Andrea Bressan, e-mail: andrea.bressan@imati.cnr.it 
%

interi = [8 16 32 64 128 256];
v_rel_err = zeros(3,6);
for p = 1:3
 for exp = 1:6
  filename = ['ST_SCHRODINGER_FURIER_1D_LUFD_PCG_ps_' num2str(p+1) '_pt_' num2str(p) '_Nt_' num2str(interi(exp)) '_final_time_2.mat']
  load(filename)
  Bt=splineBasis(d-1,T*space.tsp_trial.knots{1});
  Bs=splineBasis(d,space.xsp_trial.knots{1});  
  intd=5; n_el = 256;
  [k,w]=gauss(intd);
  els=knots2elements(Bt.knots,linspace(0,T,2*n_el+1),[0,T]);
  n_elements=numel(els)-1;
  s=diff(els);                 % element size
  m=conv(els,[.5,.5],'valid'); % element center
  t_knots=repmat(k(:),1,n_elements)*diag(s)+m(:)';
  t_weights=repmat(w(:),1,n_elements)*diag(s);
  t_integration=quadrature(t_knots,t_weights);
  els=knots2elements(Bs.knots,linspace(0,1,n_el+1),[0,1]);
  n_elements=numel(els)-1;
  s=diff(els);                 % element size
  m=conv(els,[.5,.5],'valid'); % element center
  s_knots=repmat(k(:),1,n_elements)*diag(s)+m(:)';
  s_weights=repmat(w(:),1,n_elements)*diag(s);
  s_integration=quadrature(s_knots,s_weights);

  [ts,xs]=meshgrid(t_integration.knots,s_integration.knots);

  u_ex_val=uex(xs,ts);
  f_val=f(xs,ts);

  u_l2=(s_integration.weights'*abs(u_ex_val).^2*t_integration.weights)^(1/2);
  u_sl2=(s_integration.weights'*abs(f_val).^2*t_integration.weights)^(1/2);
  u_v= norm([u_l2,u_sl2]);

  u_app=splineFunction([Bs,Bt],u,1);

  u_err=squeeze(u_app.evaluate({s_integration.knots,t_integration.knots}))-u_ex_val;
  l2_err=(s_integration.weights'*abs(u_err).^2*t_integration.weights)^(1/2);

  res=-squeeze(u_app.evaluate({s_integration.knots,t_integration.knots},[2;0]))...
    +1i*squeeze(u_app.evaluate({s_integration.knots,t_integration.knots},[0;1]))...
    -f_val;
  sl2_err=(s_integration.weights'*abs(res).^2*t_integration.weights)^(1/2);
  v_err=norm([l2_err;sl2_err]);

  v_rel_err(p,exp) = v_err/u_v
 end
end

save("errors_on_finer_mesh","v_rel_err")

%% Plot the errors on the finer mesh. This produces exactly Figure 5b as in
% the paper. 
load("errors_on_finer_mesh")

figure ('Units', 'pixels', 'Position', [150 200 500 350])
colors = ["#0072BD", "#D95319" ,"#EDB120" ,"#7E2F8E", "#77AC30"];
ord = {'$$h$$','$$ h^2 $$','$$ h^3 $$'};
symbols = {"-*", "-^", "-o"};

T = 2;  % final time
for grad = 2:4
    err_l2 = []; err_G = []; NN = [];
for i = 1 : 6
    N = T*2^(i+1);
    % load the workspace
    NN = [NN 2/N];
    % save errors: this saves the error i already computed
    err_G   = [err_G v_rel_err(grad-1,i)];
end
    deg = ['$p_s = $' num2str(grad) ', $p_t = $' num2str(grad-1) '.'];
    loglog(NN,err_G,symbols{grad-1},'Color',colors(grad-1),...
        'MarkerFaceColor',colors(grad-1),'Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    clear d s
end

legend('Location','southeast','Interpreter','latex')
title('Error convergence','Interpreter','latex')
xlabel('$$h$$','Interpreter','latex')
ylabel('$$||u-u_h||_{V}/||u||_{V}$$','Interpreter','latex')
