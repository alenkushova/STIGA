clear; close all; clc;

%==========================================================================
% ERRORE IN semi-NORMA H1 PER ALPHA = 1.5
%==========================================================================

figure ('Units', 'pixels', 'Position', [150 200 500 350])
colors = ["#0072BD", "#D95319" ,"#EDB120" ,"#7E2F8E", "#77AC30"];
ord = {'$$h^{0.5}$$','$$ h^{0.05} $$'};
symbols = {"-*", "-^", "-o","-x"};
alpha = 1.5;

for degree = 1:4
  err_l2 = []; err_h1s = []; hh = [];
  for i = 1 : 7
    N = 2^(i+1);
    filename = ['ST_ps_' num2str(degree) '_pt_' num2str(degree) '_alpha_' num2str(alpha) '_N_' num2str(N) '.mat'];
    load(filename)
    hh = [hh 1/N];
    err_l2  = [err_l2 report.rel_errl2];
    err_h1s = [err_h1s report.rel_errh1s];
  end
  deg = ['$p_s = $ ' num2str(degree) ', $p_t = $ ' num2str(degree) '.'];
  loglog(hh,err_h1s,symbols{degree},'Color',colors(degree),...
        'MarkerFaceColor',colors(degree),'Linewidth',1.5,'DisplayName',deg)
  grid on, hold on
end
loglog(hh,hh.^(alpha-1)/4,'-.','Color','black','Linewidth',1.5,'DisplayName', ord{1})
legend('Location','best','Interpreter','latex')
title('$u(x,t) = |x|^\alpha e^{-t}$, $\alpha = 1.5$', Interpreter="latex")
xlabel('$h$',Interpreter='latex')
ylabel('$||\nabla (u-u_h) ||_{L^2(\Omega)}/||\nabla u||_{L^2(\Omega)}$',Interpreter='latex')
fontsize(14,'points')

%%
%==========================================================================
% ERRORE IN semi-NORMA H1 PER ALPHA = 1.05
%==========================================================================

figure ('Units', 'pixels', 'Position', [150 200 500 350])
colors = ["#0072BD", "#D95319" ,"#EDB120" ,"#7E2F8E", "#77AC30"];
ord = {'$$h^{0.5}$$','$$ h^{0.05} $$'};
symbols = {"-*", "-^", "-o","-x"};
alpha = 1.05;

for degree = 1:4
  err_l2 = []; err_h1s = []; hh = [];
  for i = 1 : 7
    N = 2^(i+1);
    filename = ['ST_ps_' num2str(degree) '_pt_' num2str(degree) '_alpha_' num2str(alpha) '_N_' num2str(N) '.mat'];
    load(filename)
    hh = [hh 1/N];
    err_l2  = [err_l2 report.rel_errl2];
    err_h1s = [err_h1s report.rel_errh1s];
  end
  deg = ['$p_s = $ ' num2str(degree) ', $p_t = $ ' num2str(degree) '.'];
  loglog(hh,err_h1s,symbols{degree},'Color',colors(degree),...
        'MarkerFaceColor',colors(degree),'Linewidth',1.5,'DisplayName',deg)
  grid on, hold on
end
loglog(hh,hh.^(alpha-1)/(1.2),'-.','Color','black','Linewidth',1.5,'DisplayName', ord{2})
legend('Location','best','Interpreter','latex')
title('$u(x,t) = |x|^\alpha e^{-t}$, $\alpha = 1.05$', Interpreter="latex")
xlabel('$h$',Interpreter='latex')
ylabel('$||\nabla (u-u_h) ||_{L^2(\Omega)}/||\nabla u||_{L^2(\Omega)}$',Interpreter='latex')
fontsize(14,'points')
