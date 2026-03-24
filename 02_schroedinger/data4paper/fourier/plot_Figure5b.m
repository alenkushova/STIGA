%==========================================================================
% ERRORE IN NORMA GRAFICO PER SOLUZIONE DI FOURIER
%==========================================================================

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
    s = ['ST_SCHRODINGER_FURIER_1D_LUFD_PCG_ps_' num2str(grad) '_pt_' num2str(grad-1) '_Nt_' num2str(N) '_final_time_2.mat'];    
    load(s);
    NN = [NN 2/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 REL_ERR_l2];
    err_G   = [err_G REL_ERR_Graph];
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
