%figure ('Units', 'pixels', 'Position', [150 200 1000 350])
%==========================================================================
% ERRORE IN NORMA GRAFICO PER SOLUZIONE SMOOTH
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
    s = ['ST_SCHRODINGER_SMOOTH_1D_LUFD_PCG_ps_' num2str(grad) '_pt_' num2str(grad-1) '_Nt_' num2str(N) '_final_time_2.mat'];    
    load(s);
    NN = [NN 2/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 REL_ERR_l2];
    err_G   = [err_G REL_ERR_Graph];
end
    deg = ['$p_s = $ ' num2str(grad) ', $p_t = $ ' num2str(grad-1) '.'];
    loglog(NN,err_G,symbols{grad-1},'Color',colors(grad-1),...
        'MarkerFaceColor',colors(grad-1),'Linewidth',1.5,'DisplayName',deg)
    % loglog(NN,err_l2,'-s','Color',colors(grad-1),...
    %     'MarkerFaceColor',colors(grad-1),'Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    % loglog(NN,min(err_l2)*NN.^(grad+0.5)/min(NN)^(grad+0.5),'-.k','Linewidth',1.5,'DisplayName',ord)
    [m,j] = min(err_G);
    loglog(NN,m*NN.^(grad-1)/NN(i)^(grad-1)/2,'-.','Color',colors(grad-1),'Linewidth',1.5,'DisplayName', ord{grad-1})
    clear d s
end
% loglog(NN,3*NN.^(1/5),'-.','Color','black','Linewidth',1.5,'DisplayName','$$h^{0.2}$$')

legend('Location','southeast','Interpreter','latex')
% title('Error convergence','Interpreter','latex')
xlabel('$$h$$', 'Interpreter','latex')
ylabel('$$||u-u_h||_{V}/||u||_{V}$$','Interpreter','latex')
ylim([10^(-12) 1])
fontsize(14,'points')
legend(FontSize=12)

%%
%==========================================================================
% SOLO NORMA L^2 PER SOLUZIONE SMOOTH
%==========================================================================

figure ('Units', 'pixels', 'Position', [150 200 500 350])
colors = ["#0072BD", "#D95319" ,"#EDB120" ,"#7E2F8E", "#77AC30"];
ord = {'$$h^{1.5}$$','$$ h^{3.5} $$','$$ h^{4.5} $$'};
symbols = {"-*", "-^", "-o"};
for grad = 2:4
    err_l2 = []; err_G = []; NN = [];
for i = 1 : 6
    N = T*2^(i+1);
    % load the workspace
    s = ['ST_SCHRODINGER_SMOOTH_1D_LUFD_PCG_ps_' num2str(grad) '_pt_' num2str(grad-1) '_Nt_' num2str(N) '_final_time_2.mat'];    
    load(s);
    NN = [NN 2/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 REL_ERR_l2];
    err_G   = [err_G REL_ERR_Graph];
end
    deg = ['$p_s = $' num2str(grad) ', $p_t = $' num2str(grad-1) '.'];
    % loglog(NN,err_G,symbols{grad-1},'Color',colors(grad-1),...
    %     'MarkerFaceColor',colors(grad-1),'Linewidth',1.5,'DisplayName',deg)
    loglog(NN,err_l2,symbols{grad-1},'Color',colors(grad-1), 'MarkerFaceColor',colors(grad-1),'Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    m = (err_l2(1));
    if grad == 2
        loglog(NN,NN.^(grad-0.5)*m*10^(1),'-.','Color',colors(grad-1),'Linewidth',1.5,'DisplayName', ord{grad-1})        
    elseif grad == 3 
        loglog(NN,NN.^(grad+0.5)*m*10^(grad-1),'-.','Color',colors(grad-1),'Linewidth',1.5,'DisplayName', ord{grad-1})
    else 
        loglog(NN,NN.^(grad+0.5)*m*10^(grad-1.5),'-.','Color',colors(grad-1),'Linewidth',1.5,'DisplayName', ord{grad-1})
    end
    clear d s
end
%loglog(NN,3*NN.^(1/5),'-.','Color','black','Linewidth',1.5,'DisplayName','$$h^{0.2}$$')

legend('Location','southeast','Interpreter','latex')
%title('Error convergence','Interpreter','latex')
xlabel('$$h$$','Interpreter','latex')
ylabel('$$||u-u_h||_{L^2(\mathcal{Q})}/||u||_{L^2(\mathcal{Q})}$$','Interpreter','latex')
ylim([10^(-14) 10^(-2)])
fontsize(14,'points')
legend(FontSize=12)