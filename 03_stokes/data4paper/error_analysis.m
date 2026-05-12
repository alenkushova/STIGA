close all
fig=figure();
fig.WindowState = 'maximized';
colors = ["#0072BD", "#D95319" ,"#EDB120" ,"#7E2F8E", "#77AC30"];
tag_degrees = {'$$p = 1$$','$$p = 2$$','$$p = 3$$'};
tag_orders = {'$$h^2$$','$$h^3$$','$$h^4$$'};
tag_orders_lower = {'$$h^1$$','$$h^2$$','$$h^3$$'};

for d = 1:3 
ERR_PRES_L2 = []; ERR_VEL_L2  = [];
ERR_VEL_H1s = []; ERR_VEL_H1t = [];  
mesh = [];

for i = 1:4
    filname = ['Explicit_solution_TH_gmres_results_d=' num2str(d) '_n=' num2str(2^(i+1)) '_T=1.mat'];
    load(filname);

    % Errori assoluti
    % ERR_PRES_L2 = [ERR_PRES_L2 pres_errl2];
    % ERR_VEL_L2  = [ERR_VEL_L2  vel_errl2];
    % ERR_VEL_H1s = [ERR_VEL_H1s vel_errh1s];
    % ERR_VEL_H1t = [ERR_VEL_H1t vel_errh1t];  

    % Errrori relativi
    ERR_PRES_L2 = [ERR_PRES_L2 pres_errl2_rel];
    ERR_VEL_L2  = [ERR_VEL_L2  vel_errl2_rel];
    ERR_VEL_H1s = [ERR_VEL_H1s vel_errh1s_rel];
    ERR_VEL_H1t = [ERR_VEL_H1t vel_errh1t_rel];      
    mesh = [mesh 1/(2^(i+1))];
end

subplot(1,3,1)
loglog(mesh,ERR_PRES_L2,'Marker','square','LineWidth',1.5,'Color',colors(d),...
    'MarkerFaceColor',colors(d),'DisplayName',tag_degrees{d})
hold on
grid on
loglog(mesh,mesh.^(d+1)/mesh(end)^(d+1)*ERR_PRES_L2(end)*(0.9),'--',...
    'LineWidth',1.5,'Color',colors(d),'DisplayName',tag_orders{d})

subplot(1,3,2)
loglog(mesh,ERR_VEL_H1s,'Marker','square','LineWidth',1.5,'Color',colors(d),...
    'MarkerFaceColor',colors(d),'DisplayName',tag_degrees{d})
hold on
grid on
loglog(mesh,mesh.^(d+1)/mesh(end)^(d+1)*ERR_VEL_H1s(end)*(0.9),'--',...
    'LineWidth',1.5,'Color',colors(d),'DisplayName',tag_orders{d})

subplot(1,3,3)
loglog(mesh,ERR_VEL_H1t,'Marker','square','LineWidth',1.5,'Color',colors(d),...
    'MarkerFaceColor',colors(d),'DisplayName',tag_degrees{d})
hold on
grid on
loglog(mesh,mesh.^(d)/mesh(end)^(d)*ERR_VEL_H1t(end)*(0.9),'--',...
    'LineWidth',1.5,'Color',colors(d),'DisplayName',tag_orders_lower{d})

end

% Errori assoluti: 
% subplot(1,3,1)
% legend('Location','southeast','Interpreter','latex')
% title('Error convergence','Interpreter','latex')
% xlabel('$$h$$','Interpreter','latex')
% ylabel('$$||p-p_h||_{L^2(Q)}$$','Interpreter','latex')
% 
% subplot(1,3,2)
% legend('Location','southeast','Interpreter','latex')
% title('Error convergence','Interpreter','latex')
% xlabel('$$h$$','Interpreter','latex')
% ylabel('$$||\partial_t (\mathbf{v}-\mathbf{v}_h)||_{L^2(Q)}$$','Interpreter','latex')
% 
% subplot(1,3,3)
% legend('Location','southeast','Interpreter','latex')
% title('Error convergence','Interpreter','latex')
% xlabel('$$h$$','Interpreter','latex')
% ylabel('$$||\nabla (\mathbf{v}-\mathbf{v}_h)||_{L^2(Q)}$$','Interpreter','latex')

%% Error relarivi
subplot(1,3,1)
legend('Location','southeast','Interpreter','latex')
title('Error convergence','Interpreter','latex')
xlabel('$$h$$','Interpreter','latex')
ylabel('$$||p-p_h||_{L^2(Q)}/||p||_{L^2(Q)}$$','Interpreter','latex')

subplot(1,3,2)
legend('Location','southeast','Interpreter','latex')
title('Error convergence','Interpreter','latex')
xlabel('$$h$$','Interpreter','latex')
ylabel('$$||\nabla (\mathbf{v}-\mathbf{v}_h)||_{L^2(Q)}/||\nabla \mathbf{v}||_{L^2(Q)}$$','Interpreter','latex')

subplot(1,3,3)
legend('Location','southeast','Interpreter','latex')
title('Error convergence','Interpreter','latex')
xlabel('$$h$$','Interpreter','latex')
ylabel('$$||\partial_t (\mathbf{v}-\mathbf{v}_h)||_{L^2(Q)}/||\partial_t \mathbf{v}||_{L^2(Q)}$$','Interpreter','latex')
