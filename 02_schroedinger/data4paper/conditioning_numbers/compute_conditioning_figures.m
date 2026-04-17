clear; close all; clc;

% plot scatter eigs for p and nel = 32 in dim 1.
f1 = figure();
f1.WindowState = 'maximized';

cond_tab = zeros(9,5);
for p = 2:10
  for nel = [8 16 32 64 128]
    [lam, cond_num] = conditioning(p, nel);
    cond_tab(p-1,log2(nel)-2) = cond_num;
    if nel == 32
      % 1d scatter-plot
      tit = ['$p_s = $ ' num2str(p)];
      subplot(3,3,p-1)
      scatter(1:numel(lam.dim1), lam.dim1,'filled')
      title(tit,'FontSize', 14,'Interpreter','latex')
      grid on
      ylim([1 (p+1)/2]);
      xlim([1 40])
      fontsize(f1, 14, "points")
    end
  end
end

% plot of conditioning number of [(L(M)-1L)-1]*B for different p and nel.
f2 = figure();
f2.WindowState = 'maximized';
subplot(1,2,1)
plot(2.^(3:7),cond_tab(1,:),'LineWidth',1.5,'LineStyle','-')
hold on
plot(2.^(3:7),cond_tab(2,:),'LineWidth',1.5,'LineStyle','--')
plot(2.^(3:7),cond_tab(3,:),'LineWidth',1.5,'LineStyle',':')
plot(2.^(3:7),cond_tab(4,:),'LineWidth',1.5,'LineStyle','-','Marker','square')
plot(2.^(3:7),cond_tab(5,:),'LineWidth',1.5,'LineStyle','--','Marker','square')
plot(2.^(3:7),cond_tab(6,:),'LineWidth',1.5,'LineStyle',':','Marker','square')
plot(2.^(3:7),cond_tab(7,:),'LineWidth',1.5,'LineStyle','-','Marker','*')
plot(2.^(3:7),cond_tab(8,:),'LineWidth',1.5,'LineStyle','--','Marker','*')
plot(2.^(3:7),cond_tab(9,:),'LineWidth',1.5,'LineStyle',':','Marker','*')
xlabel('$n_{el}$','FontSize', 14,'interpreter','latex');
ylabel('$\kappa \left( (L_s M_s^{-1} L_s)^{-1} B_s \right)$','FontSize', 14,'interpreter','latex');
legend('$p_s=2$','$p_s=3$','$p_s=4$','$p_s=5$','$p_s=6$','$p_s=7$','$p_s=8$','$p_s=9$','$p_s=10$','interpreter','latex',Location='best')
grid('on')
hold off
subplot(1,2,2)
plot(2:10,cond_tab(:,1),'LineWidth',1.5,'LineStyle','-')
hold on
plot(2:10,cond_tab(:,2),'LineWidth',1.5,'LineStyle','--')
plot(2:10,cond_tab(:,3),'LineWidth',1.5,'LineStyle',':')
plot(2:10,cond_tab(:,4),'LineWidth',1.5,'LineStyle','-','Marker','square')
plot(2:10,cond_tab(:,5),'LineWidth',1.5,'LineStyle','--','Marker','square')
xlabel('$p_s$','FontSize', 14,'interpreter','latex');
ylabel('$ \kappa \left((L_s M_s^{-1} L_s)^{-1} B_s \right)$','FontSize', 14,'interpreter','latex');
legend('$n_{el} = 8$','$n_{el} = 16$','$n_{el} = 32$','$n_{el} = 64$','$n_{el} = 128$','interpreter','latex','Location','northwest');
grid('on')
fontsize(f2, 14, "points")

%% plot scatter eigs for p and nel = 32 in dim 2.
f3 = figure();
f3.WindowState = 'maximized';
cond_tab = zeros(9,5);
for p = 2:10
  nel = 32;
  [lam, cond_num] = conditioning(p, nel);
  % 2d scatter-plot
  tit = ['$p_s = $ ' num2str(p)];
  subplot(3,3,p-1)
  scatter(1:numel(lam.dim2), lam.dim2,'filled')
  title(tit,'FontSize', 14,'Interpreter','latex')
  grid on
  ylim([1 (p+1)/2]);
  xlim([1 1600])
  fontsize(f3, 14, "points")
end