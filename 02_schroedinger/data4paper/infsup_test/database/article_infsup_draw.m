clear all 
close all

figure
load('gal_degree_2_infsup_trial_graf_test_l2_norms.mat');
loglog(h,infsup,'LineWidth',1.5);
grid on; hold on; 
xlabel('$$ h $$',Interpreter='latex'); ylabel('$$ \alpha_h $$',Interpreter='latex');
load('gal_degree_3_infsup_trial_graf_test_l2_norms.mat');
loglog(h,infsup,'--','LineWidth',1.5);
load('gal_degree_4_infsup_trial_graf_test_l2_norms.mat');
loglog(h,infsup,'-.','LineWidth',1.5);
legend('$$ p = 2 $$','$$ p = 3 $$','$$ p = 4 $$','Interpreter','Latex','location','southeast');
hold off
fontsize(14, "points")

figure
load('ls_degree_2_infsup_trial_graf_test_l2_norms.mat');
loglog(h,infsup,'LineWidth',1.5);
grid on; hold on; 
xlabel('$$ h $$',Interpreter='latex'); ylabel('$$ \alpha_h $$',Interpreter='latex');
load('ls_degree_3_infsup_trial_graf_test_l2_norms.mat');
loglog(h,infsup,'--','LineWidth',1.5);
load('ls_degree_4_infsup_trial_graf_test_l2_norms.mat');
loglog(h,infsup,'-.','LineWidth',1.5);
legend('$$ p = 2 $$','$$ p = 3 $$','$$ p = 4 $$','Interpreter','Latex','location','southeast');
hold off
ylim([0.85, 1])
fontsize(14, "points")

