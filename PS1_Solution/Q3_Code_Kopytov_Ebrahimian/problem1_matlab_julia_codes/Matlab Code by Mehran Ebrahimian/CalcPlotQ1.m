clc
close all
clear all


load('data_Q1.mat')

Stat_Dist = mu;

%%

figure(1)
plot(k_vec,Policy_k(:,[1,(Nzg+1)/2,end]),'linewidth',1.5)
hold on
plot(0:50,0:50,'--k')
% title(['\it Q2, Part ',Part,': Policy Function'],'interpreter','Latex','fontsize',15)
xlabel('$k$','interpreter','Latex','fontsize',14)
ylabel('$k^{`}(k)$','interpreter','Latex','fontsize',14)
legend({[repmat('$z=',3,1),['low ';'med.';'high'],repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
box on
xlim([0,20])
saveas(gcf,'graph_Q1_policy.eps','epsc')
% saveas(gcf,'graph_Q1_policy.fig','fig')


%%
figure(2)
plot(k_vec,V(:,[1,(Nzg+1)/2,end]),'linewidth',1.5)
% title(['\it Q2, Part ',Part,': Policy Function'],'interpreter','Latex','fontsize',15)
xlabel('$k$','interpreter','Latex','fontsize',14)
ylabel('$V(k)$','interpreter','Latex','fontsize',14)
legend({[repmat('$z=',3,1),['low ';'med.';'high'],repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
box on
xlim([0,20])
saveas(gcf,'graph_Q1_V.eps','epsc')
% saveas(gcf,'graph_Q1_V.fig','fig')

%%


figure(3)
surf(z_vec,k_vec,Stat_Dist,'linewidth',1.5)
% title(['\it Q2, Part ',Part,': Policy Function'],'interpreter','Latex','fontsize',15)
xlabel('$z$','interpreter','Latex','fontsize',14)
ylabel('$k$','interpreter','Latex','fontsize',14)
zlabel('$\mu(k,z)$','interpreter','Latex','fontsize',14)
% legend({[repmat('$z=',3,1),num2str(z_vec([1,(Nzg+1)/2,end])'),repmat('$',3,1)]},'interpreter','Latex','location','best','fontsize',13)
box on
% saveas(gcf,'graph_Q1_Dist.eps','epsc')
% saveas(gcf,'graph_Q1_Dist.fig','fig')




%% display stuffs
Entry_Cost = sum(station_state_z'.*V(1,:));

disp('entry cost: ')
disp(Entry_Cost)

% 
disp('% of Exit Rate : ')
disp(Exit_rate*100)


%% labor demand

L_d_policy = (alpha_l*z_mat.*k_mat.^(alpha_k)/W).^(1/(1-alpha_l));
L_d = sum(sum(L_d_policy.*mu));

disp('total labor demand is: ')
disp(L_d)

B = L_d / W^0.5