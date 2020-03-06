function fig_Q2b_kDrift(k_grid,kDrift,file_name)
% Create figure for profit function

titlesize   = 24;
labelsize   = 18;
axfontsize  = 16;
legfontsize = 16; 
font        = 'Times New Roman';

f = figure('Units','inches','Position',[0,0,8,6]);

plot(k_grid, kDrift(:,1), 'LineWidth', 3); hold on
plot(k_grid, kDrift(:,2),'--', 'LineWidth', 3); hold on
plot(k_grid, kDrift(:,3), '-.','LineWidth', 3); hold on


set(gca,'FontName',font,'FontSize',axfontsize)

xlabel('Capital $k$','interpreter','latex','FontSize',labelsize,'FontName',font)
ylabel('Capital Drift $\dot{k}(a,k)$','interpreter','latex','FontSize',labelsize,'FontName',font)
grid; box off; axis tight;

title('Optimal Capital Drift Policy','interpreter','latex','FontSize',titlesize,'FontName',font)
legend('$\dot{k}(a=0.9, k)$','$\dot{k}(a=1, k)$','$\dot{k}(a=1.1, k)$','Location','best','interpreter','latex','FontSize',legfontsize)

print(f,'-depsc','-painters','-noui','-r600',file_name)


end

