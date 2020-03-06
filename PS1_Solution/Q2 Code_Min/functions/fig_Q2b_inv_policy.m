function fig_Q2b_inv_policy(k_grid,inv,file_name)
% Create figure for profit function

titlesize   = 24;
labelsize   = 18;
axfontsize  = 16;
legfontsize = 16; 
font        = 'Times New Roman';

f = figure('Units','inches','Position',[0,0,8,6]);

plot(k_grid, inv(:,1), 'LineWidth', 3); hold on
plot(k_grid, inv(:,2),'--', 'LineWidth', 3); hold on
plot(k_grid, inv(:,3), '-.','LineWidth', 3); hold on


set(gca,'FontName',font,'FontSize',axfontsize)

xlabel('Capital $k$','interpreter','latex','FontSize',labelsize,'FontName',font)
ylabel('Investment $i(a,k)$','interpreter','latex','FontSize',labelsize,'FontName',font)
grid; box off; axis tight;

title('Optimal Investment Policy','interpreter','latex','FontSize',titlesize,'FontName',font)
legend('$i(a=0.9, k)$','$i(a=1, k)$','$i(a=1.1, k)$','Location','best','interpreter','latex','FontSize',legfontsize)

print(f,'-depsc','-painters','-noui','-r600',file_name)


end

