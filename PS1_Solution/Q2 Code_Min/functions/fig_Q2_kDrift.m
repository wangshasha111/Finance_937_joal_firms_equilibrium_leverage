function fig_Q2_kDrift(k_grid,kDrift,file_name)
% Create figure for profit function

titlesize   = 24;
labelsize   = 18;
axfontsize  = 16;
legfontsize = 16; 
font        = 'Times New Roman';

f = figure('Units','inches','Position',[0,0,8,6]);

plot(k_grid, kDrift, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]); hold on


set(gca,'FontName',font,'FontSize',axfontsize)

xlabel('Capital $k$','interpreter','latex','FontSize',labelsize,'FontName',font)
ylabel('Capital policy $\dot{k}(k)$','interpreter','latex','FontSize',labelsize,'FontName',font)
grid; box off; axis tight;

title('Optimal Capital Drift Policy','interpreter','latex','FontSize',titlesize,'FontName',font)
legend('$\dot{k}(k)$ for $a=1$','Location','best','interpreter','latex','FontSize',legfontsize)

print(f,'-depsc','-painters','-noui','-r600',file_name)


end

