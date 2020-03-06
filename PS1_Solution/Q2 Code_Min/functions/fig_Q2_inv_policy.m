function fig_Q2_inv_policy(k_grid,inv,file_name)
% Create figure for profit function

titlesize   = 24;
labelsize   = 18;
axfontsize  = 16;
legfontsize = 16; 
font        = 'Times New Roman';

f = figure('Units','inches','Position',[0,0,8,6]);

plot(k_grid, inv, 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840]); hold on


set(gca,'FontName',font,'FontSize',axfontsize)

xlabel('Capital $k$','interpreter','latex','FontSize',labelsize,'FontName',font)
ylabel('Investment $i(k)$','interpreter','latex','FontSize',labelsize,'FontName',font)
grid; box off; axis tight;

title('Optimal Investment Policy','interpreter','latex','FontSize',titlesize,'FontName',font)
legend('$i(k)$ for $a=1$','Location','best','interpreter','latex','FontSize',legfontsize)

print(f,'-depsc','-painters','-noui','-r600',file_name)


end

