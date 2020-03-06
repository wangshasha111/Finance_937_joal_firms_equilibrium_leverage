function fig_Q2_value_func(k_grid,value,file_name)
% Create figure for profit function

titlesize   = 24;
labelsize   = 18;
axfontsize  = 16;
legfontsize = 16; 
font        = 'Times New Roman';

f = figure('Units','inches','Position',[0,0,8,6]);

plot(k_grid, value, 'LineWidth', 3); hold on


set(gca,'FontName',font,'FontSize',axfontsize)

xlabel('Capital $k$','interpreter','latex','FontSize',labelsize,'FontName',font)
ylabel('Value function $v(k)$','interpreter','latex','FontSize',labelsize,'FontName',font)
grid; box off; axis tight;

title('Value function','interpreter','latex','FontSize',titlesize,'FontName',font)
legend('$v(k)$ for $a=1$','Location','best','interpreter','latex','FontSize',legfontsize)

print(f,'-depsc','-painters','-noui','-r600',file_name)


end

