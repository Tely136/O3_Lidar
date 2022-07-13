function NDplot(ND,TimeInHour_avg,hkm,title_str,y_limit,z_limit,fontsize)
figure
I=uimagesc(TimeInHour_avg,hkm,ND,z_limit);
set(gca,'YDir','normal','FontSize',fontsize);
% set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(ND))
colormap('jet')
colorbar
xlabel('EDT Time (hour)')
ylabel('Altitude (km)')
title(title_str);
ylim(y_limit)