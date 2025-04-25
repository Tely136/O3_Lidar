function NDplot2(ND,TimeInHour_avg,hkm,title_str,y_limit,z_limit,fontsize) 
% use TOLNET colorbar 
figure
I=imagesc(TimeInHour_avg,hkm,ND,z_limit);
set(gca,'YDir','normal','FontSize',fontsize);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(ND));

% colormap('jet');
 cbar_tolnet=colorbar_tolnet_func;
 colormap(cbar_tolnet);
 
colorbar;
xlabel('Time (hour UTC)');
ylabel('Altitude (km)');
title(title_str);
ylim(y_limit);
set(gca, 'TickDir', 'out');