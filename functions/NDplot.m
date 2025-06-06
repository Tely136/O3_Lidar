function NDplot(ND,TimeInHour_avg,hkm,title_str,y_limit,z_limit,fontsize)
figure
I=imagesc(TimeInHour_avg,hkm,ND,z_limit);
set(gca,'YDir','normal','FontSize',fontsize);
%set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(ND))
colormap(jet);
colorbar
xlabel('Time (hour UTC)')
ylabel('Altitude (km)')
title(title_str);
ylim(y_limit);
set(gca, 'TickDir', 'out');
