
I=imagesc(TimeInHour_avg,hkm,pz2); hold on
plot(TimeInHour_avg,cldCenterZ,'o');
set(gca,'YDir','normal');
set(gca,'FontSize',14);
set(gca,'ColorScale','log');
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone lidar smoothed 299nm Pz2(a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

figure
I=imagesc(TimeInHour_avg,hkm(1:end-1),d_Pz2);
set(gca,'YDir','normal');
set(gca,'FontSize',14);
set(gca,'ColorScale','log');
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone lidar smoothed 287nm P derivative (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

id=10;
figure
plot(pz2(:,id),hkm,'b');
hold on
plot(pz2(cldCenterZ_ind(:,id),id),cldCenterZ(:,id),'o');
set(gca,'XScale','log');
legend('299nm Far')
xlabel('Range corrected signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime_avg(id))])
ylim([0,20])


