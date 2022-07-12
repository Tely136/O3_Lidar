function plot_cld_screen_result(DateTime_avg,TimeInHour_avg,hkm,pz2,d_Pz2,cldBaseZ,cldBaseZ_ind,id)
%% cloud screen test  result
figure
subplot(1,2,1)
I=imagesc(TimeInHour_avg,hkm,pz2); hold on
plot(TimeInHour_avg,cldBaseZ(1,:),'o');
set(gca,'YDir','normal','FontSize',14,'ColorScale','log');colormap('jet');colorbar
xlabel('Local Time (hour)');ylabel('Altitude (km)')
title(['Ozone lidar smoothed Pz2(a.u.) at ',datestr(DateTime_avg(1),'yy/mm/dd')]);

subplot(1,2,2)
I=imagesc(TimeInHour_avg,hkm(1:end-1),d_Pz2);
set(gca,'YDir','normal','FontSize',14,'ColorScale','log');colormap('jet');colorbar
xlabel('Local Time (hour)');ylabel('Altitude (km)')
title(['Ozone lidar smoothed P derivative (a.u.) at ',datestr(DateTime_avg(1),'yy/mm/dd')]);

% id = find(cldBaseZ(1,:)<10&cldBaseZ(1,:)>1);
% subplot(1,3,3)
% plot(pz2(:,id),hkm,'b');hold on
% plot(pz2(cldBaseZ_ind(1,id),id),cldBaseZ(1,id),'o');
% set(gca,'XScale','log');legend('299nm Far','Cloud Base');ylim([0,20]);
% xlabel('Range corrected signal (a.u.)');ylabel('Altitude (km)');
% title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime_avg(id))])