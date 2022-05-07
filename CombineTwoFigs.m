% Load saved figures
ind_gd_ozone=isbetween(ground_O3_ppb_datetime,DateTime_avg(1)-minutes(10),DateTime_avg(end)+minutes(10));
indh=find(hkm_nr>0.21,1,'first');
% Prepare subplots
figure
subplot(2,1,1)
NDplot(O3ppb_sm,TimeInHour_avg,hkm_nr,title_str,[0.35,8],[0,100],14);
subplot(2,1,2)
plot(ground_O3_ppb_datetime(ind_gd_ozone),ground_O3_ppb(ind_gd_ozone),'r-');
hold on;
plot(DateTime_avg,O3ppb(indh,:),'bo-');
set(gca,'LineWidth',1,'FontSize',14);
legend('in-situ','O_3-DIAL measurement at 210m');
xlabel('EDT Time (Hour)');
ylabel('Ozone Mixing Ratio (ppb)');
ylim([0,80]);
% title(['DIAL and Ground Ozone Measurements'])
grid on;