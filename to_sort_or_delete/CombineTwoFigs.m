% Load saved figures
ind_gd_ozone=isbetween(ground_O3_ppb_datetime,DateTime_avg(1)-minutes(10),DateTime_avg(end)+minutes(10));
ind_sh_ozone=isbetween(shed_O3_ppb_datetime,DateTime_avg(1)-minutes(5),DateTime_avg(end)+minutes(5));

indh=find(hkm_nr>0.25,1,'first');
% Prepare subplots
figure
subplot(2,1,1)
NDplot(O3ppb_sm,TimeInHour_avg,hkm_nr,title_str,[0.24,6],[0,100],14);
subplot(2,1,2)
plot(shed_O3_ppb_datetime(ind_sh_ozone),shed_O3_ppb(ind_sh_ozone),'r-');
hold on;
plot(DateTime_avg,O3ppb(indh,:),'bo-');
plot(DateTime_avg,70*ones(1,24),'g--');
set(gca,'LineWidth',1,'FontSize',14);
legend('in-situ','O_3-DIAL measurement at 250 m');
xlabel('EDT Time (Hour)');
ylabel('Ozone Mixing Ratio (ppb)');
ylim([40,100]);xlim([DateTime_avg(1) DateTime_avg(end)])
% title(['DIAL and Ground Ozone Measurements'])
grid on;