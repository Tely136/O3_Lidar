function plot_fr_nr_Prof(prof_287,prof_299,prof_287_nr,prof_299_nr,hkm_fr,hkm_nr,DateTime_avg,id,tit_str,xlb,y_lim,le)
figure
plot(prof_287(:,id),hkm_fr,'r');hold on
plot(prof_299(:,id),hkm_fr,'b');
plot(prof_287_nr(:,id),hkm_nr,'g');
plot(prof_299_nr(:,id),hkm_nr,'k');
set(gca,'XScale','log')
legend(le)
xlabel(xlb);
ylabel('Altitude (km)');
title([tit_str,datestr(DateTime_avg(id))])
ylim(y_lim)
grid on