fstruct = dir('pc_ad_*.mat');
numf=length(fstruct);
legstr287=[];
legstr299=[];
for i=1:numf
    fnames=fstruct(i).name;
    r=load(fnames);
    an287(:,i)=r.pc_ad_ratio.an287;
    an299(:,i)=r.pc_ad_ratio.an299;
    pc287(:,i)=r.pc_ad_ratio.pc287;
    pc299(:,i)=r.pc_ad_ratio.pc299;
    hv287=r.pc_ad_ratio.HV287;
    hv299=r.pc_ad_ratio.HV299;
    d287=r.pc_ad_ratio.disc287;
    d299=r.pc_ad_ratio.disc299;
    legstr287=[legstr287;'287nm PC,Far - ',num2str(hv287,'%d'),'V -disc ',num2str(d287,'%d')];
    legstr299=[legstr299;'299nm PC,Far - ',num2str(hv299,'%d'),'V -disc ',num2str(d299,'%d')];
end
figure
plot(hkm,an287(:,1));hold on;
plot(hkm,an299(:,1));
ylabel('Analog signal (mV)');
ylim([1e-4,500])
set(gca,'YScale','log')
yyaxis right
plot(hkm,pc287,'b--');hold on;
plot(hkm,pc299,'r--');
set(gca,'YScale','log')
ylabel('PC signal (MHz)');
xlabel('Altitude (km)');
xlim([0,22.5])
%ylim([1e-3,500])
% legendstr=['287nm - AD, Far, HV=800V';'299nm - AD, Far, HV=795V';legstr287;legstr299]
% legend('287nm - AD, Far, HV=800V','299nm - AD, Far, HV=795V','287nm - PC, Far, HV=800, disc=4','299nm - PC, Far,  HV=795, disc=4')
title(['Ozone lidar signal (10min-Avg) at ',datestr(DateTime_avg(id))])


