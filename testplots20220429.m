% compare analog and pc signal
id=35;
maxhid = 2054;
figure
subplot(2,1,1)
plot(hkm_fr,prof_pc_287(:,id),'r');hold on
plot(hkm_fr,prof_pc_299(:,id),'b');
plot(hkm_nr(1:maxhid),prof_pc_287_nr(1:maxhid,id),'g');
plot(hkm_nr(1:maxhid),prof_pc_299_nr(1:maxhid,id),'k');
set(gca,'YScale','log')
set(gca,'FontSize',14)
legend('287nm Far','299nm Far','287nm Near','299nm Near')
ylabel('PC Signal (MHz)');
xlabel('Altitude (km)');
title(['Ozone PC signal at ',datestr(DateTime_avg(id))])
xlim([0,15])
%ylim([0,100])
grid on
grid minor

figure
subplot(1,2,1)
plot(prof_merge_287(:,13),hkm_fr,'r');hold on
plot(prof_merge_299(:,13),hkm_fr,'b');
plot(prof_merge_287_nr(1:maxhid,13),hkm_nr(1:maxhid),'m');
plot(prof_merge_299_nr(1:maxhid,13),hkm_nr(1:maxhid),'g');
set(gca,'XScale','log')
set(gca,'FontSize',14)
legend('287nm Far','299nm Far','287nm Near','299nm Near')
xlabel('Signal (MHz)');
ylabel('Altitude (km)');
title(['Ozone AD-PC merged Signal at ',datestr(DateTime_avg(13))])
ylim([0,15])
%ylim([0,100])
grid on
grid minor



subplot(2,1,2)
plot(hkm_fr,prof_an_287(:,id),'r');hold on
plot(hkm_fr,prof_an_299(:,id),'b');
plot(hkm_nr(1:maxhid),prof_an_287_nr(1:maxhid,id),'g');
plot(hkm_nr(1:maxhid),prof_an_299_nr(1:maxhid,id),'k');
set(gca,'YScale','log')
set(gca,'FontSize',14,'LineWidth',1.2)
legend('287nm Far','299nm Far','287nm Near','299nm Near')
ylabel('AD Signal (mV)');
xlabel('Altitude (km)');
title(['Ozone AD signal at ',datestr(DateTime_avg(id))])
xlim([0,15])
%ylim([0,100])
grid on
grid minor

id=13;
figure
plot(hkm_fr,prof_pc_287(:,id)./prof_an_287(:,id),'r');hold on
plot(hkm_fr,prof_pc_299(:,id)./prof_an_299(:,id),'b');
plot(hkm_nr(1:1054),prof_pc_287_nr(1:1054,id)./prof_an_287_nr(1:1054,id),'g');
plot(hkm_nr(1:1054),prof_pc_299_nr(1:1054,id)./prof_an_299_nr(1:1054,id),'k');
% set(gca,'YScale','log')
set(gca,'FontSize',14)
legend('287nm Far','299nm Far','287nm Near','299nm Near')
ylabel('Signal (a.u.)');
xlabel('Altitude (km)');
title(['Ozone PC/AD signal ratio at ',datestr(DateTime_avg(id))])
xlim([0,12])
ylim([0,100])
grid on
grid minor

%% Compare the ozone ND of AD and PC
N_O3_fr_sm = sgolayfilt(movmean(N_O3_fr,[5,15],1),1,43);
N_O3_nr_sm = sgolayfilt(movmean(N_O3_nr,[5,15],1),1,43);
N_O3_fr_sm2 = movmean(movmean(N_O3_fr2,[10,25],1),83);
N_O3_nr_sm2 = movmean(movmean(N_O3_nr2,[10,25],1),83);

maxhid = 1054;
figure
plot(N_O3_fr_sm(:,id),hkm_fr,'b','LineWidth',1.5);hold on
plot(N_O3_nr_sm(1:maxhid,id),hkm_nr(1:maxhid),'r','LineWidth',1.5);
plot(N_O3_fr_sm2(maxhid-186:end,id),hkm_fr(maxhid-186:end),'-.','LineWidth',1.5);
plot(N_O3_nr_sm2(1:maxhid,id),hkm_nr(1:maxhid),'-.','LineWidth',1.52);
set(gca,'FontSize',14)
xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('AD-Far range','AD-Near range','PC-Far range','PC-Near range')
title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18])
%% Calculate the merge ad pc signal
[prof_merge_287 regR287]=adpc_glue_func(prof_an_287,prof_pc_287,5,10,hkm_fr,1);
[prof_merge_299 regR299]=adpc_glue_func(prof_an_299,prof_pc_299,5,10,hkm_fr,1);
[prof_merge_287_nr regR287nr]=adpc_glue_func(prof_an_287_nr,prof_pc_287_nr,3,7,hkm_nr,1);
[prof_merge_299_nr regR299nr]=adpc_glue_func(prof_an_299_nr,prof_pc_299_nr,3,7,hkm_nr,1);

% retrieve the O3 ND using adpc merged signal

frame_len1=31;% ~100 m
frame_len2=53; % ~200m
frame_len3=81;% 303.75m
 
% prof_merge_299(hkm_fr>7,:)= movmean(prof_merge_299(hkm_fr>7,:),53,1);
% prof_merge_287(hkm_fr>7,:)= movmean(prof_merge_287(hkm_fr>7,:),53,1);
% range bin sizes of the different derivative window length
h1_hkm=2; % 0-2km 1:533
h2_hkm=5;% 2-5km 534:1333
[N_O3_fr,ratio_P_fr]=retrieve_o3ND(prof_merge_287,prof_merge_299,...
                                   frame_len1,frame_len2,frame_len3,h1_hkm,h2_hkm,hkm_fr);
[N_O3_nr,ratio_P_nr]=retrieve_o3ND(prof_merge_287_nr,prof_merge_299_nr,...
                                   frame_len1,frame_len2,frame_len3,1.5,h2_hkm,hkm_nr);
% 
id=40;
NO3_nr_ppb=N_O3_nr(:,id)./NDAir_m3_prof*1e9-6.5;
NO3_fr_ppb=N_O3_fr(:,id)./NDAir_m3_prof(end-5801+1:end)*1e9-6.5;
NO3_nr_ppb(hkm_nr>4,:)=nan;
NO3_fr_ppb(hkm_fr<0.9,:)=nan;
NO3_nr_ppb(hkm_nr<4&hkm_nr>2,:)=movmean(NO3_nr_ppb(hkm_nr<4&hkm_nr>2,:),113);
NO3_nr_ppb(hkm_nr<4&hkm_nr>0.25,:)=movmean(movmean(NO3_nr_ppb(hkm_nr<4&hkm_nr>0.25,:),53),63);
NO3_fr_ppb(hkm_fr>7,:)=movmean(NO3_fr_ppb(hkm_fr>7,:),113);
NO3_fr_ppb(hkm_fr>5,:)=movmean(NO3_fr_ppb(hkm_fr>5,:),83);
NO3_fr_ppb(hkm_fr>0.9,:)=movmean(movmean(NO3_fr_ppb(hkm_fr>0.9,:),53),43);


figure 

plot(NO3_fr_ppb,hkm_fr,'LineWidth',2);hold on
plot(NO3_nr_ppb,hkm_nr,'LineWidth',2);
set(gca,'LineWidth',1,'FontSize',14);
xlabel('Ozone mixing ratio (ppb)');ylabel('Altitude (km)');legend('Far range','Near range')
title(['Ozone mixing ratio (ppb) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,150]);ylim([0.2,10]);grid minor                     
%% merge the near range and far range ozone number density
start_merge_hkm =0.8;
end_merge_hkm =0.9;
N_O3_merge = merge_nr_fr_o3ND(N_O3_nr,N_O3_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr);

% the aerosol result above 5.5 are not reliable 
N_O3_bsc(hkm_nr>5.5,:)= 0;
D_aext(hkm_nr>5.5,:)= 0;
ND_O3_corr = (N_O3_merge - N_O3_bsc - D_aext - D_molex);
O3ppb = ND_O3_corr./NDAir_m3*1e9;

ND_O3_corr(hkm_nr<2&hkm_nr>0.25,:)=movmean(movmean(ND_O3_corr(hkm_nr<2&hkm_nr>0.25,:),53),63,1);
ND_O3_corr(hkm_fr>7,:)=movmean(ND_O3_corr(hkm_fr>7,:),113,1);
ND_O3_corr(hkm_fr>5,:)=movmean(ND_O3_corr(hkm_fr>5,:),83,1);
ND_O3_corr_sm(hkm_fr>0.9,:)=movmean(ND_O3_corr(hkm_fr>0.9,:),43,1);

O3ppb_sm = ND_O3_corr_sm./NDAir_m3*1e9;
% O3ppb_sm(hkm_nr>8.6,TimeInHour_avg<15)=nan;
title_str = [datestr(DateTime_avg(1),'yyyy/mm/dd'),' CCNY-DIAL ozone mixing ratio (ppb)'];
NDplot(O3ppb_sm,TimeInHour_avg,hkm_nr,title_str,[0.25,10],[0,120],14);

