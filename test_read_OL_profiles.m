clear all; close all;
%% Load the data files, get time average, dead-time correction
% and background subtraction 
folder_path='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar/20220531/';
save_path='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar_results/20220531/';
nbin=8000;
dzm=3.75;
bgbins=100;
td=280;
nAvg=10;
% [OLfileName,save_path]=read_OL_profiles(folder_path,save_path,nbin,dzm,bgbins,td,nAvg);

[OLfileName,save_path]=read_OL_profiles2(folder_path,save_path,nbin,dzm,bgbins,td,nAvg);
load(OLfileName)
% OLfileName='ol220429_1113_1804.mat';
% load([save_path,OLfileName]);
%% Remove the signals which collect befor the gate opens
start_bin_fr=218;
start_bin_nr=14;
[new_sigprof,hkm_fr,hkm_nr] = gate_prof(sigprof,start_bin_fr,start_bin_nr,hkm);
%% Calculate the merge ad pc signal
[new_sigprof.merge287fr regR287]=adpc_glue_func(new_sigprof.an287,new_sigprof.pc287,5,10,hkm_fr,1);
[new_sigprof.merge299fr regR299]=adpc_glue_func(new_sigprof.an299,new_sigprof.pc299,5,10,hkm_fr,1);
[new_sigprof.merge287nr regR287nr]=adpc_glue_func(new_sigprof.an287nr,new_sigprof.pc287nr,3,7,hkm_nr,1);
[new_sigprof.merge299nr regR299nr]=adpc_glue_func(new_sigprof.an299nr,new_sigprof.pc299nr,3,7,hkm_nr,1);

%% Cloud screening 
% take derivative of the signal profile from the start bin
% the clouds are the place:1. the signal at cloud base is larger than the
% signal below the cloud 2. negative derivative < Threshold
% smooth the profiles

sgwin_len_cld=31;% must be odd
cld_screen_prof=sgolayfilt(new_sigprof.an299,1,sgwin_len_cld);

cld_start_bin = 1;% 1km
cld_end_bin = 10;% 10km
[cldFlag, cldBaseZ, cldCenterZ,cldBaseZ_ind, cldCenterZ_ind,cld_mask,pz2,d_Pz2]=cld_detect(cld_screen_prof,cld_start_bin,cld_end_bin,hkm_fr);
%% cloud screen test
id=30;
plot_cld_screen_result(DateTime_avg,TimeInHour_avg,hkm_fr,pz2,d_Pz2,cldBaseZ,cldBaseZ_ind,id)

%% cloud screened signal
% Use analog signal
prof_an_287 = new_sigprof.an287; %prof_an_287(cld_mask)=nan;
prof_an_299 = new_sigprof.an299; %prof_an_299(cld_mask)=nan;
prof_an_287_nr = new_sigprof.an287nr; 
prof_an_299_nr = new_sigprof.an299nr; 

% Use PC signal
prof_pc_287 = new_sigprof.pc287; %prof_pc_287(cld_mask)=nan;
prof_pc_299 = new_sigprof.pc299; %prof_pc_299(cld_mask)=nan;
prof_pc_287_nr = new_sigprof.pc287nr; 
prof_pc_299_nr = new_sigprof.pc299nr; 

prof_merge_287_fr=new_sigprof.merge287fr; %prof_merge_287(cld_mask)=nan;
prof_merge_299_fr=new_sigprof.merge299fr; %prof_merge_299(cld_mask)=nan;
prof_merge_287_nr=new_sigprof.merge287nr; 
prof_merge_299_nr=new_sigprof.merge299nr; 

start_merge_hkm = 0.9;
end_merge_hkm =1;
z1 = 0.9;
z2 =1.2;
prof_merge_287 = merge_nr_fr_prof(prof_merge_287_nr,prof_merge_287_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr,z1,z2);
prof_merge_299 = merge_nr_fr_prof(prof_merge_299_nr,prof_merge_299_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr,z1,z2);
%% Calculate the SNR
bg299 =repmat(sigprof.bgpc299,length(hkm_nr),1);
bg287 = repmat(sigprof.bgpc287,length(hkm_nr),1);
snr_299 = sqrt(12000/40*(prof_merge_299).^2./(prof_merge_299+bg299));
snr_287 = sqrt(12000/40*(prof_merge_287).^2./(prof_merge_287+bg287));

id=20;
figure
plot(snr_287(:,id),hkm_nr,'g');hold on
plot(snr_299(:,id),hkm_nr,'k');
set(gca,'XScale','log')
legend('287 nm','299 nm')
xlabel('SNR');
ylabel('Altitude (km)');
title(['SNR ',datestr(DateTime_avg(id))])
ylim([0,15])
grid on
%
%% ploting the cloud screened smoothed profiles
id=20;
y_lim =[0,10];
tit_str='CCNY-O3-DIAL AD-PC Merged Signal ';
xlb='Merged Signal (MHz)';
le={'287nm Far','299nm Far','287nm Merge','299nm Merge'};
plot_fr_nr_Prof(prof_merge_287_fr,prof_merge_299_fr,prof_merge_287,prof_merge_299,...
                hkm_fr,hkm_nr,DateTime_avg,id,tit_str,xlb,y_lim,le);
%% calculating the ozone number density from the pon and poff
% SG filter derivative filter window length
frame_len1=31;% ~100 m
frame_len2=53; % ~200m
frame_len3=81;% 303.75m
 
% range bin sizes of the different derivative window length
h1_hkm=2; % 0-2km 1:533
h2_hkm=5;% 2-5km 534:1333

[N_O3,ratio_P]=retrieve_o3ND(prof_merge_287,prof_merge_299,...
                                   frame_len1,frame_len2,frame_len3,h1_hkm,h2_hkm,hkm_nr);
[N_O3_nr,ratio_P_nr]=retrieve_o3ND(prof_an_287_nr,prof_an_299_nr,...
                                   frame_len1,frame_len2,frame_len3,1.5,h2_hkm,hkm_nr);
[N_O3_fr,ratio_P_fr]=retrieve_o3ND(prof_merge_287_fr,prof_merge_299_fr,...
                                   frame_len1,frame_len2,frame_len3,1.5,h2_hkm,hkm_fr);
                               
%% Plot the ozone number density profile for near range and far range
id=4;
% figure
% plot(ratio_P_fr(:,id),hkm_fr,'LineWidth',1.2);hold on;
% plot(ratio_P_nr(:,id),hkm_nr,'LineWidth',1.2);
% set(gca,'XScale','log');xlabel('ratio_P');ylabel('Altitude (km)');legend('Far range','Near range')
% title(['ratio_P ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);grid on

figure
plot(N_O3_fr(:,id),hkm_fr,'LineWidth',1.2);hold on
plot(N_O3_nr(:,id),hkm_nr,'LineWidth',1.2);
plot(N_O3(:,id),hkm_nr,'LineWidth',1.2);
xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('Far','Near','Merge')
title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18])

%% Vertical smoothing of the ozone number density
% N_O3_fr_sm = sgolayfilt(movmean(N_O3_fr,[5,15],1),1,43);
% N_O3_nr_sm = sgolayfilt(movmean(N_O3_nr,[5,15],1),1,43);
% id=13;
% figure
% plot(N_O3_fr_sm(:,id),hkm_fr,'LineWidth',1.2);hold on
% plot(N_O3_nr_sm(:,id),hkm_nr,'LineWidth',1.2);
% xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('Far range','Near range')
% title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
% grid on;xlim([0,2.5e18])

%% merge the near range and far range ozone number density
start_merge_hkm =1;
end_merge_hkm =1.1;
N_O3_merge = merge_nr_fr_o3ND(N_O3_nr,N_O3_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr);

figure
plot(N_O3_merge(:,id),hkm_nr,'LineWidth',1.2);hold on
xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('merged')
title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18])

 %% Vertical smoothing the ozone number density
h1_hkm=find(hkm_nr>2,1,'first'); % index at 2km
h2_hkm=find(hkm_nr>5,1,'first'); % index at 5km

movnum1=53;% 200m 0-2km 1:533
movnum2=81;% 300m  2-5km 534:1333
movnum3=133;% 500m  >5km  1334:end 
N_O3_merge_sm = sglp_smooth(N_O3_merge,1,h1_hkm,h2_hkm,movnum1,movnum2,movnum3);

% figure
% plot(N_O3_merge_sm(:,id),hkm_nr,'LineWidth',1.2);
% xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('merge')
% title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
% grid on;xlim([0,2.5e18])
%% calculate air density
% option 1: using International Standard Atmosphere (ISA)
% option 2: using micro-ratiometer data (10-min resolution/daily average)
% option 3: using radiosonde 

sigmaOn=203.4*10^(-20);% O3 cross section at 287.2 (cm^2/molecule)
sigmaOff=45.51*10^(-20);% O3 cross section at 299.1 (cm^2/molecule)
d_sigma=(sigmaOn-sigmaOff)*(10^-2)^2; % delta cross section (m^2/molecule)

%% option 1: Standard T and P profile,
p0 = 1013.25*1e2; % surface pressure 1013.25 hpa = 101325pa
T0 = 288.15;% surface temp (K)
[NDAir_m3_prof,D_molex_prof,am287,am299] = isa_NDair_m3(hkm_nr,p0,T0,d_sigma);
len_t=length(TimeInHour_avg);
NDAir_m3= repmat(NDAir_m3_prof, 1,len_t);
D_molex = repmat(D_molex_prof, 1,len_t);

%% option 2: using radiosonde data
sondefile='/Users/Tinker/Documents/MATLAB/ozonelidar/sondeData/20220429_OKX_z12.txt';
[NDAir_m3_prof,D_molex_prof,am287,am299] = sonde_NDair_m3(sondefile,hkm_nr,d_sigma);
len_t=length(TimeInHour_avg);
NDAir_m3 = repmat(NDAir_m3_prof, 1,len_t);
D_molex = repmat(D_molex_prof, 1,len_t);
%% Make a nice plote of the fr and nr ozone ppb 
% id=13;
% NO3_nr_ppb=N_O3_nr(:,id)./NDAir_m3_prof*1e9;
% NO3_fr_ppb=N_O3_fr(:,id)./NDAir_m3_prof(end-5801+1:end)*1e9;
% NO3_nr_ppb(hkm_nr>2,:)=nan;
% NO3_fr_ppb(hkm_fr<0.9,:)=nan;
% NO3_nr_ppb(hkm_nr<2&hkm_nr>0.25,:)=movmean(movmean(NO3_nr_ppb(hkm_nr<2&hkm_nr>0.25,:),53),43);
% NO3_fr_ppb(hkm_fr>0.9,:)=movmean(movmean(NO3_fr_ppb(hkm_fr>0.9,:),53),43);
% 
% 
% figure
% plot(NO3_fr_ppb,hkm_fr,'LineWidth',2);hold on
% plot(NO3_nr_ppb,hkm_nr,'LineWidth',2);
% set(gca,'LineWidth',1,'FontSize',14);
% xlabel('Ozone mixing ratio (ppb)');ylabel('Altitude (km)');legend('Far range','Near range')
% title(['Ozone mixing ratio (ppb) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
% grid on;xlim([0,100]);ylim([0.2,7.5])

%% option 3: using micro-ratiometer data (10-min resolution/daily average)
mwrfile='~\mwrdata\2022-04-15_00-04-06_lv2';
tdiff = 4; % tdiff is the time difference between the local timezone and utc
[NDAir_m3,T,P] = mwr_NDair_m3(mwrfile,hkm,TimeInHour_avg,tdiff);

%% Aerosol backscatter retrieval and correction

% [absc_299_retri,ND_O3,N_O3_bsc,D_aext]=OL299_aero_retrieval(N_O3,NDAir_m3_prof,prof_merge_299,hkm_nr)

% option2 using ceilometer data 
aodfile='/Users/Tinker/Documents/MATLAB/ozonelidar/AOD/20220301_20220602_CCNY.txt';
timediff = 4; % time difference between the local time and utc time (EDT = 4, EST =5)
[chm15k,N_O3_bsc,D_aext]=aero_ext_bsc_retrieval(NDAir_m3_prof,d_sigma,hkm_nr,TimeInHour_avg,timediff,aodfile);

h1_hkm=find(hkm_nr>2,1,'first'); % index at 2km
h2_hkm=find(hkm_nr>5,1,'first'); % index at 5km

movnum1=53;% 200m 0-2km 1:533
movnum2=81;% 300m  2-5km 534:1333
movnum3=133;% 500m  >5km  1334:end 
N_O3_bsc_sm = sglp_smooth(N_O3_bsc,1,h1_hkm,h2_hkm,movnum1,movnum2,movnum3);

figure
plot(N_O3_bsc_sm(:,id),hkm_nr,'LineWidth',1.2);
xlabel('Aerosol correction terms (molecule / m^3)');ylabel('Altitude (km)');legend('merge')
title(['Aerosol correction terms  (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18]);
% 

%% Plot the ozone number density or ppb
% ozone number density, no correction
y_limit = [0.27,8]; %km
z_limit = [0,2e18];
fontsize= 13;
chm_t=chm15k.time_aeroprof - timediff;
ind_t = chm_t<TimeInHour_avg(end)&chm_t>TimeInHour_avg(1); 
% ozone ND, no correction
title_str = ['Ozone number density (m^{-3}),no correction ',datestr(DateTime_avg(1),'yy/mm/dd')];
NDplot(N_O3_merge_sm,TimeInHour_avg,hkm_nr,title_str,y_limit,z_limit,fontsize);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');legend('PBLH')

% ozone ppb, no correction
title_str = ['Ozone mixing ratio (ppb),no correction ',datestr(DateTime_avg(1),'yy/mm/dd')];
NDplot(N_O3_merge_sm./NDAir_m3*1e9,TimeInHour_avg,hkm_nr,title_str,y_limit,[0,100],fontsize);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');legend('PBLH')


% molecular extinction correction term
title_str = ['Molecular correction (ext) terms (ppb)',datestr(DateTime_avg(1),'yy/mm/dd')];
NDplot(D_molex./NDAir_m3*1e9,TimeInHour_avg,hkm_nr,title_str,y_limit,[0,10],fontsize);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');legend('PBLH')

% aerosol extinction correction term
title_str = ['Aerosol correction (ext) terms (ppb)',datestr(DateTime_avg(1),'yy/mm/dd')];
NDplot(D_aext./NDAir_m3*1e9,TimeInHour_avg,hkm_nr,title_str,y_limit,[0,10],fontsize);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');legend('PBLH')

% aerosol bsc correction term
title_str = ['Aerosol correction (bsc) terms (ppb)',datestr(DateTime_avg(1),'yy/mm/dd')];
NDplot(N_O3_bsc_sm./NDAir_m3*1e9,TimeInHour_avg,hkm_nr,title_str,y_limit,[-10,10],fontsize);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');legend('PBLH')

% ozone ppb, with molecular , aerosol correction
% the aerosol result above 5.5 are not reliable 
N_O3_bsc(hkm_nr>5.5,:)= 0;
D_aext(hkm_nr>5.5,:)= 0;
ND_O3_corr = (N_O3_merge - N_O3_bsc - D_aext- D_molex);
O3ppb = ND_O3_corr./NDAir_m3*1e9;
movnum1=53;% 200m 0-2km 1:533
movnum2=81;% 300m  2-5km 534:1333
movnum3=133;% 500m  >5km  1334:end 
ND_O3_corr_sm = sglp_smooth(ND_O3_corr,1,h1_hkm,h2_hkm,movnum1,movnum2,movnum3);
O3ppb_sm = ND_O3_corr_sm./NDAir_m3*1e9;

title_str = [datestr(DateTime_avg(1),'yyyy/mm/dd'),' CCNY-DIAL ozone mixing ratio (ppb)'];
NDplot(O3ppb_sm,TimeInHour_avg,hkm_nr,title_str,[0.27,6],[0,100],14);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');
legend('PBLH')
%% plot the aerosol backscatter at the same time as ozone dial
chm_t=chm15k.time_aeroprof - timediff;
chm_h=chm15k.height;
ind_t = chm_t<TimeInHour_avg(end)&chm_t>TimeInHour_avg(1);
title_str = [datestr(DateTime_avg(1),'yyyy/mm/dd'),' CCNY-Ceilometer aerosol backscatter (km^{-1}sr^{-1})'];
NDplot(movmean(chm15k.aero_bsa1064(:,ind_t),3,2),chm_t(ind_t),chm_h,title_str,[0.25,6],[0,2e-3],12);
hold on
plot(chm_t(ind_t),chm15k.pbl_aero(ind_t),'*');
legend('PBLH')
%% Read in the ground ozone file and compare with the ozone dial value at 300m
gd_ozone_file='/Users/Tinker/Documents/MATLAB/ozonelidar/gd_ozone_measur/20220521_gdO3.mat';
load(gd_ozone_file);
ind_gd_ozone=isbetween(ground_O3_ppb_datetime,DateTime_avg(1)-minutes(10),DateTime_avg(end)+minutes(10));
% ind_sh_ozone=isbetween(shed_O3_ppb_datetime,DateTime_avg(1)-minutes(10),DateTime_avg(end)+minutes(10));

indh=find(hkm_nr>0.26,1,'first');
figure;
% plot(shed_O3_ppb_datetime(ind_sh_ozone),shed_O3_ppb(ind_sh_ozone),'r-');hold on;
plot(ground_O3_ppb_datetime(ind_gd_ozone),ground_O3_ppb(ind_gd_ozone),'r-');hold on;
plot(DateTime_avg,O3ppb(indh,:),'bo-');
set(gca,'LineWidth',1,'FontSize',14);
legend('in-situ','O_3-DIAL measurement at 250m');
xlabel('EDT Time (Hour)');
ylabel('Ozone Mixing Ratio (ppb)');
ylim([20,120]);
title(['DIAL and Ground Ozone Measurements'])
grid on;

