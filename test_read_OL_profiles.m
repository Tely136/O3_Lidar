clear all; close all;
%% Load the data files, get time average, dead-time correction
% and background subtraction 
folder_path='C:\Users\OzoneLidar\Documents\Lidar Data\20220415\txt';
save_path='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\ozone lidar data\ozone lidar results\';
nbin=6000;
dzm=3.75;
bgbins=100;
td=280;
nAvg=10;
% [OLfileName,save_path]=read_OL_profiles(folder_path,save_path,nbin,dzm,bgbins,td,nAvg);
% load(OLfileName)
[OLfileName,save_path]=read_OL_profiles2(folder_path,save_path,nbin,dzm,bgbins,td,nAvg);
load(OLfileName)
%% Cloud screening 
% take derivative of the signal profile from the start bin
% the clouds are the place:1. the signal at cloud base is larger than the
% signal below the cloud 2. negative derivative < Threshold
% smooth the profiles

start_bin_fr=160;
sgwin_len_cld=31;% must be odd
cld_screen_prof=nan(size(sigprof.an299));
cld_screen_prof(start_bin_fr:end,:)=sgolayfilt(sigprof.an299(start_bin_fr:end,:),1,sgwin_len_cld);

cld_start_bin = 267;% 1km
cld_end_bin = 5333;% 10km
[cldFlag, cldBaseZ, cldCenterZ,cldBaseZ_ind, cldCenterZ_ind,cld_mask,pz2,d_Pz2]=cld_detect(cld_screen_prof,cld_start_bin,cld_end_bin,hkm);
%% cloud screen test
% figure
% I=imagesc(TimeInHour_avg,hkm,pz2); hold on
% plot(TimeInHour_avg,cldBaseZ(1,:),'o');
% set(gca,'YDir','normal','FontSize',14,'ColorScale','log');colormap('jet');colorbar
% xlabel('Local Time (hour)');ylabel('Altitude (km)')
% title(['Ozone lidar smoothed 299nm Pz2(a.u.) at ',datestr(DateTime_avg(1),'yy/mm/dd')]);
% figure
% I=imagesc(TimeInHour_avg,hkm(1:end-1),d_Pz2);
% set(gca,'YDir','normal','FontSize',14,'ColorScale','log');colormap('jet');colorbar
% xlabel('Local Time (hour)');ylabel('Altitude (km)')
% title(['Ozone lidar smoothed 287nm P derivative (a.u.) at ',datestr(DateTime_avg(1),'yy/mm/dd')]);
% id=13;figure
% plot(pz2(:,id),hkm,'b');hold on
% plot(pz2(cldBaseZ_ind(1,id),id),cldBaseZ(1,id),'o');
% set(gca,'XScale','log');legend('299nm Far','Cloud Base');ylim([0,20]);
% xlabel('Range corrected signal (a.u.)');ylabel('Altitude (km)');
% title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime_avg(id))])

%% cloud screened signal
prof_an_287 = sigprof.an287; prof_an_287(cld_mask)=nan;
prof_an_299 = sigprof.an299; prof_an_299(cld_mask)=nan;
prof_an_287_nr = sigprof.an287nr; prof_an_287_nr(cld_mask)=nan;
prof_an_299_nr = sigprof.an299nr; prof_an_299_nr(cld_mask)=nan;

%% vertical smoothing
start_bin_fr=200;% the start bins when the far range gate opens
prof_an_287(start_bin_fr:end,:) = sgolayfilt(sgolayfilt(prof_an_287(start_bin_fr:end,:),1,31),1,31);
prof_an_299(start_bin_fr:end,:) = sgolayfilt(sgolayfilt(prof_an_299(start_bin_fr:end,:),1,31),1,31);

start_bin_nr=20; % the start bins when the near range gate opens
prof_an_287_nr(start_bin_nr:end,:) = sgolayfilt(sgolayfilt(prof_an_287_nr(start_bin_nr:end,:),1,31),1,31);
prof_an_299_nr(start_bin_nr:end,:) = sgolayfilt(sgolayfilt(prof_an_299_nr(start_bin_nr:end,:),1,31),1,31);

%% ploting the cloud screened smoothed profiles
% id=3;
% figure
% plot(prof_an_287(:,id),hkm,'r');hold on
% plot(prof_an_299(:,id),hkm,'b');
% plot(prof_an_287_nr(:,id),hkm,'g');
% plot(prof_an_299_nr(:,id),hkm,'k');
% set(gca,'XScale','log')
% legend('287nm Far','299nm Far','287nm Near','299nm Near')
% xlabel('Signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar signal at ',datestr(DateTime_avg(id))])
% ylim([0,3])
% grid on
%% calculating the ozone number density from the pon and poff
% SG filter derivative filter window length
frame_len1=31;% ~100 m
frame_len2=53; % ~200m
frame_len3=81;% 303.75m
% window length of smoothing the signal ratio
sm_len_fr=31;
sm_len_nr=15;
% range bin sizes of the different derivative window length
h1=533; % 0-2km 1:533
h2=1333;% 2-5km 534:1333
[N_O3_fr,ratio_P_fr]=retrieve_o3ND(prof_an_287,prof_an_299,start_bin_fr,...
                                   frame_len1,frame_len2,frame_len3,h1,h2,sm_len_fr);
[N_O3_nr,ratio_P_nr]=retrieve_o3ND(prof_an_287_nr,prof_an_299_nr,start_bin_nr,...
                                   frame_len1,frame_len2,frame_len3,h1,h2,sm_len_nr);
%% Plot the ozone number density profile for near range and far range
 id=4;
figure
plot(ratio_P_fr(:,id),hkm,'LineWidth',1.2);hold on;
plot(ratio_P_nr(:,id),hkm,'LineWidth',1.2);
set(gca,'XScale','log');xlabel('ratio_P');ylabel('Altitude (km)');legend('Far range','Near range')
title(['ratio_P ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);grid on

figure
plot(N_O3_fr(:,id),hkm,'LineWidth',1.2);hold on
plot(N_O3_nr(:,id),hkm,'LineWidth',1.2);
xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('Far range','Near range')
title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18])
%% merge the near range and far range ozone number density
start_merge_hkm =0.9;
end_merge_hkm =1;
N_O3_merge = merge_nr_fr_o3ND(N_O3_nr,N_O3_fr,start_merge_hkm,end_merge_hkm,hkm);

figure
plot(N_O3_merge(:,id),hkm,'LineWidth',1.2);hold on
xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('merged')
title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18])

%% Vertical smoothing the ozone number density
h1=533;h2=1334;
movnum1=53;% 200m 0-2km 1:533
movnum2=81;% 300m  2-5km 534:1333
movnum3=133;% 500m  >5km  1334:end 
N_O3_merge_sm = sglp_smooth(N_O3_merge,start_bin_fr,h1,h2,movnum1,movnum2,movnum3);

figure
plot(N_O3_merge_sm(:,id),hkm,'LineWidth',1.2);
xlabel('Ozone number density (molecule / m^3)');ylabel('Altitude (km)');legend('merge')
title(['Ozone number density (molecule / m^3) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
grid on;xlim([0,2.5e18])
%% calculate air density
% option 1: using International Standard Atmosphere (ISA)
% option 2: using micro-ratiometer data (10-min resolution/daily average)
% option 3: using radiosonde 


% option 1: Standard T and P profile,
% p0 = 1013.25*1e2; % surface pressure 1013.25 hpa = 101325pa
% T0 = 288.15;% surface temp (K)
[NDAir_m3,D_molex] = isa_NDair_m3(hkm,TimeInHour_avg);

% option 2: using micro-ratiometer data (10-min resolution/daily average)
mwrfile='~\mwrdata\2022-04-15_00-04-06_lv2';
tdiff = 4; % tdiff is the time difference between the local timezone and utc
[NDAir_m3,T,P] = mwr_NDair_m3(mwrfile,hkm,TimeInHour_avg,tdiff);
