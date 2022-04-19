clear all; close all;
% Get a list of all txt files in the current folder, or subfolders of it.
folderPath='C:\Users\OzoneLidar\Documents\Lidar Data\20211101_HV_PCtest\txt';
mwrfilename='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\mwrdata\2021-11-01_00-04-06_lv2.csv';
aodfile='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\AeronetAoddata\20210101_20211231_CCNY.txt';


fds = fileDatastore(folderPath, 'ReadFcn', @importdata,'FileExtensions',{'.txt'});
fullFileNames = fds.Files;
numFiles = length(fullFileNames);

% Loop over all files reading them in and plotting them.

% te2192013.483376
% CCNY     20/09/2021 13:47:34 20/09/2021 13:48:33 0100 -073.950000 0040.821000 00.0 00.0
% 0001201 0020 0000000 0020 04 0000000 0020 0000000 0000
% 1 0 1 04000 1 0000 3.75 00287.o 0 0 01 000 16 001201 0.500 BT0
% 1 1 1 04000 1 0000 3.75 00287.o 0 0 00 000 00 001201 3.1746 BC0
% 1 0 1 04000 1 0000 3.75 00299.o 0 0 01 000 16 001201 0.500 BT1
% 1 1 1 04000 1 0000 3.75 00299.o 0 0 00 000 00 001201 3.1746 BC1
% %
nbin=6000;
profile_287_an1=nan(nbin,numFiles);
profile_299_an1=nan(nbin,numFiles);
profile_287_an_raw=nan(nbin,numFiles);
profile_299_an_raw=nan(nbin,numFiles);
height=3.75*[1:1:nbin]';
DateTime=NaT(1,numFiles);
TimeInHour=nan(1,numFiles);
for k = 1 : numFiles
    %fprintf('Now reading file %s\n', fullFileNames{k});
    % read both data and metadata from the filestore
    Data=read(fds);
    profile_287_an_raw(:,k)= Data.data(:,1);
    profile_299_an_raw(:,k)= Data.data(:,3);
    profile_287_an1(:,k)= profile_287_an_raw(:,k)-median(profile_287_an_raw(nbin-500:nbin,k));
    profile_299_an1(:,k)= profile_299_an_raw(:,k)-median(profile_299_an_raw(nbin-500:nbin,k));
    temp=cell2mat(Data.textdata(2,1));
    infmt='dd/MM/yyyy HH:mm:ss';
    DateTime(k)=datetime(temp(10:28),'InputFormat',infmt);
    TimeInHour(k)=hour(DateTime(k))+minute(DateTime(k))/60+second(DateTime(k))/3600;
end
profile_287_an=movmedian(profile_287_an1,21,1);
profile_299_an=movmedian(profile_299_an1,21,1);
%% Plot the profile
% id=3;
% figure
% plot(profile_287_an(:,id),height/1000,'r');hold on;
% plot(profile_299_an(:,id),height/1000,'b');
% legend('287nm','299nm')
% xlabel('Signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar signal profile at ',datestr(DateTime(id))])
% ylim([0.05,5])

%% Plot the signal color image
figure
subplot(1,2,1)
I=imagesc(TimeInHour,height./1000,profile_287_an);
set(gca,'YDir','normal');
set(gca,'FontSize',14);
set(gca,'ColorScale','log');
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone lidar 287nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

subplot(1,2,2)
I=imagesc(TimeInHour,height./1000,profile_299_an);
set(gca,'YDir','normal');
set(gca,'ColorScale','log');
set(gca,'FontSize',14);
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone lidar 299nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

% 
%% Plot the range corrected
[m,n]=size(profile_287_an_raw);
pz2_287=nan(size(profile_287_an));
pz2_299=nan(size(profile_299_an));
for i=1:n
   pz2_287(:,i)= profile_287_an(:,i).*(height./1000).^2;
   pz2_299(:,i)= profile_299_an(:,i).*(height./1000).^2;
end
% 
% id=48;
% figure
% plot(pz2_287(:,id),height/1000,'r');hold on;
% plot(pz2_299(:,id),height/1000,'b');
% legend('287nm','299nm')
% xlabel('Range corrected signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime(id))])
% ylim([0.05,15])
% 
% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour,height./1000,pz2_287);
% set(gca,'YDir','normal');
% set(gca,'FontSize',14);
% set(gca,'ColorScale','log');
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar 287nm P\cdotz^2 (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
% 
% subplot(1,2,2)
% I=imagesc(TimeInHour,height./1000,pz2_299);
% set(gca,'YDir','normal');
% set(gca,'ColorScale','log');
% set(gca,'FontSize',14);
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar 299nm P\cdotz^2 (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

%% Time averaging of the lidar signal
nAvg=5;% 5min average

N=floor(n/nAvg);
profile_287_avgT=nan(m,N);
profile_299_avgT=nan(m,N);
profile_287_avg=nan(m,N);
profile_299_avg=nan(m,N);
TimeInHour_avg=nan(1,N);
DateTime_avg=NaT(1,N);

for i=1:N
    ind_t=(i-1)*nAvg;
    profile_287_avgT(:,i)=mean(profile_287_an(:,1+ind_t:nAvg+ind_t),2,'omitnan');
    profile_299_avgT(:,i)=mean(profile_299_an(:,1+ind_t:nAvg+ind_t),2,'omitnan');
    TimeInHour_avg(i)=mean(TimeInHour(1+ind_t:nAvg+ind_t));
    DateTime_avg(i)=DateTime(5+ind_t);
end 
%% Vertical filtering of the lidar signal using 1st-order Savitzky-Golay (boxcar smoothing)
% order = 1;% 
% framelen = 81;% smooth bins
% 
% profile_287_avg=sgolayfilt(profile_287_avgT,order,framelen);
% profile_299_avg=sgolayfilt(profile_299_avgT,order,framelen); 
% 0-1km  78.5m runing mean --> noise filtering
temp=movmean(profile_287_avgT,21,1);
profile_287_avg(1:266,:)=temp(1:266,:);
% 1-2.5 km 198.75m running mean
temp=movmean(profile_287_avgT,53,1);
profile_287_avg(267:666,:)=temp(267:666,:);
% 2.5-5 km 300.75m running mean
temp=movmean(profile_287_avgT,81,1);
profile_287_avg(667:1334,:)=temp(667:1334,:);
% 5-22.5 km 498.75m running mean
temp=movmean(profile_287_avgT,133,1);
profile_287_avg(1335:end,:)=temp(1335:end,:);

temp=movmean(profile_299_avgT,21,1);
profile_299_avg(1:266,:)=temp(1:266,:);
% 1-2.5 km 198.75m running mean
temp=movmean(profile_299_avgT,53,1);
profile_299_avg(267:666,:)=temp(267:666,:);
% 2.5-5 km 300.75m running mean
temp=movmean(profile_299_avgT,81,1);
profile_299_avg(667:1334,:)=temp(667:1334,:);
% 5-22.5 km 498.75m running mean
temp=movmean(profile_299_avgT,133,1);
profile_299_avg(1335:end,:)=temp(1335:end,:);


figure
subplot(1,2,1)
I=imagesc(TimeInHour_avg,height./1000,profile_287_avg);
set(gca,'ColorScale','log');
set(gca,'YDir','normal','FontSize',14);
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Smoothed Ozone lidar 287nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

subplot(1,2,2)
I=imagesc(TimeInHour_avg,height./1000,profile_299_avg);
set(gca,'ColorScale','log');
set(gca,'YDir','normal','FontSize',14);
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Smoothed Ozone lidar 299nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
%% Plot the smoothed range corrected signal
pz2_287_avg=nan(size(profile_287_avg));
pz2_299_avg=nan(size(profile_299_avg));
for i=1:N
   pz2_287_avg(:,i)= profile_287_avg(:,i).*(height./1000).^2;
   pz2_299_avg(:,i)= profile_299_avg(:,i).*(height./1000).^2;
end


id=22;
figure
semilogx(pz2_287_avg(:,id),height/1000,'r');hold on;
semilogx(pz2_299_avg(:,id),height/1000,'b');
legend('287nm','299nm')
xlabel('Range corrected signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar smoothed P\cdotz^2 at ',datestr(DateTime(id))])
ylim([0.05,15])

% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour_avg,height./1000,pz2_287_avg,[0,8]);
% set(gca,'YDir','normal');
% set(gca,'FontSize',14);
% %set(gca,'ColorScale','log');
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar smoothed 287nm P\cdotz^2 (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
% 
% subplot(1,2,2)
% I=imagesc(TimeInHour_avg,height./1000,pz2_299_avg,[0,8]);
% set(gca,'YDir','normal');
% %set(gca,'ColorScale','log');
% set(gca,'FontSize',14);
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar smoothed 299nm P\cdotz^2 (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
% 
%% Calculate ozone number density using equation
% the cross section is given by (Molina and Molina 1986)
%N=1/(2*O3CrossSec*dr)ln(Poff(r+dr)/Poff*(Pon(r)/Pon(r+dr)))
sigmaOn=203.4*10^(-20);% O3 cross section at 287 (cm^2/molecule)
sigmaOff=45.51*10^(-20);% O3 cross section at 299 (cm^2/molecule)
d_sigma=(sigmaOn-sigmaOff)*(10^-2)^2; % delta cross section (m^2/molecule)
dr=1;%
dR=dr*3.75; % range increment (m)
Poff=profile_299_avg;
Pon=profile_287_avg;
N_O3=nan(m-1,N);
for i=1:m-1
    diff=(Poff(i+1,:)./Poff(i,:)).*(Pon(i,:)./Pon(i+1,:));
    diff(diff<1)=nan; % set all negative derivative to nan
    N_O3(i,:)=1/(2*d_sigma*dR).*log(diff);
    
end
% Using covolution with the SG filter to calculate the first order
% derivative 
[b,g] = sgolay(2,21);
ratio_P=Pon./Poff;
ratio_P(ratio_P<0)=nan;
Ln_ratio_P=log(ratio_P);
dR=3.75;
N_O3_2=nan(size(ratio_P));
for i=1:N
  diff= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g(:,2), 'same');
  N_O3_2(:,i)=1/(2*d_sigma).*diff;
end
figure
I=imagesc(TimeInHour_avg,height./1000,N_O3,[1e17,1e20]);
set(gca,'YDir','normal','FontSize',14);
set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density using noraml difference (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,5])

figure
I=imagesc(TimeInHour_avg,height./1000,N_O3_2,[1e17,1e20]);
set(gca,'YDir','normal','FontSize',14);
set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_2))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density using SG filter (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,5])

%% molecular extinction correction 
% mwr=loadMWRdata4lidar(mwrfilename,N_O3_2,DateTime_avg,height);
mwr=loadMWRdata(mwrfilename);
chm15kRetrieval=CHM15kForwardRetrieval_new(mwr,aodfile);
%% interpolate the time and vertical of mwr and chm15k
t_chm15k=chm15kRetrieval.time_array_utc;% utc hour
h_chm15k=chm15kRetrieval.height; % in km
t_mwr=mwr.time_array_utc;% utc hour
h_mwr=mwr.height;% km

timediff=hours(4); % time difference between the local time and utc time (EDT = 4; EST =5)

%% calculate the molecular number density and extinction
lamda_1=287.2; %% unit: nm
WAVE_1=lamda_1/(1e+3);  %% wavelength unit: micro, um
lamda_2=299.1; %% unit: nm
WAVE_2=lamda_2/(1e+3);  %% wavelength unit: micro, um
PI=3.1415926; 
ns= (7.247249e+18)* 1013.5/288.15; %% /* unit: molecule/cm3 */
    % 7.247249e+18=Ns*Ys/Ps (S: surface, Ns=2.546899e+19 /cm3, P3=1013.25mbar, Ts=288.15K % 
    
ms1=1.0+(1e-8)*(6432.8 + 2949810/(146- 1.0/WAVE_1/WAVE_1) + 25540/(41-1/WAVE_1/WAVE_1));
sigma1=(6+3*0.035)/(6-7*0.035)*8*PI*PI*PI*(ms1*ms1-1)*(ms1*ms1-1)/3.0/WAVE_1/WAVE_1/WAVE_1/WAVE_1/ns/ns*(1e+16);
               %%% Rayleigh scattering cross section, unit: cm2 
               %%% Corrected by depolarization ratio delta=0.035
molex_1=sigma1* (1e-6*mwr.NDAir_m3)*1e2;   %% unit: m^-1 

delta_molex=(1-(lamda_2/lamda_1)^(-4))* molex_1;
D_molex=delta_molex./d_sigma;

% aerosol total backscatter correction

S1= 8*pi/3;
ae=chm15kRetrieval.ae_1020_340;
absc_299=(1064/299)^ae*chm15kRetrieval.aero_bsa1064;
absc_287=(1064/287)^ae*chm15kRetrieval.aero_bsa1064;
mbsc_299=(1064/299)^4/S1*chm15kRetrieval.am_interp;
mbsc_287=(1064/287)^4/S1*chm15kRetrieval.am_interp;
totbsc_299off=absc_299+mbsc_299;
totbsc_287on=absc_287+mbsc_287;
ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;

[b,g] = sgolay(2,21);
ln_bsc=log(ratio_totbsc_onoff);
dz=15;
N_O3_bsc=nan(size(ln_bsc));
for i=1:length(t_chm15k)
  diff= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g(:,2), 'same');
  N_O3_bsc(:,i)=1/(2*d_sigma).*diff;
end
% aerosol extinction correction term

delta_aext=1e-3*((1064/287.2)^ae-(1064/299.1)^ae)*chm15kRetrieval.aero_ext1064;
D_aext=delta_aext./d_sigma;

% convert to the same time base of ozone lidar and vertical res as CHM15k
[~,N]=size(N_O3_2);
M=length(h_chm15k);

N_O3_interp=nan(M,N);
NDAir_m3_interp=nan(M,N);
D_molex_interp=nan(M,N);
N_O3_bsc_interp=nan(M,N);
D_aext_interp=nan(M,N);
pbl_interp=nan(1,N);


for i=1:N  
    N_O3_interp(:,i)=interp1(height./1000,N_O3_2(:,i),h_chm15k);
    
    indt=isbetween(t_mwr-timediff,DateTime_avg(i)-minutes(5),DateTime_avg(i)+minutes(5));
    temp=mean(mwr.NDAir_m3(:,indt),2,'omitnan')';
    NDAir_m3_interp(:,i)=interp1(h_mwr,temp,h_chm15k);
    
    indt=isbetween(t_mwr-timediff,DateTime_avg(i)-minutes(5),DateTime_avg(i)+minutes(5));
    temp=mean(D_molex(:,indt),2,'omitnan')';
    D_molex_interp(:,i)=interp1(h_mwr,temp,h_chm15k);
    
    indt=isbetween(t_chm15k-timediff,DateTime_avg(i)-minutes(5),DateTime_avg(i)+minutes(5));
    N_O3_bsc_interp(:,i)=mean(N_O3_bsc(:,indt),2,'omitnan')';
    D_aext_interp(:,i)=mean(D_aext(:,indt),2,'omitnan')';
    pbl_interp(1,i)=mean(chm15kRetrieval.pbl_aero(indt),2,'omitnan')';
    
end 
%% convert the number density to mixing ratio by volume (ppbv)
% ppbv_D_molex=D_molex./mwr.NDAir_m3*10^9;
% ppbv_O3=(N_O3_2)./mwr.NDAir_m3*10^9;
% 

% Simple converstion: molecule/m^3 to ppm (parts per million)
% Concentration = Mass of one molecule * Number of molecule
% Mass of one molecule = molecular weight * AMU (atomic mass unit)
% O3 molecular weight= 48 AMU = 48*1.66065*10^-27 kg 
% 1 g O3 / m3 = 467 PPM O3
% 1 PPM O3 = 2.14 mg O3/m3 @ Standard conditions P = 1013.25 MB, T = 273.3 K
% M_O3= 48* 1.66065*10^-27; % mass of one Ozone molecule (kg)
% C_O3=M_O3.*N_O3;% ozone concentration kg/m^3
% ppm_O3=C_O3./(2.14*10^-6); 


%%
figure
subplot(2,2,1)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp))
colormap('jet')
colorbar
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3})',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,h_chm15k,D_molex_interp);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_molex_interp))
colormap('jet')
colorbar
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} correction term ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,h_chm15k,-N_O3_bsc_interp,[-0.5,0.5]*1e18);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_bsc_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \beta correction term ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,h_chm15k,D_aext_interp);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_aext_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{a} correction term',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])
%%
figure
subplot(2,2,1)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv),no correction ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,h_chm15k,(D_molex_interp)./NDAir_m3_interp*1e9,[0,10]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_molex_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} term ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,h_chm15k,(-N_O3_bsc_interp)./NDAir_m3_interp*1e9,[-10,10]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_bsc_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \beta term ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,h_chm15k,(D_aext_interp)./NDAir_m3_interp*1e9,[0,10]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_aext_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{a} term ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

%%
figure
subplot(2,2,1)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp,[1e17,4e18]);
set(gca,'YDir','normal','FontSize',13);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp))
colormap('jet')
colorbar
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}),no correction ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp-D_molex_interp,[1e17,4e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp-D_molex_interp))
colormap('jet')
colorbar
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} corrected ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp-D_molex_interp-N_O3_bsc_interp,[1e17,4e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp-D_molex_interp-N_O3_bsc_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} \beta corrected ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp-D_molex_interp-N_O3_bsc_interp-D_aext_interp,[1e17,4e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp-D_molex_interp-N_O3_bsc_interp-D_aext_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m}, \beta and \alpha_{a} corrected ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

%%
figure
subplot(2,2,1)
I=imagesc(TimeInHour_avg,h_chm15k,N_O3_interp./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv),no correction ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,h_chm15k,(N_O3_interp-D_molex_interp)./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp-D_molex_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} corrected ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,h_chm15k,(N_O3_interp-D_molex_interp-N_O3_bsc_interp)./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp-D_molex_interp-N_O3_bsc_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} and \beta corrected ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,h_chm15k,(N_O3_interp-D_molex_interp-N_O3_bsc_interp-D_aext_interp)./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_interp-D_molex_interp-N_O3_bsc_interp-D_aext_interp))
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m}, \beta and \alpha_{a} corrected ',datestr(DateTime(1),'yy/mm/dd')],'FontSize',11);
ylim([0.45,5])


O3_retrieval.profile_287_an_raw=profile_287_an_raw;
O3_retrieval.profile_299_an_raw=profile_299_an_raw;
O3_retrieval.profile_287_an=profile_287_an;
O3_retrieval.profile_299_an=profile_287_an;
O3_retrieval.profile_299_avg=profile_299_avg;
O3_retrieval.profile_287_avg=profile_287_avg;
O3_retrieval.Date=datestr(DateTime(1),'yyyymmdd');
O3_retrieval.DateTime=DateTime;
O3_retrieval.TimeInHour=TimeInHour;
O3_retrieval.height=height;
O3_retrieval.TimeInHour_avg=TimeInHour_avg;
O3_retrieval.DateTime_avg=DateTime_avg;
% O3_retrieval.SmoothBinSize=framelen;
O3_retrieval.TimeAvgNum=nAvg;

O3_retrieval.sigmaOn=sigmaOn;% O3 cross section at 287 (cm^2/molecule)
O3_retrieval.sigmaOff=sigmaOff;% O3 cross section at 299 (cm^2/molecule)
O3_retrieval.d_sigma=d_sigma; % delta cross section (m^2/molecule)

O3_retrieval.ND_O3=N_O3_2;
O3_retrieval.height_avg=h_chm15k;
O3_retrieval.N_O3_interp=N_O3_interp;
O3_retrieval.D_molex_interp=D_molex_interp;
O3_retrieval.N_O3_bsc_interp=N_O3_bsc_interp;
O3_retrieval.D_aext_interp=D_aext_interp;
O3_retrieval.NDAir_m3_interp=NDAir_m3_interp;
save(['OL retrieval',O3_retrieval.Date,'.mat'],'O3_retrieval','mwr','chm15kRetrieval');