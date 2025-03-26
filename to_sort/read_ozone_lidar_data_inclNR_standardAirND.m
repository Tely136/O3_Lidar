 clear all; close all;
% Get a list of all txt files in the current folder, or subfolders of it.
% Path of the folder that stores all the profile data of a selected date(.txt)
folderPath='C:\Users\OzoneLidar\Documents\Lidar Data\20220415\txt';
savepath='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\ozone lidar data\ozone lidar results\';
% microradiometer data file of the selected date 
mwrfilename='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\mwrdata\2021-12-16_00-04-06_lv2.csv';
% Aeronet AOD data file contains the AOD of the selected date 
aodfile='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\AeronetAoddata\20211201_20211215_CCNY.txt';

fds = fileDatastore(folderPath, 'ReadFcn', @importdata,'FileExtensions',{'.txt'});
fullFileNames = fds.Files;
numFiles = length(fullFileNames);

%----------Sample data file looks like below-------------------
% te2192013.483376
% CCNY     20/09/2021 13:47:34 20/09/2021 13:48:33 0100 -073.950000 0040.821000 00.0 00.0
% 0001201 0020 0000000 0020 04 0000000 0020 0000000 0000
% 1 0 1 04000 1 0000 3.75 00287.o 0 0 01 000 16 001201 0.500 BT0
% 1 1 1 04000 1 0000 3.75 00287.o 0 0 00 000 00 001201 3.1746 BC0
% 1 0 1 04000 1 0000 3.75 00299.o 0 0 01 000 16 001201 0.500 BT1
% 1 1 1 04000 1 0000 3.75 00299.o 0 0 00 000 00 001201 3.1746 BC1
%--------------------------------------------------------------

% profile raw data
nbin=6000;
profile_287_an_raw=nan(nbin,numFiles);
profile_299_an_raw=nan(nbin,numFiles);
profile_287_nr_an_raw=nan(nbin,numFiles);
profile_299_nr_an_raw=nan(nbin,numFiles);
% background subtracted data
profile_287_an=nan(nbin,numFiles);% 
profile_299_an=nan(nbin,numFiles);
profile_287_nr_an=nan(nbin,numFiles);% 
profile_299_nr_an=nan(nbin,numFiles);
% height in m
height=3.75*[1:1:nbin]';
hkm=height/1000;
% datetime array (local time)
DateTime=NaT(1,numFiles);
% time in hour (local time)
TimeInHour=nan(1,numFiles);
bgbin=200;
% Loop over all files reading them in
for k = 1 : numFiles
    %fprintf('Now reading file %s\n', fullFileNames{k});
    % read both data and metadata from the filestore
    Data=read(fds);
    profile_287_an_raw(:,k)= Data.data(:,1);
    profile_299_an_raw(:,k)= Data.data(:,3);
    profile_287_an(:,k)= profile_287_an_raw(:,k)-mean(profile_287_an_raw(nbin-bgbin:nbin,k));
    profile_299_an(:,k)= profile_299_an_raw(:,k)-mean(profile_299_an_raw(nbin-bgbin:nbin,k));
    
    profile_287_nr_an_raw(:,k)= Data.data(:,5);
    profile_299_nr_an_raw(:,k)= Data.data(:,7);
    profile_287_nr_an(:,k)= profile_287_nr_an_raw(:,k)-mean(profile_287_nr_an_raw(nbin-bgbin:nbin,k));
    profile_299_nr_an(:,k)= profile_299_nr_an_raw(:,k)-mean(profile_299_nr_an_raw(nbin-bgbin:nbin,k));

    temp=cell2mat(Data.textdata(2,1));
    infmt='dd/MM/yyyy HH:mm:ss';
    DateTime(k)=datetime(temp(10:28),'InputFormat',infmt);
    TimeInHour(k)=hour(DateTime(k))+minute(DateTime(k))/60+second(DateTime(k))/3600;
end
% profile_287_an= movmedian(profile_287_an,21,1);
% profile_299_an= movmedian(profile_299_an,21,1);
%% Plot a single profile
id=1;
figure
plot(profile_287_an(:,id)./profile_287_nr_an(:,id),hkm,'r');hold on;
plot(profile_299_an(:,id)./profile_299_nr_an(:,id),hkm,'b');
legend('287nm','299nm')
xlabel('Signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar signal ratio (Far/Near) at ',datestr(DateTime(id))])
ylim([0,2])
xlim([0,200])
% 
% id=1;
% figure
% plot(profile_287_an(:,id),hkm,'r');hold on;
% plot(profile_299_an(:,id),hkm,'b');
% plot(profile_287_nr_an(:,id),hkm,'r--');
% plot(profile_299_nr_an(:,id),hkm,'b--');
% set(gca,'XScale','log')
% legend('287nm Far','299nm Far','287nm Near','299nm Near')
% xlabel('Signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar signal at ',datestr(DateTime(id))])
% ylim([0,2])
% xlim([0,200])

%% Plot the signal color image
% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour,hkm,profile_287_an);
% set(gca,'YDir','normal');
% set(gca,'FontSize',14);
% set(gca,'ColorScale','log');
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar 287nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
% 
% subplot(1,2,2)
% I=imagesc(TimeInHour,hkm,profile_299_an);
% set(gca,'YDir','normal');
% set(gca,'ColorScale','log');
% set(gca,'FontSize',14);
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar 299nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);


%% Plot the range corrected
[m,n]=size(profile_287_an_raw);
pz2_287=nan(size(profile_287_an));
pz2_299=nan(size(profile_299_an));
pz2_287_nr=nan(size(profile_287_nr_an));
pz2_299_nr=nan(size(profile_299_nr_an));
for i=1:n
   pz2_287(:,i)= profile_287_an(:,i).*(hkm).^2;
   pz2_299(:,i)= profile_299_an(:,i).*(hkm).^2;
   pz2_287_nr(:,i)= profile_287_nr_an(:,i).*(hkm).^2;
   pz2_299_nr(:,i)= profile_299_nr_an(:,i).*(hkm).^2;
end
% 
% id=1;
% figure
% plot(pz2_287(:,id),hkm,'r');hold on;
% plot(pz2_299(:,id),hkm,'b');
% plot(pz2_287_nr(:,id),hkm,'r--');
% plot(pz2_299_nr(:,id),hkm,'b--');
% set(gca,'XScale','log');
% legend('287nm Far','299nm Far','287nm Near','299nm Near')
% xlabel('Range corrected signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime(id))])
%ylim([0,2])

% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour,hkm,pz2_287);
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
% I=imagesc(TimeInHour,hkm,pz2_299);
% set(gca,'YDir','normal');
% set(gca,'ColorScale','log');
% set(gca,'FontSize',14);
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone lidar 299nm P\cdotz^2 (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);

%% Time averaging of the lidar signal
nAvg=10;% 10min average

N=floor(n/nAvg);
profile_287_avgT=nan(m,N);
profile_299_avgT=nan(m,N);
profile_287_avg=nan(m,N);
profile_299_avg=nan(m,N);

profile_287_nr_avgT=nan(m,N);
profile_299_nr_avgT=nan(m,N);
profile_287_nr_avg=nan(m,N);
profile_299_nr_avg=nan(m,N);
TimeInHour_avg=nan(1,N);
DateTime_avg=NaT(1,N);

for i=1:N
    ind_t=(i-1)*nAvg;
    profile_287_avgT(:,i)=mean(profile_287_an(:,1+ind_t:nAvg+ind_t),2,'omitnan');
    profile_299_avgT(:,i)=mean(profile_299_an(:,1+ind_t:nAvg+ind_t),2,'omitnan');
    profile_287_nr_avgT(:,i)=mean(profile_287_nr_an(:,1+ind_t:nAvg+ind_t),2,'omitnan');
    profile_299_nr_avgT(:,i)=mean(profile_299_nr_an(:,1+ind_t:nAvg+ind_t),2,'omitnan');

    TimeInHour_avg(i)=mean(TimeInHour(1+ind_t:nAvg+ind_t));
    DateTime_avg(i)=DateTime(5+ind_t);
end 
% id=1;
% figure
% plot(profile_287_nr_avgT(:,id),hkm,'r');hold on
% plot(profile_299_nr_avgT(:,id),hkm,'b');
% plot(profile_287_avgT(:,id),hkm,'r--');
% plot(profile_299_avgT(:,id),hkm,'b--');
% set(gca,'XScale','log')
% legend('287nm Near','299nm Near','287nm Far','299nm Far')
% xlabel('Signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar signal at ',datestr(DateTime_avg(id))])
% ylim([0,3])
% grid on

% prof_287nr6=profile_287_nr_avgT(:,id);
% prof_299nr6=profile_299_nr_avgT(:,id);
% prof_287fr6=profile_287_avgT(:,id);
% prof_299fr6=profile_299_avgT(:,id);


% id1=1;
% id2=2;
% figure
% plot(profile_287_avgT(:,id2)./profile_287_avgT(:,id1),hkm,'r');hold on;
% plot(profile_299_avgT(:,id2)./profile_299_avgT(:,id1),hkm,'b');
% plot(profile_287_nr_avgT(:,id2)./profile_287_nr_avgT(:,id1),hkm,'r--');hold on
% plot(profile_299_nr_avgT(:,id2)./profile_299_nr_avgT(:,id1),hkm,'b--');
% legend('287nm Far signal ratio','299nm Far signal ratio','287nm Near signal ratio','299nm Near signal ratio')
% xlabel('Signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar signal ratio (high vol/low vol) at ',datestr(DateTime(id2))])
% ylim([0,15])
%% Vertical filtering of the lidar signal using 1st-order Savitzky-Golay (boxcar smoothing)
% 
% framelen = 21;% smooth bins
% 
% profile_287_avg=sgolayfilt(profile_287_avgT,order,framelen);
% profile_299_avg=sgolayfilt(profile_299_avgT,order,framelen); 

order = 1;% 
movnum1=31;% 78.75m 0-2km 1:533
start_bin_nr=find(hkm>0.05,1);
start_bin=find(hkm>0.78,1);


profile_287_nr_avg(start_bin_nr:end,:)=sgolayfilt(sgolayfilt(profile_287_nr_avgT(start_bin_nr:end,:),order,movnum1),order,movnum1);
profile_299_nr_avg(start_bin_nr:end,:)=sgolayfilt(sgolayfilt(profile_299_nr_avgT(start_bin_nr:end,:),order,movnum1),order,movnum1);



movnum1=31;% 78.5m 0-2km 1:533
movnum2=53;% 200m  2-5km 534:1333
movnum3=81;% 300m  >5km  1334:end 
temp=nan(size(profile_287_avgT));
temp(start_bin:end,:)=sgolayfilt(profile_287_avgT(start_bin:end,:),order,movnum1);
profile_287_avg(1:533,:)=temp(1:533,:);
temp=sgolayfilt(profile_287_avgT,order,movnum2);
profile_287_avg(534:1333,:)=temp(534:1333,:);
temp=sgolayfilt(profile_287_avgT,order,movnum3);
profile_287_avg(1334:end,:)=temp(1334:end,:);

temp(start_bin:end,:)=sgolayfilt(profile_299_avgT(start_bin:end,:),order,movnum1); 
profile_299_avg(1:533,:)=temp(1:533,:);
temp=sgolayfilt(profile_299_avgT,order,movnum2); 
profile_299_avg(534:1333,:)=temp(534:1333,:);
temp=sgolayfilt(profile_299_avgT,order,movnum3); 
profile_299_avg(1334:end,:)=temp(1334:end,:);

id=1;
figure
plot(profile_287_nr_avg(:,id),hkm,'r');hold on
plot(profile_299_nr_avg(:,id),hkm,'b');
plot(profile_287_avg(:,id),hkm,'r--');
plot(profile_299_avg(:,id),hkm,'b--');
set(gca,'XScale','log')
legend('287nm Near','299nm Near','287nm Far','299nm Far')
xlabel('Signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar signal at ',datestr(DateTime_avg(id))])
ylim([0,3])
grid on

% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour_avg,hkm,profile_287_avg);
% set(gca,'ColorScale','log');
% set(gca,'YDir','normal','FontSize',14);
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Smoothed Ozone lidar 287nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
% 
% subplot(1,2,2)
% I=imagesc(TimeInHour_avg,hkm,profile_299_avg);
% set(gca,'ColorScale','log');
% set(gca,'YDir','normal','FontSize',14);
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Smoothed Ozone lidar 299nm signal (a.u.) at ',datestr(DateTime(1),'yy/mm/dd')]);
%% Plot the smoothed range corrected signal
pz2_287_avg=nan(size(profile_287_avg));
pz2_299_avg=nan(size(profile_299_avg));
pz2_287_nr_avg=nan(size(profile_287_avg));
pz2_299_nr_avg=nan(size(profile_299_avg));
for i=1:N
   pz2_287_avg(:,i)= profile_287_avg(:,i).*(hkm).^2;
   pz2_299_avg(:,i)= profile_299_avg(:,i).*(hkm).^2;
   pz2_287_nr_avg(:,i)= profile_287_nr_avg(:,i).*(hkm).^2;
   pz2_299_nr_avg(:,i)= profile_299_nr_avg(:,i).*(hkm).^2;
end

% id=1;
% figure
% semilogx(pz2_287_avg(:,id),hkm,'r');hold on;
% semilogx(pz2_299_avg(:,id),hkm,'b');
% semilogx(pz2_287_nr_avg(:,id),hkm,'r--');
% semilogx(pz2_299_nr_avg(:,id),hkm,'b--');
% legend('287nm-Far range','299nm-Far range','287nm-Near range','299nm-Near range')
% xlabel('Range corrected signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar smoothed P\cdotz^2 at ',datestr(DateTime(id))])
% ylim([0,2])
% 
% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour_avg,hkm,pz2_287_avg,[0,8]);
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
% I=imagesc(TimeInHour_avg,hkm,pz2_299_avg,[0,8]);
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
sigmaOn=203.4*10^(-20);% O3 cross section at 287.2 (cm^2/molecule)
sigmaOff=45.51*10^(-20);% O3 cross section at 299.1 (cm^2/molecule)
d_sigma=(sigmaOn-sigmaOff)*(10^-2)^2; % delta cross section (m^2/molecule)


%dR=frameLen*3.75; % range resolution (m)
Poff=profile_299_avg;
Pon=profile_287_avg;
ratio_P=Pon./Poff;
ratio_P(1:start_bin-1,:)=nan;
ratio_P(start_bin:end,:)=sgolayfilt(movmean(ratio_P(start_bin:end,:),movnum3,1),order,movnum2);
ratio_P(ratio_P<0)=nan;
Ln_ratio_P=log(ratio_P);

order = 1;% 
movnum1=31;% 101.25m 0-2km 1:533
Poff_nr=profile_299_nr_avg;
Pon_nr=profile_287_nr_avg;
ratio_P_nr=Pon_nr./Poff_nr;
ratio_P_nr(1:start_bin_nr-1,:)=nan;
ratio_P_nr=sgolayfilt(movmean(ratio_P_nr,movnum3,1),order,movnum2);
ratio_P_nr(ratio_P_nr<0)=nan;
Ln_ratio_P_nr=log(ratio_P_nr);


 id=1;
% figure
% semilogx(ratio_P(:,id),hkm,'r');hold on;
% semilogx(ratio_P_nr(:,id),hkm,'b');
% legend('P_{287}/P_{299}-Far range','P_{287}/P_{299}-Near range')
% xlabel('P_{on}/P_{off} ratio (a.u.)');
% ylabel('Altitude (km)');
% title(['P_{off}/P_{off} ratio at ',datestr(DateTime(id))])
% ylim([0,2])
% % 
figure
plot(Ln_ratio_P(:,id),hkm,'r');hold on;
plot(Ln_ratio_P_nr(:,id),hkm,'b');
legend('Log(P_{287}/P_{299})-Far range','Log(P_{287}/P_{299})-Near range')
xlabel('Log(P_{on}/P_{off}) ratio (a.u.)');
ylabel('Altitude (km)');
title(['Log(P_{off}/P_{off}) ratio at ',datestr(DateTime(id))])
ylim([0,2])

% Using covolution with the SG filter to calculate the first order
% derivative 


frameLen1=31;% 78.5m 0-2km 1:533
dh1=533;
frameLen2=53;% 198.75m  2-5km 534:1333
dh2=1333;
frameLen3=81;% 303.75m  >5km  1334:end 
[b,g1] = sgolay(2,frameLen1);
[b,g2] = sgolay(2,frameLen2);
[b,g3] = sgolay(2,frameLen3);
dR=3.75;
N_O3_SG_fr=nan(size(ratio_P));
diff0=nan(size(ratio_P,1),1);

for i=1:N
  diff1= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g1(:,2), 'same');
  diff2= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g2(:,2), 'same');
  diff3= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g3(:,2), 'same');
  diff0(1:dh1)=diff1(1:dh1);
  diff0(dh1+1:dh2)=diff2(dh1+1:dh2);
  diff0(dh2+1:end)=diff3(dh2+1:end);
  N_O3_SG_fr(:,i)=1/(2*d_sigma).*diff0;
end


[b,g] = sgolay(2,31);
N_O3_SG_nr=nan(size(ratio_P_nr));
diff0=nan(size(ratio_P_nr,1),1);

for i=1:N
  diff0= conv(Ln_ratio_P_nr(:,i), factorial(1)/(dR)^1 * g(:,2), 'same');
  N_O3_SG_nr(:,i)=1/(2*d_sigma).*diff0;
end


% id=2
% figure
% plot(N_O3_SG_fr(:,id)/1e6,hkm,'LineWidth',1.2);hold on
% plot(N_O3_SG_nr(:,id)/1e6,hkm,'LineWidth',1.2);
% xlabel('Ozone number density (molecule/cm^3)')
% ylabel('Altitude (km)')
% legend('Far range','Near range')
% title(['Ozone number density using SG-filtering (molecule/cm^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,5])
% xlim([0,2e12])

% Concanation of two channel
hh=1.5;% the highest altitude of the concanation range 
hl=(start_bin+2*movnum1)*3.75/1000;% the lowest altitude of the concanation range
ind_hh=hkm<hh&hkm>hl;
N_O3_SG_concat=N_O3_SG_fr;
N_O3_SG_concat(~(hkm>hl),:)=N_O3_SG_nr(~(hkm>hl),:);
% merge range
h_concat=hkm(ind_hh);
% weight of high channel
w=(h_concat-h_concat(1))./(h_concat(end)-h_concat(1));
N_O3_SG_concat(ind_hh,:)=w.*(N_O3_SG_fr(ind_hh,:))+(1-w).*(N_O3_SG_nr(ind_hh,:));

N_O3_SG=N_O3_SG_concat;
% Plot the retrieved ozone number density (no correction) 
figure
I=imagesc(TimeInHour_avg,hkm,N_O3_SG_fr,[1e17,2.5e18]);
set(gca,'YDir','normal','FontSize',14);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density(Far Range) using SG-filtering (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0,1])

figure
I=imagesc(TimeInHour_avg,hkm,N_O3_SG_nr,[1e17,2.5e18]);
set(gca,'YDir','normal','FontSize',14);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG_nr))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density(Near range) using SG-filtering (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0,1])

figure
I=imagesc(TimeInHour_avg,hkm,N_O3_SG_concat,[1e17,2.5e18]);
set(gca,'YDir','normal','FontSize',10);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG_concat))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density(merged) using SG-filtering (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0,1])
% 


%% using Standard Atmosphere Modeling
% Pressure p0=101 325 N/m2 = 1013.25 hPa
% Density ρ0= 1.225 kg/m3
% Temperature T0= 288.15°K (15°C)
% Speed of sound a0= 340.294 m/sec
% Acceleration of gravity g0=9.80665 m/sec 2
% T= T0-6.5h(m)/1000 [0, 11000m] troposphere
% T=216.65K [11000,20000]
% P= P0 (1-0.0065 h/T0)^5.2561 [0,11000m]
% P = 226.32hPa*exp(-g/(R*216.65K)*(h-11000)) [11000m, 20000m]
% ρ = p/RT

% Constant
R=8.31446; %ideal gas constant(J/mol.K)
g=9.80665;% gravitational acceleration (m/s^2)
M=0.0289652;% molar mass of dry air 0.0289652 (kg/mol)
Na=6.02214*10^23;% [/mol]
% Standard T and P profile
p0=1013.25*1e2; % surface pressure 1013.25 hpa = 101325pa
T0=288.15;% surface temp (K)
T=T0-6.5*height/1000;
T(height>11000)=216.65;
P=((1-0.0065*height/T0).^5.2561).*p0;
P(height>11000)=226.32*1e2*exp(-g/(R*216.65)*(height(height>11000)-11000));
% figure
% plot(T,height/1000)
% standar air number density
NDAir_m3=Na/R.*P./T;

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
molex_1=sigma1*1e-4*NDAir_m3;   %% unit: m^-1 

delta_molex=(1-(lamda_2/lamda_1)^(-4))* molex_1;
D_molex=delta_molex./d_sigma;

id=4;
figure
plot((N_O3_SG_fr(:,id)-D_molex)./NDAir_m3*1e9,hkm,'LineWidth',1.2);hold on
plot((N_O3_SG_nr(:,id)-D_molex)./NDAir_m3*1e9,hkm,'LineWidth',1.2);
xlabel('Ozone concentration (ppb)')
ylabel('Altitude (km)')
legend('Far range','Near range')
title(['Ozone concentration (ppb) ',datestr(DateTime_avg(id),'yy/mm/dd HH:MM:ss')]);
ylim([0.2,5])
xlim([0,100])
grid on


%% molecular extinction correction
[rows,cols]=size(N_O3_SG_concat);
NDAir_m3_mat = repmat(NDAir_m3, 1,cols);
D_molex_mat = repmat(D_molex, 1,cols);

% figure
% subplot(1,2,1)
% I=imagesc(TimeInHour_avg,hkm,1e9*N_O3_SG_concat./NDAir_m3_mat,[0,100]);
% set(gca,'YDir','normal','FontSize',10);
% % set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(N_O3_SG_concat))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone concentration (ppb) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,1])
% 
% subplot(1,2,2)
% I=imagesc(TimeInHour_avg,hkm,1e9*D_molex_mat./NDAir_m3_mat,[0,10]);
% set(gca,'YDir','normal','FontSize',10);
% % set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(D_molex_mat))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone concentration molecular correction (ppb) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,1])

figure
I=imagesc(TimeInHour_avg,hkm,1e9*(N_O3_SG_concat-D_molex_mat)./NDAir_m3_mat,[0,100]);
set(gca,'YDir','normal','FontSize',10);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG_concat))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone concentration (ppb) molecular corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0,1])

% figure
% I=imagesc(TimeInHour_avg(13:16),hkm,1e9*(N_O3_SG_nr(:,13:16)-D_molex_mat(:,13:16))./NDAir_m3_mat(:,13:16),[0,100]);
% set(gca,'YDir','normal','FontSize',14);
% % set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(N_O3_SG_nr(:,13:16)))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone concentration (ppb) molecular corrected ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,1])
