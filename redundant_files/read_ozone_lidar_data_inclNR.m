clear all; close all;
% Get a list of all txt files in the current folder, or subfolders of it.
% Path of the folder that stores all the profile data of a selected date(.txt)
folderPath='C:\Users\OzoneLidar\Documents\Lidar Data\20220318\txt';
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
% Loop over all files reading them in
for k = 1 : numFiles
    %fprintf('Now reading file %s\n', fullFileNames{k});
    % read both data and metadata from the filestore
    Data=read(fds);
    profile_287_an_raw(:,k)= Data.data(:,1);
    profile_299_an_raw(:,k)= Data.data(:,3);
    profile_287_an(:,k)= profile_287_an_raw(:,k)-median(profile_287_an_raw(nbin-500:nbin,k));
    profile_299_an(:,k)= profile_299_an_raw(:,k)-median(profile_299_an_raw(nbin-500:nbin,k));
    
    profile_287_nr_an_raw(:,k)= Data.data(:,5);
    profile_299_nr_an_raw(:,k)= Data.data(:,7);
    profile_287_nr_an(:,k)= profile_287_nr_an_raw(:,k)-median(profile_287_nr_an_raw(nbin-500:nbin,k));
    profile_299_nr_an(:,k)= profile_299_nr_an_raw(:,k)-median(profile_299_nr_an_raw(nbin-500:nbin,k));

    temp=cell2mat(Data.textdata(2,1));
    infmt='dd/MM/yyyy HH:mm:ss';
    DateTime(k)=datetime(temp(10:28),'InputFormat',infmt);
    TimeInHour(k)=hour(DateTime(k))+minute(DateTime(k))/60+second(DateTime(k))/3600;
end
% profile_287_an= movmedian(profile_287_an,21,1);
% profile_299_an= movmedian(profile_299_an,21,1);
%% Plot a single profile
% id=1;
% figure
% plot(profile_287_an(:,id)./profile_287_nr_an(:,id),hkm,'r');hold on;
% plot(profile_299_an(:,id)./profile_299_nr_an(:,id),hkm,'b');
% legend('287nm','299nm')
% xlabel('Signal (a.u.)');
% ylabel('Altitude (km)');
% title(['Ozone lidar signal ratio (Far/Near) at ',datestr(DateTime(id))])
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
% 

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
id=1;
figure
plot(pz2_287(:,id),hkm,'r');hold on;
plot(pz2_299(:,id),hkm,'b');
plot(pz2_287_nr(:,id),hkm,'r--');
plot(pz2_299_nr(:,id),hkm,'b--');
set(gca,'XScale','log');
legend('287nm Far','299nm Far','287nm Near','299nm Near')
xlabel('Range corrected signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime(id))])
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

id=1;figure
plot(profile_287_nr_avgT(:,id),hkm,'r');hold on
plot(profile_299_nr_avgT(:,id),hkm,'b');
plot(profile_287_avgT(:,id),hkm,'r--');
plot(profile_299_avgT(:,id),hkm,'b--');
set(gca,'XScale','log')
legend('287nm Near','299nm Near','287nm Far','299nm Far')
xlabel('Signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar signal at ',datestr(DateTime(id))])
ylim([0,3])
grid on

% prof_287nr6=profile_287_nr_avgT(:,id);
% prof_299nr6=profile_299_nr_avgT(:,id);
% prof_287fr6=profile_287_avgT(:,id);
% prof_299fr6=profile_299_avgT(:,id);

id=1;
figure
plot(profile_287_avgT(:,id)./profile_287_nr_avgT(:,id),hkm,'r');hold on;
plot(profile_299_avgT(:,id)./profile_299_nr_avgT(:,id),hkm,'b');
legend('287nm','299nm')
xlabel('Signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar signal ratio (Far/Near) at ',datestr(DateTime(id))])
ylim([0,1])



%% Vertical filtering of the lidar signal using 1st-order Savitzky-Golay (boxcar smoothing)

% framelen = 21;% smooth bins
% 
% profile_287_avg=sgolayfilt(profile_287_avgT,order,framelen);
% profile_299_avg=sgolayfilt(profile_299_avgT,order,framelen); 

 order = 1;% 
 movnum1=21;% 78.75m 0-2km 1:533

profile_287_nr_avg=sgolayfilt(sgolayfilt(profile_287_nr_avgT,order,movnum1),order,movnum1);
profile_299_nr_avg=sgolayfilt(sgolayfilt(profile_299_nr_avgT,order,movnum1),order,movnum1);



movnum1=21;% 78.5m 0-2km 1:533
movnum2=53;% 200m  2-5km 534:1333
movnum3=81;% 300m  >5km  1334:end 

temp=sgolayfilt(profile_287_avgT,order,movnum1);
profile_287_avg(1:533,:)=temp(1:533,:);
temp=sgolayfilt(profile_287_avgT,order,movnum2);
profile_287_avg(534:1333,:)=temp(534:1333,:);
temp=sgolayfilt(profile_287_avgT,order,movnum3);
profile_287_avg(1334:end,:)=temp(1334:end,:);

temp=sgolayfilt(profile_299_avgT,order,movnum1); 
profile_299_avg(1:533,:)=temp(1:533,:);
temp=sgolayfilt(profile_299_avgT,order,movnum2); 
profile_299_avg(534:1333,:)=temp(534:1333,:);
temp=sgolayfilt(profile_299_avgT,order,movnum3); 
profile_299_avg(1334:end,:)=temp(1334:end,:);


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
ratio_P(ratio_P<0)=nan;
Ln_ratio_P=log(ratio_P);

order = 1;% 
movnum1=21;% 101.25m 0-2km 1:533
Poff_nr=profile_299_nr_avg;
Pon_nr=profile_287_nr_avg;
ratio_P_nr=Pon_nr./Poff_nr;
ratio_P_nr=sgolayfilt(movmean(ratio_P_nr,movnum1,1),order,movnum1);
ratio_P_nr(ratio_P_nr<0)=nan;
Ln_ratio_P_nr=log(ratio_P_nr);

id=4;
figure
semilogx(ratio_P(:,id),hkm,'r');hold on;
semilogx(ratio_P_nr(:,id),hkm,'b');
legend('P_{287}/P_{299}-Far range','P_{287}/P_{299}-Near range')
xlabel('P_{on}/P_{off} ratio (a.u.)');
ylabel('Altitude (km)');
title(['P_{off}/P_{off} ratio at ',datestr(DateTime(id))])
ylim([0,2])
% 
% figure
% plot(Ln_ratio_P(:,id),hkm,'r');hold on;
% plot(Ln_ratio_P_nr(:,id),hkm,'b');
% legend('Log(P_{287}/P_{299})-Far range','Log(P_{287}/P_{299})-Near range')
% xlabel('Log(P_{on}/P_{off}) ratio (a.u.)');
% ylabel('Altitude (km)');
% title(['Log(P_{off}/P_{off}) ratio at ',datestr(DateTime(id))])
% ylim([0,2])

% N_O3_PF=nan(m,N);
% frameLen1=21;% 78.5m 0-2km 1:533
% dh1=533;
% frameLen2=53;% 198.75m  2-5km 534:1333
% dh2=1333;
% frameLen3=81;% 303.75m  >5km  1334:end 
% for j=1:N
% for i=1:m
%     % calculate the coef of 2nd order poly fit 
%     if i<dh1+1
%     p_coef=polyfit(height(max(1,i-(frameLen1-1)/2):min(i+(frameLen1-1)/2,m)),...
%         Ln_ratio_P(max(1,i-(frameLen1-1)/2):min(i+(frameLen1-1)/2,m),j),1);
%     end
%     if i>dh1 && i<dh2+1
%      p_coef=polyfit(height(max(1,i-(frameLen2-1)/2):min(i+(frameLen2-1)/2,m)),...
%         Ln_ratio_P(max(1,i-(frameLen2-1)/2):min(i+(frameLen2-1)/2,m),j),1);
%     end
%     if i>dh2
%      p_coef=polyfit(height(max(1,i-(frameLen3-1)/2):min(i+(frameLen3-1)/2,m)),...
%         Ln_ratio_P(max(1,i-(frameLen3-1)/2):min(i+(frameLen3-1)/2,m),j),1);    
%     end
%     
%     diff=p_coef(1);
%     N_O3_PF(i,j)=1/(2*d_sigma)*(-diff);
%     
% end
% end
% Using covolution with the SG filter to calculate the first order
% derivative 


frameLen1=21;% 78.5m 0-2km 1:533
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


[b,g] = sgolay(2,21);
N_O3_SG_nr=nan(size(ratio_P_nr));
diff0=nan(size(ratio_P_nr,1),1);

for i=1:N
  diff0= conv(Ln_ratio_P_nr(:,i), factorial(1)/(dR)^1 * g(:,2), 'same');
  N_O3_SG_nr(:,i)=1/(2*d_sigma).*diff0;
end


id=1;
figure
plot(N_O3_SG_fr(:,id)/1e6,hkm,'LineWidth',1.2);hold on
plot(N_O3_SG_nr(:,id)/1e6,hkm,'LineWidth',1.2);
xlabel('Ozone number density (molecule/cm^3)')
ylabel('Altitude (km)')
legend('Far range','Near range')
title(['Ozone number density using SG-filtering (molecule/cm^3) ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0,5])
xlim([0,2e12])

% Concanation of two channel
hh=0.6;% the highest altitude of the concanation range 
hl=0.5;% the lowest altitude of the concanation range
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
I=imagesc(TimeInHour_avg,hkm,N_O3_SG_fr,[1e17,2e18]);
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
I=imagesc(TimeInHour_avg,hkm,N_O3_SG_nr,[1e17,2e18]);
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
I=imagesc(TimeInHour_avg,hkm,N_O3_SG_concat,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',10);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG_concat))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density(merged) using SG-filtering (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0,1])


%% molecular extinction correction 
% mwr=loadMWRdata4lidar(mwrfilename,N_O3_2,DateTime_avg,height);
mwr=loadMWRdata(mwrfilename);
chm15kRetrieval=CHM15kForwardRetrieval_new(mwr,aodfile);
%% interpolate the time and vertical of mwr and chm15k
t_chm15k=chm15kRetrieval.time_array_utc;% utc hour
h_chm15k=chm15kRetrieval.height; % in km
t_mwr=mwr.time_array_utc;% utc hour
h_mwr=mwr.height;% km

timediff=hours(5); % time difference between the local time and utc time (EDT = 4; EST =5)

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
molex_1=sigma1*1e-4*mwr.NDAir_m3;   %% unit: m^-1 

delta_molex=(1-(lamda_2/lamda_1)^(-4))* molex_1;
D_molex=delta_molex./d_sigma;

% aerosol total backscatter correction

S1= 8*pi/3;
prompt = 'Do you want to use the Aeronet Angstrom Exponent? (if no aeronet AE available or choosing no, ae=1.5 will be used)? Y/N [N]:';
str = input(prompt,'s');
if isempty(str)
    str = 'N';
end
if str=='N'|str=='n'|isnan(chm15kRetrieval.ae_1020_340)|chm15kRetrieval.ae_1020_340<1
    ae=1.5;
else
    ae=chm15kRetrieval.ae_1020_340;
end
absc_299=(1064/299)^ae*chm15kRetrieval.aero_bsa1064;
absc_287=(1064/287)^ae*chm15kRetrieval.aero_bsa1064;
mbsc_299=(1064/299)^4/S1*chm15kRetrieval.am_interp;
mbsc_287=(1064/287)^4/S1*chm15kRetrieval.am_interp;
totbsc_299off=absc_299+mbsc_299;
totbsc_287on=absc_287+mbsc_287;
ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
ln_bsc=log(ratio_totbsc_onoff);

% Using second order polyfit to calculate the derivative of log backscatter on-off ratio
bscframeLen1=5;% 75m 0-2km 1:133
bscdh1=133;
bscframeLen2=13;% 195m  2-5km 134:333
bscdh2=333;
bscframeLen3=21;% 315m  >5km  334:end
% N_O3_bsc_PF=nan(size(ln_bsc));
% [Mbsc,Nbsc]=size(ln_bsc);
% for j=1:Nbsc
% for i=1:Mbsc
%     % calculate the coef of 2nd order poly fit 
%     if i<bscdh1+1
%     p_coef=polyfit(h_chm15k(max(1,i-(bscframeLen1-1)/2):min(i+(bscframeLen1-1)/2,Mbsc)),...
%         ln_bsc(max(1,i-(bscframeLen1-1)/2):min(i+(bscframeLen1-1)/2,Mbsc),j),1);
%     end
%     if i>bscdh1 && i<bscdh2+1
%      p_coef=polyfit(h_chm15k(max(1,i-(bscframeLen2-1)/2):min(i+(bscframeLen2-1)/2,Mbsc)),...
%         ln_bsc(max(1,i-(bscframeLen2-1)/2):min(i+(bscframeLen2-1)/2,Mbsc),j),1);
%     end
%     if i>bscdh2
%      p_coef=polyfit(h_chm15k(max(1,i-(bscframeLen3-1)/2):min(i+(bscframeLen3-1)/2,Mbsc)),...
%         ln_bsc(max(1,i-(bscframeLen3-1)/2):min(i+(bscframeLen3-1)/2,Mbsc),j),1);    
%     end
%     
%     diff=p_coef(1);
%     N_O3_bsc_PF(i,j)=1/(2*d_sigma)*(-diff/1000);
%     
% end
% end
% Using SG filter to calculate the derivative of log backscatter on-off ratio

[b,g1] = sgolay(2,bscframeLen1);
[b,g2] = sgolay(2,bscframeLen2);
[b,g3] = sgolay(2,bscframeLen3);
dz=15;
N_O3_bsc_SG=nan(size(ln_bsc));
diff0=nan(size(ln_bsc,1),1);
for i=1:length(t_chm15k)
  diff1= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g1(:,2), 'same');
  diff2= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g2(:,2), 'same');
  diff3= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g3(:,2), 'same');
  diff0(1:bscdh1)=diff1(1:bscdh1);
  diff0(bscdh1+1:bscdh2)=diff2(bscdh1+1:bscdh2);
  diff0(bscdh2+1:end)=diff3(bscdh2+1:end);
  N_O3_bsc_SG(:,i)=1/(2*d_sigma).*diff0;
end

% figure
% I=imagesc(chm15kRetrieval.time_aeroprof,h_chm15k,-N_O3_bsc_PF,[-1e18,1e18]);
% set(gca,'YDir','normal','FontSize',14);
% %set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(N_O3_bsc_SG))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['O_3 number density \beta correction terms using 1st-polyfit (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% % ylim([0,5])
% 
% figure
% I=imagesc(chm15kRetrieval.time_aeroprof,h_chm15k,-N_O3_bsc_SG,[-1e18,1e18]);
% set(gca,'YDir','normal','FontSize',14);
% %set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(N_O3_bsc_SG))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['O_3 number density \beta correction terms using SG-filtering (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% % ylim([0,5])


% aerosol extinction correction term

delta_aext=1e-3*((1064/287.2)^ae-(1064/299.1)^ae)*chm15kRetrieval.aero_ext1064;% aerosol extinction in (/m)
D_aext=delta_aext./d_sigma;

% convert to the same time base of ozone lidar and vertical res as CHM15k
[M,N]=size(N_O3_SG);


NDAir_m3_interp=nan(M,N);
D_molex_interp=nan(M,N);
N_O3_bsc_interp=nan(M,N);
D_aext_interp=nan(M,N);
pbl_interp=nan(1,N);
cbh_interp=nan(1,N);
chm15k_absc_interp=nan(M,N);
for i=1:N  
       
    indt=isbetween(t_mwr-timediff,DateTime_avg(i)-minutes(5),DateTime_avg(i)+minutes(5));
    temp=mean(mwr.NDAir_m3(:,indt),2,'omitnan');
    NDAir_m3_interp(:,i)=interp1(h_mwr,temp,hkm,'pchip',nan);
    NDAir_m3_interp(hkm>h_mwr(end),i)=nan;
       
    temp=mean(D_molex(:,indt),2,'omitnan');
    D_molex_interp(:,i)=interp1(h_mwr,temp,hkm,'pchip',0);
    D_molex_interp(hkm>h_mwr(end),i)=0;
    
    indt=isbetween(t_chm15k-timediff,DateTime_avg(i)-minutes(5),DateTime_avg(i)+minutes(5));
    temp=mean(N_O3_bsc_SG(:,indt),2,'omitnan');
    N_O3_bsc_interp(:,i)=interp1(h_chm15k,temp,hkm,'pchip',0);
    N_O3_bsc_interp(hkm>h_chm15k(end),i)=0;
    if nnz(isnan(temp))>0
       N_O3_bsc_interp(hkm>h_chm15k(find(isnan(temp),1,'first')),i)=nan;
    end
    
    
    temp=mean(D_aext(:,indt),2,'omitnan');
    D_aext_interp(:,i)=interp1(h_chm15k,temp,hkm,'pchip',0);
    D_aext_interp(hkm>h_chm15k(end),i)=0;
    if nnz(isnan(temp))>0
    D_aext_interp(hkm>h_chm15k(find(isnan(temp),1,'first')),i)=nan;
    end
    
    pbl_interp(i)=mean(chm15kRetrieval.pbl_aero(indt),2,'omitnan');
    cbh_interp(i)=min(chm15kRetrieval.cbh_aero(indt));
    
    temp=mean(chm15kRetrieval.aero_bsa1064(:,indt),2,'omitnan');
    chm15k_absc_interp(:,i)=interp1(h_chm15k,temp,hkm,'pchip',0);
    chm15k_absc_interp(hkm>h_chm15k(end),i)=0;
    if nnz(isnan(temp))>0
    chm15k_absc_interp(hkm>h_chm15k(find(isnan(temp),1,'first')),i)=nan;
    end
    
end 

% Plot 299nm aerosol backscatter profile
% 
% id=87;
% figure
% plot(chm15kRetrieval.am_interp(:,id),h_chm15k,'r');hold on;
% legend('1064 nm')
% xlabel('molecular extinction (km^{-1}');
% ylabel('Altitude (km)');
% title(['Molecular Backscatter profile at ',datestr(t_chm15k(id))])
% ylim([0.05,5])
% 
% figure
% subplot(1,3,1)
% plot(totbsc_299off(:,id),h_chm15k,'r');hold on;
% plot(totbsc_287on(:,id),h_chm15k,'b');
% legend('287nm','299nm')
% xlabel('Totoal Backscatter (km^{-1}sr^{-1})');
% ylabel('Altitude (km)');
% title(['Totoal Backscatter profile at ',datestr(t_chm15k(id))])
% ylim([0.05,5])
% 
% 
% subplot(1,3,2)
% plot(ln_bsc(:,id),h_chm15k);
% xlabel('Ln(\beta_{299}/\beta_{287}) (a.u.)');
% ylabel('Altitude (km)');
% title(['ln(\beta_{299}/\beta_{287}) profile at ',datestr(t_chm15k(id))])
% ylim([0.05,5])
% 
% subplot(1,3,3)
% plot(N_O3_bsc_SG(:,id),h_chm15k);
% xlabel('$\frac{1}{2\delta\sigma_{O_3}}\cdot\frac{d ln(\beta_{299}/\beta_{287})}{dz} (/m^3)$','Interpreter','latex');
% ylabel('Altitude (km)');
% title(['$\frac{1}{2\delta\sigma_{O_3}}\cdot\frac{d ln(\beta_{299}/\beta_{287})}{dz} (/m^3)$ profile at ',datestr(t_chm15k(id))],'Interpreter','latex')
% ylim([0.05,5])



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

% figure
% I=imagesc(TimeInHour_avg,hkm,N_O3,[1e17,2e18]);
% set(gca,'YDir','normal','FontSize',14);
% set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(N_O3))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone number density using 1st-polyfit (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,5])

% figure
% I=imagesc(TimeInHour_avg,hkm,N_O3_SG,[1e17,2e18]);
% set(gca,'YDir','normal','FontSize',14);
% set(gca,'ColorScale','log');
% set(I,'AlphaData',~isnan(N_O3_2))
% colormap('jet')
% colorbar
% xlabel('Local Time (hour)')
% ylabel('Altitude (km)')
% title(['Ozone number density using SG-filtering (molecule/m^3) ',datestr(DateTime(1),'yy/mm/dd')]);
% ylim([0,5])

%%
f1=figure;
subplot(2,2,1)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}),no correction ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,hkm,D_molex_interp,[1e16,1e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_molex_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} correction term ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,hkm,-N_O3_bsc_interp,[-1e18,1e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_bsc_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \beta correction term ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,hkm,D_aext_interp,[1e16,1e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_aext_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{a} correction term ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

saveas(f1,[savepath,datestr(DateTime(1),'yyyymmdd'),'O3ND_corterm.fig']);

f2=figure;
subplot(2,2,1)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv),no correction ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,hkm,(D_molex_interp)./NDAir_m3_interp*1e9,[-20,20]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_molex_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} correction term ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,hkm,(-N_O3_bsc_interp)./NDAir_m3_interp*1e9,[-20,20]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
colormap('jet')
set(I,'AlphaData',~isnan(N_O3_bsc_interp))
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \beta correction term ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,hkm,(D_aext_interp)./NDAir_m3_interp*1e9,[-20,20]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(D_aext_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{a} correction term ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

saveas(f2,[savepath,datestr(DateTime(1),'yyyymmdd'),'O3ppbv_corterm.fig']);


%%
f3=figure;
subplot(2,2,1)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), no correction ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG-D_molex_interp,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG-D_molex_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG-D_molex_interp-N_O3_bsc_interp,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_bsc_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} and \beta corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG-D_molex_interp-N_O3_bsc_interp-D_aext_interp,[1e17,2e18]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_bsc_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone number density (m^{-3}), \alpha_{m} \beta and \alpha_{a} corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])
saveas(f3,[savepath,datestr(DateTime(1),'yyyymmdd'),'O3ND_corr.fig']);



f4=figure;
subplot(2,2,1)
I=imagesc(TimeInHour_avg,hkm,N_O3_SG./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',11);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv),no correction ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,2)
I=imagesc(TimeInHour_avg,hkm,(N_O3_SG-D_molex_interp)./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',11);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG-D_molex_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,3)
I=imagesc(TimeInHour_avg,hkm,(N_O3_SG-D_molex_interp-N_O3_bsc_interp)./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',11);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG-D_molex_interp-N_O3_bsc_interp))
colormap('jet')
colorbar
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} \beta corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])

subplot(2,2,4)
I=imagesc(TimeInHour_avg,hkm,(N_O3_SG-D_molex_interp-N_O3_bsc_interp-D_aext_interp)./NDAir_m3_interp*1e9,[0,100]);
set(gca,'YDir','normal','FontSize',11);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(N_O3_SG-D_molex_interp-N_O3_bsc_interp))
colormap('jet')
colorbar
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
legend('PBLH')
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['Ozone mixing ratio (ppbv), \alpha_{m} \beta and \alpha_{a} corrected ',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])
saveas(f4,[savepath,datestr(DateTime(1),'yyyymmdd'),'O3ppbv_corr.fig']);
%% Plot the aerosol retrieval with the same time
f5=figure;
I=imagesc(TimeInHour_avg,hkm,chm15k_absc_interp,[0,1e-3]);
set(gca,'YDir','normal','FontSize',12);
% set(gca,'ColorScale','log');
set(I,'AlphaData',~isnan(chm15k_absc_interp))
colormap('jet')
colorbar
hold on;
plot(TimeInHour_avg,pbl_interp,'m+')
plot(TimeInHour_avg,cbh_interp,'go')
legend('PBLH')
xlabel('Local Time (hour)')
ylabel('Altitude (km)')
title(['CHM15k aerosol backscatter at 1064nm (km^{-1}sr^{-1})',datestr(DateTime(1),'yy/mm/dd')]);
ylim([0.45,5])



% O3_retrieval.profile_287_an_raw=profile_287_an_raw;
% O3_retrieval.profile_299_an_raw=profile_299_an_raw;
% O3_retrieval.profile_287_an=profile_287_an;
% O3_retrieval.profile_299_an=profile_287_an;

O3_retrieval.profile_299_avg=profile_299_avg;
O3_retrieval.profile_287_avg=profile_287_avg;
O3_retrieval.Date=datestr(DateTime(1),'yyyymmdd');
O3_retrieval.DateTime=DateTime;
O3_retrieval.TimeInHour=TimeInHour;
O3_retrieval.hkm=hkm;
O3_retrieval.TimeInHour_avg=TimeInHour_avg;
O3_retrieval.DateTime_avg=DateTime_avg;
% O3_retrieval.SmoothBinSize=framelen;
O3_retrieval.TimeAvgNum=nAvg;

O3_retrieval.sigmaOn=sigmaOn;% O3 cross section at 287 (cm^2/molecule)
O3_retrieval.sigmaOff=sigmaOff;% O3 cross section at 299 (cm^2/molecule)
O3_retrieval.d_sigma=d_sigma; % delta cross section (m^2/molecule)

O3_retrieval.ND_O3=N_O3_SG;
O3_retrieval.hkm=hkm;
O3_retrieval.D_molex_interp=D_molex_interp;
O3_retrieval.N_O3_bsc_interp=N_O3_bsc_interp;
O3_retrieval.D_aext_interp=D_aext_interp;
O3_retrieval.NDAir_m3_interp=NDAir_m3_interp;


O3_retrieval.O3_ppbv=N_O3_SG./NDAir_m3_interp*1e9;
O3_retrieval.O3_ppbv_amcorr=(N_O3_SG-D_molex_interp)./NDAir_m3_interp*1e9;
O3_retrieval.O3_ppbv_am_ba_corr=(N_O3_SG-D_molex_interp-N_O3_bsc_interp)./NDAir_m3_interp*1e9;
O3_retrieval.O3_ppbv_am_ba_aa_corr=(N_O3_SG-D_molex_interp-N_O3_bsc_interp-D_aext_interp)./NDAir_m3_interp*1e9;
O3_retrieval.pbl_interp=pbl_interp;

save(['OL retrieval',O3_retrieval.Date,'.mat'],'O3_retrieval','mwr','chm15kRetrieval');


[~,indh]=min(abs(hkm-0.55));
O3_550m_no_correction=O3_retrieval.O3_ppbv(indh,:);
O3_550m_molext_corr=O3_retrieval.O3_ppbv_amcorr(indh,:);
O3_550m_molext_aerobsc_corr=O3_retrieval.O3_ppbv_am_ba_corr(indh,:);
O3_550m_molext_aerobsc_aeroext_corr=O3_retrieval.O3_ppbv_am_ba_aa_corr(indh,:);
O3_localTime=TimeInHour_avg;
O3_localdatetime=DateTime_avg;
O3ppbv550m=table(O3_localdatetime',O3_localTime',O3_550m_no_correction',O3_550m_molext_corr',O3_550m_molext_aerobsc_corr',O3_550m_molext_aerobsc_aeroext_corr');
save(['O3ppbv_550m_',O3_retrieval.Date,'.mat'],'O3ppbv550m');

[~,indh]=min(abs(hkm-0.2));
O3_200m_no_correction=O3_retrieval.O3_ppbv(indh,:);
O3_200m_molext_corr=O3_retrieval.O3_ppbv_amcorr(indh,:);
O3_200m_molext_aerobsc_corr=O3_retrieval.O3_ppbv_am_ba_corr(indh,:);
O3_200m_molext_aerobsc_aeroext_corr=O3_retrieval.O3_ppbv_am_ba_aa_corr(indh,:);
O3_localTime=TimeInHour_avg;
O3_localdatetime=DateTime_avg;
O3ppbv200m=table(O3_localdatetime',O3_localTime',O3_200m_no_correction',O3_200m_molext_corr',O3_200m_molext_aerobsc_corr',O3_200m_molext_aerobsc_aeroext_corr');
save(['O3ppbv_200m_',O3_retrieval.Date,'.mat'],'O3ppbv200m');
