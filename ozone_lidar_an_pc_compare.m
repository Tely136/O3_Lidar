% Ozone lidar retrieval quick look, using standard atmosphere input 
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
pc_ad_ratio.disc287=8;
pc_ad_ratio.HV287=750;
pc_ad_ratio.disc299=8;
pc_ad_ratio.HV299=795;

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
profile_287_an=nan(nbin,numFiles);
profile_299_an=nan(nbin,numFiles);
profile_287_pc=nan(nbin,numFiles);
profile_299_pc=nan(nbin,numFiles);

profile_287_an_nr= nan(nbin,numFiles);
profile_287_pc_nr= nan(nbin,numFiles);
profile_299_an_nr= nan(nbin,numFiles);
profile_299_pc_nr= nan(nbin,numFiles);
% height
height=3.75*[1:1:nbin]';
hkm=height/1000;
% datetime array (local time)
DateTime=NaT(1,numFiles);
% time in hour (local time)
TimeInHour=nan(1,numFiles);
td=280;
bgbin=100;
% Loop over all files reading them in
for k = 1 : numFiles
    %fprintf('Now reading file %s\n', fullFileNames{k});
    % read both data and metadata from the filestore
    Data=read(fds);
    profile_287_an(:,k)= Data.data(:,1)-mean(Data.data(end-bgbin:end,1));
    profile_299_an(:,k)= Data.data(:,3)-mean(Data.data(end-bgbin:end,3));
     % dead time correation
    profile_287_pc(:,k)= Data.data(:,2)-mean(Data.data(end-bgbin:end,2));
    profile_287_pc(:,k)= profile_287_pc(:,k)./(1-profile_287_pc(:,k)/td);
    profile_299_pc(:,k)= Data.data(:,4)-mean(Data.data(end-bgbin:end,4));
    profile_299_pc(:,k)= profile_299_pc(:,k)./(1-profile_299_pc(:,k)/td);
    
    profile_287_an_nr(:,k)= Data.data(:,5)-mean(Data.data(end-bgbin:end,5));
    profile_287_pc_nr(:,k)= Data.data(:,6)-mean(Data.data(end-bgbin:end,6));
    profile_299_an_nr(:,k)= Data.data(:,7)-mean(Data.data(end-bgbin:end,7));
    profile_299_pc_nr(:,k)= Data.data(:,8)-mean(Data.data(end-bgbin:end,8));
   
    temp=cell2mat(Data.textdata(2,1));
    infmt='dd/MM/yyyy HH:mm:ss';
    DateTime(k)=datetime(temp(10:28),'InputFormat',infmt);
    TimeInHour(k)=hour(DateTime(k))+minute(DateTime(k))/60+second(DateTime(k))/3600;
end
%% Sigle Signal Profile
id=1;
figure
plot(hkm,profile_287_pc(:,id));hold on;
plot(hkm,profile_299_pc(:,id));
set(gca,'YScale','log')
ylabel('Analog signal (mV)');
% ylim([1e-3,200])
xlabel('Altitude (km)');
% xlim([0,8])
% ylim([1e-3,200])
legend('287nm - Analog, Far','299nm - Analog, Far','287nm - PC, Far','299nm - PC, Far')
title(['Ozone lidar signal (1min) at ',datestr(DateTime(id))])


%% % Plot a single profile ratio
% far range
% id=1;
% figure
% plot(hkm,profile_287_pc(:,id)./profile_287_an(:,id),'b');hold on;
% plot(hkm,profile_299_pc(:,id)./profile_299_an(:,id),'r');
% legend('287nm','299nm')
% ylabel('log(Signal ratio) (a.u.)');
% xlabel('Altitude (km)');
% title(['Ozone lidar far range signal ratio (Analog/PC -1min) at ',datestr(DateTime(id))])
% ylim([0,200])
% xlim([0,8])
% 
%% near range
% id=1;
% figure
% plot(hkm,profile_287_pc_nr(:,id)./profile_287_an_nr(:,id),'b');hold on;
% plot(hkm,profile_299_pc_nr(:,id)./profile_299_an_nr(:,id),'r');
% legend('287nm','299nm')
% ylabel('log(Signal ratio) (a.u.)');
% xlabel('Altitude (km)');
% title(['Ozone lidar near range signal ratio (Analog/PC-1min) at ',datestr(DateTime(id))])
% ylim([0,200])
% xlim([0,3])
% 
%% Time averaging of the lidar signal
[len_height,len_time]=size(profile_287_an);
nAvg=10;% 5min average

[profile_287_an_avgT,~]=MoveTimeAve(profile_287_an,len_height,len_time,nAvg);
[profile_287_pc_avgT,~]=MoveTimeAve(profile_287_pc,len_height,len_time,nAvg);
[profile_299_an_avgT,~]=MoveTimeAve(profile_299_an,len_height,len_time,nAvg);
[profile_299_pc_avgT,~]=MoveTimeAve(profile_299_pc,len_height,len_time,nAvg);

[profile_287_an_nr_avgT,~]=MoveTimeAve(profile_287_an_nr,len_height,len_time,nAvg);
[profile_287_pc_nr_avgT,~]=MoveTimeAve(profile_287_pc_nr,len_height,len_time,nAvg);
[profile_299_an_nr_avgT,~]=MoveTimeAve(profile_299_an_nr,len_height,len_time,nAvg);
[profile_299_pc_nr_avgT,TimeIndArr]=MoveTimeAve(profile_299_pc_nr,len_height,len_time,nAvg);
TimeInHour_avg=TimeInHour(TimeIndArr);
DateTime_avg=DateTime(TimeIndArr);
[m,N]=size(profile_287_an_avgT);

movmeanLen=40;
profile_287_an_avg=movmean(profile_287_an_avgT,movmeanLen,1);
profile_287_pc_avg=movmean(profile_287_pc_avgT,movmeanLen,1);
profile_299_an_avg=movmean(profile_299_an_avgT,movmeanLen,1);
profile_299_pc_avg=movmean(profile_299_pc_avgT,movmeanLen,1);

profile_287_an_nr_avg=movmean(profile_287_an_nr_avgT,movmeanLen,1);
profile_287_pc_nr_avg=movmean(profile_287_pc_nr_avgT,movmeanLen,1);
profile_299_an_nr_avg=movmean(profile_299_an_nr_avgT,movmeanLen,1);
profile_299_pc_nr_avg=movmean(profile_299_pc_nr_avgT,movmeanLen,1);
%% Plot Time averaged AD and PC signal far channel 
id=1;
figure
plot(hkm,profile_287_an_avg(:,id));hold on;
plot(hkm,profile_299_an_avg(:,id));
ylabel('Analog signal (mV)');
ylim([1e-4,500])
set(gca,'YScale','log')
yyaxis right
plot(hkm,profile_287_pc_avg(:,id),'b--');hold on;
plot(hkm,profile_299_pc_avg(:,id),'r--');
set(gca,'YScale','log')
ylabel('PC signal (MHz)');
xlabel('Altitude (km)');
xlim([0,22.5])
%ylim([1e-3,500])
legend('287nm - AD, Far, HV=800V','299nm - AD, Far, HV=795V','287nm - PC, Far, HV=800, disc=4','299nm - PC, Far,  HV=795, disc=4')
title(['Ozone lidar signal (10min-Avg) at ',datestr(DateTime_avg(id))])

%% Plot Time averaged PC/AD signal ratio far channel 


% id=1;
% pc_ad_ratio.an287=profile_287_an_avg(:,id);
% pc_ad_ratio.an299=profile_299_an_avg(:,id);
% pc_ad_ratio.pc287=profile_287_pc_avg(:,id);
% pc_ad_ratio.pc299=profile_299_pc_avg(:,id);

% figure
% plot(hkm,profile_287_pc_avgT(:,id)./profile_287_an_avgT(:,id));hold on;
% plot(hkm,profile_299_pc_avgT(:,id)./profile_299_an_avgT(:,id));
% set(gca,'YScale','log')
% ylabel('PC/AD ratio');
% xlabel('Altitude (km)');
% xlim([0,10])
% ylim([0,200])
% grid on
% legend('287nm - Far','299nm - Far')
% title(['Ozone lidar Far Range PC/AD signal ratio (10min-Avg) at ',datestr(DateTime_avg(id))])
% 
% savefilename=['pc_ad_r',datestr(DateTime_avg(id),'yyyymmddTHHMM'),'.mat']
% save(savefilename,'pc_ad_ratio');

% %% Plot Time averaged PC/AD signal ratio Near channel
% id=1;
% figure
% plot(hkm,profile_287_an_nr_avgT(:,id),'b-');hold on;
% plot(hkm,profile_299_an_nr_avgT(:,id),'r-');
% ylabel('Analog signal (mV)');
% ylim([1e-3,200])
% set(gca,'YScale','log')
% yyaxis right
% plot(hkm,profile_287_pc_nr_avgT(:,id),'b-.');hold on;
% plot(hkm,profile_299_pc_nr_avgT(:,id),'r-.');
% set(gca,'YScale','log')
% ylabel('PC signal (MHz)');
% xlabel('Altitude (km)');
% xlim([0,3])
% %ylim([1e-3,200])
% legend('287nm - AD, Near','299nm - AD, Near','287nm - PC, Near','299nm - PC, Near')
% title(['Ozone lidar near range signal (10min-Avg) at ',datestr(DateTime_avg(id))])

%% Plot time averaged signal ratio(PC/AD) different time
% id=1;
% figure
% plot(hkm,profile_287_pc_nr_avgT(:,id)./profile_287_an_nr_avgT(:,id));hold on;
% plot(hkm,profile_299_pc_nr_avgT(:,id)./profile_299_an_nr_avgT(:,id));
% % set(gca,'YScale','log')
% ylabel('PC/AD ratio');
% xlabel('Altitude (km)');
% xlim([0,3])
% ylim([0,200])
% grid on
% legend('287nm - Near','299nm - Near')
% title(['Ozone lidar near range PC/AD signal ratio (10min-Avg) at ',datestr(DateTime_avg(id))])
% %% Plot time averaged PC AD signal different time
% id=1;
% hold on
% % plot(hkm,profile_287_an_nr_avgT(:,id),'b-');
% % plot(hkm,profile_299_an_nr_avgT(:,id),'r-');
% plot(hkm,profile_287_pc_nr_avgT(:,id),'b--');
% plot(hkm,profile_299_pc_nr_avgT(:,id),'r--');
% legend('287nm - AD, Near, disc =6','299nm - AD, Near, disc =6','287nm - PC, Near, disc =6','299nm - PC, Near, disc =6')
%% Line relationship of AD and PC signal 

% id=1;
% z1=floor(4000/3.75);
% z2=floor(16000/3.75);
% [p,S]=polyfit(profile_287_an_avgT(z1:z2,id),profile_287_pc_avgT(z1:z2,id),1);
% figure
% plot(profile_287_pc_avgT(z1:z2,id),profile_287_an_avgT(z1:z2,id),'.');hold on;
% %plot(profile_299_pc_avgT(z1:z2,id),profile_299_an_avgT(z1:z2,id),'.');
% % set(gca,'YScale','log')
% ylabel('Analog (mV)');
% xlabel('PC (MHz)');
% legend('287nm - Far','299nm - Far')
% title(['Ozone lidar far range PC vs. AD signal (10min-Avg) at ',datestr(DateTime_avg(id))])
