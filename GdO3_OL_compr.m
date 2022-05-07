% ground ozone measurement and ozone lidar measurement
% 1-st column: O3 mixing ratio in ppb
% 8th Column: day/mon/year
% 9th column: hour:min:sec (EDT)local time
clear all; close all;
% Get a list of all txt files in the current folder, or subfolders of it.
% Path of the folder that stores all the profile data of a selected date(.txt)
folderPath='/Users/Tinker/Documents/MATLAB/ozonelidar/gd_ozone_measur/';
savepath='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar_results/20220429/';
FileName='20220429_O3_Lab.txt';
fullFileName = [folderPath,FileName];
%fprintf('Now reading file %s\n', fullFileName);
% read data into a table
T = readtable(fullFileName,'Format','%f %f %f %f %f  %{dd/MM/yy}D %{HH:mm:ss}D','Delimiter',',','HeaderLines', 0, 'ReadVariableNames', false);
ground_O3_ppb=T{:,2};% 1st col ozone ppbv
%last two cols are the datetime (local time)
T.Var6.Format='dd/MM/yy HH:mm:ss';
ground_O3_ppb_datetime=T{:,end-1} + timeofday(T{:,end});
%% load ozone lidar measurement 
disp('Select the ozone lidar retrieval at 550m file(.mat)')
instruction='Select the ozone lidar retrieval file (.mat)';
[file, path] = uigetfile( ...
{  '*.mat', 'Ozone lidar retrieval result selector'; ...
   '*.*',  'All Files (*.*)', ... % Add to this to support other files (change filterindex number below)
   },...
   instruction);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
   load([path file]);
end

ind_gd_ozone=isbetween(ground_O3_ppb_datetime,O3ppbv550m{1,1}-minutes(10),O3ppbv550m{end,1}+minutes(10));
f1=figure;
plot(ground_O3_ppb_datetime(ind_gd_ozone),ground_O3_ppb(ind_gd_ozone),'g-');
hold on;
plot(O3ppbv550m{:,1},O3ppbv550m{:,3},'bo-');
set(gca,'LineWidth',1,'FontSize',13);
legend('in-situ','O_3-Lidar measurement at 550m (no correction)');
xlabel('Local Time (Hour)');
ylabel('Ozone Mixing Ratio (ppbv)');
ylim([0,100]);
title(datestr(O3ppbv550m{1,1}, 'yyyy-mmm-dd'))
grid on;
saveas(f1,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_noCorr.fig']);
saveas(f1,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_noCorr.jpg']);

f2=figure;
plot(ground_O3_ppb_datetime(ind_gd_ozone),ground_O3_ppb(ind_gd_ozone),'g-');
hold on;
plot(O3ppbv550m{:,1},O3ppbv550m{:,4},'bo-');
set(gca,'LineWidth',1,'FontSize',13);
legend('in-situ','O_3-Lidar measurement at 550m (\alpha_m corr)');
xlabel('Local Time (Hour)');
ylabel('Ozone Mixing Ratio (ppbv)');
ylim([0,100]);
title(datestr(O3ppbv550m{1,1}, 'yyyy-mmm-dd'));
grid on;
saveas(f2,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_mcorr.fig']);
saveas(f2,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_mcorr.jpg']);

f3=figure;
plot(ground_O3_ppb_datetime(ind_gd_ozone),ground_O3_ppb(ind_gd_ozone),'g-');
hold on;
plot(O3ppbv550m{:,1},O3ppbv550m{:,5},'bo-');
set(gca,'LineWidth',1,'FontSize',13);
legend('in-situ','O_3-Lidar measurement at 550m (\alpha_m and \beta_a corr)');
xlabel('Local Time (Hour)');
ylabel('Ozone Mixing Ratio (ppbv)');
title(datestr(O3ppbv550m{1,1}, 'yyyy-mmm-dd'));
ylim([0,100]);
grid on;
saveas(f3,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_m_absc_corr.fig']);
saveas(f3,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_m_absc_corr.jpg']);


f4=figure;
plot(ground_O3_ppb_datetime(ind_gd_ozone),ground_O3_ppb(ind_gd_ozone),'g-');
hold on;
plot(O3ppbv550m{:,1},O3ppbv550m{:,6},'bo-');
set(gca,'LineWidth',1,'FontSize',13);
legend('in-situ','O_3-Lidar measurement at 550m (\alpha_m, \beta_a and \alpha_a corr)');
xlabel('Local Time (Hour)');
ylabel('Ozone Mixing Ratio (ppbv)');
ylim([0,100]);
title(datestr(O3ppbv550m{1,1}, 'yyyy-mmm-dd'));
grid on;
saveas(f4,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_m_absc_aext_corr.fig']);
saveas(f4,[savepath,datestr(O3ppbv550m{1,1},'yyyymmdd'),'_O3ppbv_GdCompr_m_absc_aext_corr.jpg']);

