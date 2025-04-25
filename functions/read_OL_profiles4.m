%% Function: read_OL_profiles()
% Read ozone lidar signal profiles from the folder which contains a list of
%% Input: 
% folder_path: string, the folder path which saves the data .txt files 
% save_path: string, the path to save the results
% nbins: int 1x1, number of bins of each profile
% dz: float 1x1 the height in meter of each bin, for example, dz= 3.75m a 40MHz system
% bgbins: int 1x1, number of bins to be averaged at the end of the profile as background
% td: float 1x1 dead time constant, to calculate the true count rate Ct=Cm/(1-td*Cm)
% nAvg: int 1x1 number of bins to perform time average 

%% Output: 
%  
%  DateTime: datetime, Dim-(1, numFiles), datetime of each profiles
%  TimeInHour: float, Dim-(1, numFiles), time in hour of each profiles
%  Height: float, Dim- (nbin,1), height in meters
%  hkm: float, Dim- (nbin,1), height in km.
%  rawprof: structure, 8 fields, the raw signal profiles
% 
%  rawprof.an287: float, Dim-(nbin, numFiles), 287nm-Far channel analog raw data profile.
%  rawprof.an299: float, Dim-(nbin, numFiles), 299nm-Far channel analog raw data profile.
%  rawprof.an287nr: float, Dim-(nbin, numFiles), 287nm-Near channel analog raw data profile.
%  rawprof.an299nr: float, Dim-(nbin, numFiles), 299nm-Near channel analog raw data profile.
%  rawprof.pc287: float, Dim-(nbin, numFiles), 287nm-Far channel PC raw data profile.
%  rawprof.pc299: float, Dim-(nbin, numFiles), 299nm-Far channel PC raw data profile.
%  rawprof.pc287nr: float, Dim-(nbin, numFiles), 287nm-Near channel PC raw data profile.
%  rawprof.pc299nr: float, Dim-(nbin, numFiles), 299nm-Near channel PC raw data profile.
% 
%  intermprof: structure, 12 fields, 8 time-averaged profiles and
%  4 dead time corr PC signal profiles (after time-averaging)
%  sigprof: structure, 8 fields, the background subtracted profiles
%  DateTime_avg: datetime, Dim-(1, floor(numFiles/nAvg)), datetime of averaged profiles
%  TimeInHour_avg: float, Dim-(1, floor(numFiles/nAvg)), time in hour of averaged profiles


%%
function [OLfileName,save_path]=read_OL_profiles4(folder_path,save_path,nbin,dzm,bgbins,td,nAvg,sib_hkm,sib_hkm_nr)
% 
%
% Get a list of all txt files in the current folder, or subfolders of it.
% Path of the folder that stores all the profile data of a selected date(.txt)

fds = fileDatastore(folder_path, 'ReadFcn', @importdata,'FileExtensions',{'.txt'});
fullFileNames = fds.Files;
numFiles = length(fullFileNames);

% profile raw data
%% AD signal
profile_287_an_raw=nan(nbin,numFiles);
profile_299_an_raw=nan(nbin,numFiles);
profile_287_nr_an_raw=nan(nbin,numFiles);
profile_299_nr_an_raw=nan(nbin,numFiles);

%% PC signal
profile_287_pc_raw=nan(nbin,numFiles);
profile_299_pc_raw=nan(nbin,numFiles);
profile_287_nr_pc_raw=nan(nbin,numFiles);
profile_299_nr_pc_raw=nan(nbin,numFiles);

% height in m
height=dzm*[1:1:nbin]';
hkm=height/1000;
%bgbins=100;
% td=280;
% datetime array (local time)
DateTime=NaT(1,numFiles);
% time in hour (local time)
TimeInHour=nan(1,numFiles);
%% Loop over all files to load all raw signal
for k = 1 : numFiles
    %fprintf('Now reading file %s\n', fullFileNames{k});
    % read both data and metadata from the filestore
% AD signal
    Data=read(fds);
%     if k == 1
%         data_str=Data.textdata(2,1);
%         expression = ['"Far (?<f>\d*\.\d*)\D+(?<n>\d*\.\d*)\D."'];
%         token= regexp(data_str,expression,'names');
%         start_bin_fr=ceil(str2num(token{1}.f)*150/3.75)+6;
%         start_bin_nr=ceil(str2num(token{1}.n)*150/3.75)+6;
%     end 
    profile_287_an_raw(:,k)= Data.data(:,1);
    profile_299_an_raw(:,k)= Data.data(:,3);
    profile_287_nr_an_raw(:,k)= Data.data(:,5);
    profile_299_nr_an_raw(:,k)= Data.data(:,7);
% PC signal
    profile_287_pc_raw(:,k)= Data.data(:,2);
    profile_299_pc_raw(:,k)= Data.data(:,4);
    profile_287_nr_pc_raw(:,k)= Data.data(:,6);
    profile_299_nr_pc_raw(:,k)= Data.data(:,8);   
    temp=cell2mat(Data.textdata(2,1));
    infmt='dd/MM/yyyy HH:mm:ss';
    DateTime(k)=datetime(temp(10:28),'InputFormat',infmt);
    TimeInHour(k)=hour(DateTime(k))+minute(DateTime(k))/60+second(DateTime(k))/3600;
end
% save the raw signal in the structure "raw"
rawprof.an287=profile_287_an_raw;
rawprof.an299=profile_299_an_raw;
rawprof.an287nr=profile_287_nr_an_raw;
rawprof.an299nr=profile_299_nr_an_raw;
rawprof.pc287=profile_287_pc_raw;
rawprof.pc299=profile_299_pc_raw;
rawprof.pc287nr=profile_287_nr_pc_raw;
rawprof.pc299nr=profile_299_nr_pc_raw;

%% Time averaging of the raw signal
[len_height,len_time]=size(profile_287_an_raw);
[profile_287_an_avgT,~]=MoveTimeAve(profile_287_an_raw,len_height,len_time,nAvg);
[profile_299_an_avgT,~]=MoveTimeAve(profile_299_an_raw,len_height,len_time,nAvg);
[profile_287_nr_an_avgT,~]=MoveTimeAve(profile_287_nr_an_raw,len_height,len_time,nAvg);
[profile_299_nr_an_avgT,~]=MoveTimeAve(profile_299_nr_an_raw,len_height,len_time,nAvg);

[profile_287_pc_avgT,~]=MoveTimeAve(profile_287_pc_raw,len_height,len_time,nAvg);
[profile_299_pc_avgT,~]=MoveTimeAve(profile_299_pc_raw,len_height,len_time,nAvg);
[profile_287_nr_pc_avgT,~]=MoveTimeAve(profile_287_nr_pc_raw,len_height,len_time,nAvg);
[profile_299_nr_pc_avgT,TimeIndArr]=MoveTimeAve(profile_299_nr_pc_raw,len_height,len_time,nAvg);
TimeInHour_avg=TimeInHour(TimeIndArr);
DateTime_avg=DateTime(TimeIndArr);
[~,len_time_avg]=size(profile_287_an_avgT);

%% Dead time correction of the time averaged signal
% dead time corrected signal
    profile_287_pc_corr= profile_287_pc_avgT./(1-profile_287_pc_avgT/td);
    profile_299_pc_corr= profile_299_pc_avgT./(1-profile_299_pc_avgT/td);
    profile_287_nr_pc_corr= profile_287_nr_pc_avgT./(1-profile_287_nr_pc_avgT/td);
    profile_299_nr_pc_corr= profile_299_nr_pc_avgT./(1-profile_299_nr_pc_avgT/td);
% save the intermediate level signal in intermprof

intermprof.an287avg=profile_287_an_avgT;
intermprof.an299avg=profile_299_an_avgT;
intermprof.an287nravg=profile_287_nr_an_avgT;
intermprof.an299nravg=profile_299_nr_an_avgT;
intermprof.pc287avg=profile_287_pc_avgT;
intermprof.pc299avg=profile_299_pc_avgT;
intermprof.pc287nravg=profile_287_nr_pc_avgT;
intermprof.pc299nravg=profile_299_nr_pc_avgT;
intermprof.pc287corr=profile_287_pc_corr;
intermprof.pc299corr=profile_299_pc_corr;
intermprof.pc287nrcorr=profile_287_nr_pc_corr;
intermprof.pc299nrcorr=profile_299_nr_pc_corr;
%% Background subtraction signal
bg_287_an= nan(1,len_time_avg);
bg_299_an= nan(1,len_time_avg);
bg_287_nr_an= nan(1,len_time_avg);
bg_299_nr_an= nan(1,len_time_avg);

bg_287_pc= nan(1,len_time_avg);
bg_299_pc= nan(1,len_time_avg);
bg_287_nr_pc= nan(1,len_time_avg);
bg_299_nr_pc= nan(1,len_time_avg);

% background subtracted ad
profile_287_an=nan(size(profile_287_an_avgT));% 
profile_299_an=nan(size(profile_287_an_avgT));
profile_287_nr_an=nan(size(profile_287_an_avgT));% 
profile_299_nr_an=nan(size(profile_287_an_avgT));

% background subtracted pc
profile_287_pc=nan(size(profile_287_an_avgT));% 
profile_299_pc=nan(size(profile_287_an_avgT));
profile_287_nr_pc=nan(size(profile_287_an_avgT));% 
profile_299_nr_pc=nan(size(profile_287_an_avgT));

ind_sib= (hkm>sib_hkm);
sib_p = nan(2,len_time_avg);
sib_fitpc =nan(size(profile_287_an_avgT));
profile_287_pc_sib = nan(size(profile_287_an_avgT));

sib_fitan =nan(size(profile_287_an_avgT));
profile_287_an_sib = nan(size(profile_287_an_avgT));

ind_sib_nr= (hkm>sib_hkm_nr);
sib_fitan_nr =nan(size(profile_287_nr_an_avgT));
profile_287_nr_an_sib = nan(size(profile_287_nr_an_avgT));

for k=1:len_time_avg
    bg_287_an(k)= mean(profile_287_an_avgT(nbin-bgbins:nbin,k));
    bg_299_an(k)= mean(profile_299_an_avgT(nbin-bgbins:nbin,k));
    bg_287_nr_an(k)= mean(profile_287_nr_an_avgT(nbin-bgbins:nbin,k));
    bg_299_nr_an(k)= mean(profile_299_nr_an_avgT(nbin-bgbins:nbin,k));
    
    bg_287_pc(k)= mean(profile_287_pc_corr(nbin-bgbins:nbin,k));
    bg_299_pc(k)= mean(profile_299_pc_corr(nbin-bgbins:nbin,k));
    bg_287_nr_pc(k)= mean(profile_287_nr_pc_corr(nbin-bgbins:nbin,k));
    bg_299_nr_pc(k)= mean(profile_299_nr_pc_corr(nbin-bgbins:nbin,k));

    
    profile_287_an(:,k)= profile_287_an_avgT(:,k)-mean(profile_287_an_avgT(nbin-bgbins:nbin,k));
    profile_299_an(:,k)= profile_299_an_avgT(:,k)-mean(profile_299_an_avgT(nbin-bgbins:nbin,k));
    profile_287_nr_an(:,k)= profile_287_nr_an_avgT(:,k)-mean(profile_287_nr_an_avgT(nbin-bgbins:nbin,k));
    profile_299_nr_an(:,k)= profile_299_nr_an_avgT(:,k)-mean(profile_299_nr_an_avgT(nbin-bgbins:nbin,k));
    
    profile_287_pc(:,k)= profile_287_pc_corr(:,k)-mean(profile_287_pc_corr(nbin-bgbins:nbin,k));
    profile_299_pc(:,k)= profile_299_pc_corr(:,k)-mean(profile_299_pc_corr(nbin-bgbins:nbin,k));
    profile_287_nr_pc(:,k)= profile_287_nr_pc_corr(:,k)-mean(profile_287_nr_pc_corr(nbin-bgbins:nbin,k));
    profile_299_nr_pc(:,k)= profile_299_nr_pc_corr(:,k)-mean(profile_299_nr_pc_corr(nbin-bgbins:nbin,k));

    sib_p(:,k) = polyfit(hkm(ind_sib),log(movmean(profile_287_pc_corr(ind_sib,k),10)),1);
    sib_fitpc(:,k)=exp(polyval(sib_p(:,k),hkm));
    profile_287_pc_sib(:,k)= profile_287_pc_corr(:,k)-sib_fitpc(:,k);
    
    sib_p(:,k) = polyfit(hkm(ind_sib),log(movmean(profile_287_an_avgT(ind_sib,k),10)),1);
    sib_fitan(:,k)=exp(polyval(sib_p(:,k),hkm));
    profile_287_an_sib(:,k)= profile_287_an_avgT(:,k)-sib_fitan(:,k);
    
    sib_p(:,k) = polyfit(hkm(ind_sib_nr),log(movmean(profile_287_nr_an_avgT(ind_sib_nr,k),10)),1);
    sib_fitan_nr(:,k)=exp(polyval(sib_p(:,k),hkm));
    profile_287_nr_an_sib(:,k)= profile_287_nr_an_avgT(:,k)-sib_fitan_nr(:,k);
    
 
end

%% Test the SIB
figure
plot(hkm,profile_287_nr_an_avgT(:,k),'b--','LineWidth',1);hold on
plot(hkm, sib_fitan_nr(:,k),'g-','LineWidth',1); 
plot(hkm,profile_287_nr_an_sib(:,k),'k-','LineWidth',0.7);
plot(hkm,profile_287_nr_an(:,k),'m-','LineWidth',0.7)
set(gca,'YScale','log')
xlabel('Altitude (km)');
ylabel('PC signal (MHz)');
legend('287 nr an-not bg subtr','exponential fit bg','287 nr an-exp bg subtr','287nr-const bg subtr')
title(datestr(DateTime_avg(k),'yyyymmdd HH:MM:ss'))
grid on;
ylim([1e-4,1e4]);xlim([1,30]);

%%
% save the background subtracted signal in the structure "sigprof"
sigprof.an287=profile_287_an;
sigprof.an299=profile_299_an;
sigprof.an287nr=profile_287_nr_an;
sigprof.an299nr=profile_299_nr_an;
sigprof.pc287=profile_287_pc;
sigprof.pc299=profile_299_pc;
sigprof.pc287nr=profile_287_nr_pc;
sigprof.pc299nr=profile_299_nr_pc;
sigprof.pc287sib=profile_287_pc_sib;
sigprof.pcsibfit=sib_fitpc;
sigprof.an287sib=profile_287_an_sib;
sigprof.ansibfit=sib_fitan;
sigprof.an287sibnr=profile_287_nr_an_sib;
sigprof.ansibfitnr=sib_fitan_nr;

sigprof.bgan287=bg_287_an;
sigprof.bgan299=bg_299_an;
sigprof.bgan287nr=bg_287_nr_an;
sigprof.bgan299nr=bg_299_nr_an;
sigprof.bgpc287=bg_287_pc;
sigprof.bgpc299=bg_299_pc;
sigprof.bgpc287nr=bg_287_nr_pc;
sigprof.bgpc299nr=bg_299_nr_pc;


OLfileName = [save_path,'ol',datestr(DateTime(1),'yymmdd'),'.mat'];
save(OLfileName,'rawprof','sigprof','intermprof',...
    'DateTime','TimeInHour','DateTime_avg','TimeInHour_avg',...
    'height','hkm','folder_path','save_path','nbin','dzm','bgbins','td','nAvg')