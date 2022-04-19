%% Function: read_OL_profiles()
% Read ozone lidar signal profiles from the folder which contains a list of
% Input: 
% folder_path: string 
% save_path: string
% nbins: int 1x1, number of bins of each profile
% dz: float 1x1 vertical resolution
% bgbins: int 1x1, number of bins at the end of the profile used as background
% td: float 1x1 system dead time, the dead time correction parameter for photocounting signal
% nAvg: int 1x1 number of bins to perform time average 

% Output: netcdf file ('netcdf4_classic')
%  DateTime: datetime, Dim-(1, numFiles), datetime of each profiles
%  TimeInHour: float, Dim-(1, numFiles), time in hour of each profiles
%  Height: float, Dim- (nbin,1), height in meters
%  hkm: float, Dim- (nbin,1), height in km.
%  profile_287_an_raw: float, Dim-(nbin, numFiles), 287nm-Far channel analog raw data profile.
%  profile_299_an_raw: float, Dim-(nbin, numFiles), 299nm-Far channel analog raw data profile.
%  profile_287_nr_an_raw: float, Dim-(nbin, numFiles), 287nm-Near channel analog raw data profile.
%  profile_299_nr_an_raw: float, Dim-(nbin, numFiles), 299nm-Near channel analog raw data profile.
% 
% background subtracted data
% profile_287_an: float, Dim-(nbin, numFiles), 287nm-Far channel analog bg subtracted data profile.
% profile_299_an: float, Dim-(nbin, numFiles), 299nm-Far channel analog bg subtracted data profile.
% profile_287_nr_an: float, Dim-(nbin, numFiles), 287nm-Near channel analog bg subtracted data profile.
% profile_299_nr_an: float, Dim-(nbin, numFiles), 299nm-Near channel analog bg subtracted data profile.
% 
% PC signal
%  profile_287_pc_raw: float, Dim-(nbin, numFiles), 287nm-Far channel photocounting raw data profile.
%  profile_299_pc_raw: float, Dim-(nbin, numFiles), 299nm-Far channel photocounting raw data profile.
%  profile_287_nr_pc_raw: float, Dim-(nbin, numFiles), 287nm-Near channel photocounting raw data profile.
%  profile_299_nr_pc_raw: float, Dim-(nbin, numFiles), 299nm-Near channel photocounting raw data profile.
% 
% background subtracted data
% profile_287_pc: float, Dim-(nbin, numFiles), 287nm-Far channel photocounting bg subtracted data profile.
% profile_299_pc: float, Dim-(nbin, numFiles), 299nm-Far channel photocounting bg subtracted data profile.
% profile_287_nr_pc: float, Dim-(nbin, numFiles), 287nm-Near channel photocounting bg subtracted data profile.
% profile_299_nr_pc: float, Dim-(nbin, numFiles), 299nm-Near channel photocounting bg subtracted data profile.
% 
% dead time corrected signal
% profile_287_pc_corr: float, Dim-(nbin, numFiles), 287nm-Far channel photocounting dead time corrected data profile.
% profile_299_pc_corr: float, Dim-(nbin, numFiles), 299nm-Far channel photocounting dead time corrected data profile.
% profile_287_nr_pc_corr: float, Dim-(nbin, numFiles), 287nm-Near channel photocounting dead time corrected data profile.
% profile_299_nr_pc_corr: float, Dim-(nbin, numFiles), 299nm-Near channel photocounting dead time corrected data profile.
%
% Time averaged profile:
%  DateTime_avg: datetime, Dim-(1, floor(numFiles/nAvg)), datetime of averaged profiles
%  TimeInHour_avg: float, Dim-(1, floor(numFiles/nAvg)), time in hour of averaged profiles
% profile_287_an_avgT
% profile_299_an_avgT
% profile_287_an_nr_avgT
% profile_299_an_nr_avgT
% 
% profile_287_pc_avgT
% profile_299_pc_avgT
% profile_287_pc_nr_avgT
% profile_299_pc_nr_avgT


%%
function [OLfileName,save_path]=read_OL_profiles2(folder_path,save_path,nbin,dzm,bgbins,td,nAvg)
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
for k=1:len_time_avg
    profile_287_an(:,k)= profile_287_an_avgT(:,k)-mean(profile_287_an_avgT(nbin-bgbins:nbin,k));
    profile_299_an(:,k)= profile_299_an_avgT(:,k)-mean(profile_299_an_avgT(nbin-bgbins:nbin,k));
    profile_287_nr_an(:,k)= profile_287_nr_an_avgT(:,k)-mean(profile_287_nr_an_avgT(nbin-bgbins:nbin,k));
    profile_299_nr_an(:,k)= profile_299_nr_an_avgT(:,k)-mean(profile_299_nr_an_avgT(nbin-bgbins:nbin,k));
    
    profile_287_pc(:,k)= profile_287_pc_corr(:,k)-mean(profile_287_pc_corr(nbin-bgbins:nbin,k));
    profile_299_pc(:,k)= profile_299_pc_corr(:,k)-mean(profile_299_pc_corr(nbin-bgbins:nbin,k));
    profile_287_nr_pc(:,k)= profile_287_nr_pc_corr(:,k)-mean(profile_287_nr_pc_corr(nbin-bgbins:nbin,k));
    profile_299_nr_pc(:,k)= profile_299_nr_pc_corr(:,k)-mean(profile_299_nr_pc_corr(nbin-bgbins:nbin,k));
end 
    % save the background subtracted signal in the structure "sigprof"


sigprof.an287=profile_287_an;
sigprof.an299=profile_299_an;
sigprof.an287nr=profile_287_nr_an;
sigprof.an299nr=profile_299_nr_an;
sigprof.pc287=profile_287_pc;
sigprof.pc299=profile_299_pc;
sigprof.pc287nr=profile_287_nr_pc;
sigprof.pc299nr=profile_299_nr_pc;

OLfileName = [save_path,'ol',datestr(DateTime(1),'yymmdd_HHMM'),'_',datestr(DateTime(end),'HHMM'),'.mat'];
save(OLfileName,'rawprof','sigprof','intermprof',...
    'DateTime','TimeInHour','DateTime_avg','TimeInHour_avg',...
    'height','hkm','folder_path','save_path','nbin','dzm','bgbins','td','nAvg')