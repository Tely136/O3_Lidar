function [NC,pathname]=newgetbschm15k(instruction,startingfolder)

if nargin==0  % nargin is the number of arguments
    instruction='Pick a NC file to load';
end

if ~(nargin==2)
    startingfolder= cd;
end
%Set up to look for .nc files
[filename, pathname, filterindex] = uigetfile( ...
{  '*.nc', 'WinSpec/32 Data File, version 2.5 header'; ...
   '*.*',  'All Files (*.*)', ... % Add to this to support other files (change filterindex number below)
   },...
   instruction, ...
   startingfolder,...
   'MultiSelect', 'on');

nfiles=length(filename);% the number of files is selected
NC.filenames=filename; % Create cell for filename string
Height=double(ncread([pathname filename{1}],'range'))./1000;% read range, in km.
NC.Height=Height(:);%upto 4.5km
m=length(Height);
for i=1:nfiles 
%     fprintf(1, 'Now reading %s\n', filename{i});
    raw=double(ncread([pathname filename{i}],'beta_raw'))./1e6; % read raw data, a.u.   
%     raw=raw(1:400,:);
    %m=length(NC.Height);
    UniTime=ncread([pathname filename{i}],'time');% read the UTC
    n=length(UniTime);
    DateTime=datetime(UniTime, 'ConvertFrom', 'epochtime','Epoch','1904-01-01','TimeZone','America/New_York' ); %
    mon=month(DateTime);
    Time=datevec(DateTime);
    Date=datestr(DateTime);
    NC.Date{i}=Date(1,1:11);
%     skc=int8(double(ncread([pathname filename{i}],'sci')));% the sky condition index: (0: nothing,1: rain,2: fog,3: snow,4:precipitation or particles on window)'
    TimeInHour=Time(:,4)+(Time(:,5)+Time(:,6)./60)./60;
    %% change at 8/21 2019 --- remove all pbltop when cbh is 0-3km
 pbl=double(ncread([pathname filename{i}],'pbl'))./1000;% the range of the signal in km   
 cbh=double(ncread([pathname filename{i}],'cbh'))./1000;% the cloud base height in km
    k= pbl<0;  
    pbl(k)=NaN;% change all pbl=-999 to NaN
    pbltop=pbl(1,:);% set the lowest pblh as pbltop
    k2=cbh(1,:)>0&cbh(1,:)<3;% if the cloud is between 0-3km
    pbltop(k2)=NaN;% pbltop can't be used
    
    k3= cbh(:)<0; 
    cbh(k3)=NaN; 
%% First check the pbltop and then average
    if ((mon<3)|(mon>9))% winter
        pbltfp=newqt_pbltop(TimeInHour,pbltop,cbh(1,:),'sunrise',14,'afternoon',20,...
        'meannight',0.5,'stdnight',0.1,'meansunrise',0.5,'stdsunrise',0.15,'meannoon',1.5,...
        'stdnoon',0.3,'meanafternoon',1.5,'stdafternoon',0.3,'meansunset',1.1,'stdsunset',0.2,...
        'meanevening',0.6,'stdevening',0.1); 
    else
        pbltfp=newqt_pbltop(TimeInHour,pbltop,cbh(1,:),'meannoon',1.8,'meanafternoon',2); 
    end
    pblh_qt=pbltop;
    pblh_qt(pbltfp)=nan;
% Time average and height move average
    T=4;% Time average number, 1 min average
   [DataAveTime,TimeIndArr]=MoveTimeAve(raw,m,n,T);%MoveTimeAve(DataMat,lenHeight,TimeLen,T)
    DataAve=movmean(DataAveTime,1,1);%average over height 150m
   [pbltopAveTime,~]=MoveTimeAve(pblh_qt,1,n,T);
   [cbhAveTime,~]=MoveTimeMin(cbh(1,:),1,n,T);
%     [skcAveTime,TimeIndArr]=MoveTimeAve(skc(:,1)',1,n,T);
   TimeArr=TimeInHour(TimeIndArr);% the new time array after time average
    if ((mon<3)|(mon>9))% winter
        pbltfp2=newlv2qt_pbltop(TimeArr,pbltopAveTime);
    else
        pbltfp2=newlv2qt_pbltop(TimeArr,pbltopAveTime,'ThresholdMT',1.82,'movnum2',10,'T2',0.2);
    end
    %pbltfp2=newlv2qt_pbltop(TimeArr,pbltopAveTime,'ThresholdNight',0.8,'ThresholdMT',1.2); 
    % change when over 30% of the data are not clear
%     k=(skcAveTime>0);
%     skcAveTime(k)=1;
    
    ranglen=30;% above 200m below 450m if bs >5 then is raining
    [~,len_time]=size(DataAveTime);% first dimention: height  second dimenstion: time
    skc=zeros(1,len_time);
    for ii=1:len_time
        if mean(DataAveTime(12:ranglen,ii))>5
            skc(ii)=1;
%             pbltfp2(ii)=1;
        end
    end
%%
   NC.Avedata{i}=DataAve; 
   NC.skc{i}=skc;
   NC.Time{i}=TimeInHour;
   NC.AveTime{i}=TimeArr;
   NC.AveTimeInd{i}=TimeIndArr;
   NC.pbltop{i}=pbltopAveTime;
   NC.pbl{i}=pbl;
   NC.cbh{i}=cbhAveTime;
   NC.pbltfp{i}=pbltfp;   
   NC.pblh_raw{i}=pbltop;
   NC.pblfp2{i}=pbltfp2;
 end 
