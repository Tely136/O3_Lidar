function [NC,pathname]=newgetbschm15k_longrange(instruction,startingfolder)

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
NC.Height=Height;%upto 3km
for i=1:nfiles 
%     fprintf(1, 'Now reading %s\n', filename{i});
    raw=double(ncread([pathname filename{i}],'beta_raw'))./1e6; % read raw data, a.u.   [km^2]
%     raw=raw(1:200,:);
    %m=length(NC.Height);
    UniTime=ncread([pathname filename{i}],'time');% read the UTC
    n=length(UniTime);
    DateTime=datetime(UniTime, 'ConvertFrom', 'epochtime','Epoch','1904-01-01'); %
    mon=month(DateTime);
    Time=datevec(DateTime);
    Date=datestr(DateTime);
    NC.cs{i}=ncread([pathname filename{i}],'scaling');
    NC.p_calc{i}=ncread([pathname filename{i}],'p_calc');
    NC.Date{i}=Date(1,1:11);
%     skc=int8(double(ncread([pathname filename{i}],'sci')));% the sky condition index: (0: nothing,1: rain,2: fog,3: snow,4:precipitation or particles on window)'
    TimeInHour=Time(:,4)+(Time(:,5)+Time(:,6)./60)./60;
    %% change at 8/21 2019 --- remove all pbltop when cbh is 0-3km
 pbl=double(ncread([pathname filename{i}],'pbl'))./1000;% the range of the signal in km   
 cbh=double(ncread([pathname filename{i}],'cbh'))./1000;% the cloud base height in km
 cdp=double(ncread([pathname filename{i}],'cdp'))./1000;% the cloud depth in km

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
    T=40;% Time average number, 15*T/60 min average
   [DataAveTime,TimeIndArr]=MoveTimeAve(raw,length(Height),n,T);%MoveTimeAve(DataMat,lenHeight,TimeLen,T)
    DataAve=movmean(DataAveTime,1,1);%average over height 150m
   [pbltopAveTime,~]=MoveTimeAve(pblh_qt,1,n,T);
   [cbhAveTime,~]=MoveTimeMin(cbh(1,:),1,n,T);
   [cdpAveTime,~]=MoveTimeAve(cdp(1,:),1,n,T);
%     [skcAveTime,TimeIndArr]=MoveTimeAve(skc(:,1)',1,n,T);
   TimeArr=TimeInHour(TimeIndArr);% the new time array after time average
   datetime_array=DateTime(TimeIndArr); % the new datetime array after time average
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
   NC.AveDateTime{i}=datetime_array;
   NC.pbltop{i}=pbltopAveTime;
   NC.pbl{i}=pbl;
   NC.cbh{i}=cbhAveTime;
   NC.cdp{i}=cdpAveTime;
   NC.pbltfp{i}=pbltfp;   
   NC.pblh_raw{i}=pbltop;
   NC.pblfp2{i}=pbltfp2;
end 

% % Global Attributes:
% %            title            = 'CHM15k Nimbus'
% %            source           = 'CHM180124'
% %            device_name      = 'CHM180124'
% %            serlom           = 'TUB180032'
% %            day              = 5
% %            month            = 12
% %            year             = 2018
% %            location         = 'NYC'
% %            institution      = 'CCNY'
% %            wmo_id           = 0
% %            software_version = '17.05.1 2.13 0.754 0'
% %            comment          = ''
% %            overlap_file     = 'TUB180032 (2018-07-18 15:12:33)'
% % Dimensions:
% %            time     = 5760  (UNLIMITED)
% %            range    = 1024
% %            range_hr = 32
% %            layer    = 3
% % Variables:
% %     time          
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   double
% %            Attributes:
% %                        units     = 'seconds since 1904-01-01 00:00:00.000 00:00'
% %                        long_name = 'time UTC'
% %                        axis      = 'T'
% %     range         
% %            Size:       1024x1
% %            Dimensions: range
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'distance from lidar'
% %                        axis      = 'Z'
% %     range_hr      
% %            Size:       32x1
% %            Dimensions: range_hr
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'high resolution distance from lidar'
% %     layer         
% %            Size:       3x1
% %            Dimensions: layer
% %            Datatype:   int32
% %            Attributes:
% %                        long_name = 'layer index'
% %     latitude      
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'degrees_north'
% %                        long_name = 'latitude of site'
% %     longitude     
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'degrees_east'
% %                        long_name = 'longitude of site'
% %     azimuth       
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'degree'
% %                        long_name = 'laser direction of site'
% %     zenith        
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'degree'
% %                        long_name = 'laser direction of site'
% %     altitude      
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'altitude of ceilometer
% %                                    above mean sea level'
% %     wavelength    
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'nm'
% %                        long_name = 'laser wavelength'
% %     life_time     
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int32
% %            Attributes:
% %                        units     = 'h'
% %                        long_name = 'laser life time'
% %     range_gate    
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'length of range gate, binwidth'
% %     range_gate_hr 
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'length of range gate with high resolution, binwidth'
% %     average_time  
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int32
% %            Attributes:
% %                        units     = 'ms'
% %                        long_name = 'average time per record'
% %     laser_pulses  
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int32
% %            Attributes:
% %                        long_name = 'number of laser pulses per record (lp)'
% %     error_ext     
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int32
% %            Attributes:
% %                        long_name = '32 bit service code'
% %     temp_int      
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units        = 'K'
% %                        long_name    = 'internal temperature'
% %                        scale_factor = 0.1
% %     temp_ext      
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units        = 'K'
% %                        long_name    = 'external temperature'
% %                        scale_factor = 0.1
% %     temp_det      
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units        = 'K'
% %                        long_name    = 'detector temperature'
% %                        scale_factor = 0.1
% %     temp_lom      
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units        = 'K'
% %                        long_name    = 'laser optic module temperature'
% %                        scale_factor = 0.1
% %     state_laser   
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int8
% %            Attributes:
% %                        units     = 'percent'
% %                        long_name = 'laser quality index'
% %     state_detector
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int8
% %            Attributes:
% %                        units     = 'percent'
% %                        long_name = 'quality of detector signal'
% %     state_optics  
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int8
% %            Attributes:
% %                        units     = 'percent'
% %                        long_name = 'transmission of optics'
% %     base          
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'counts'
% %                        long_name = 'baseline raw signal in photons per shot (b)'
% %     stddev        
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   single
% %            Attributes:
% %                        units     = 'counts'
% %                        long_name = 'standard deviation raw signal in photons per shot'
% %     p_calc        
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units        = 'counts'
% %                        long_name    = 'calibration pulse in photons per shot'
% %                        scale_factor = 1e-05
% %     scaling       
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   single
% %            Attributes:
% %                        long_name = 'scaling factor (cs)'
% %     beta_raw      
% %            Size:       1024x5760
% %            Dimensions: range,time
% %            Datatype:   single
% %            Attributes:
% %                        units     = ''
% %                        long_name = 'normalized range corrected signal
% %                                    ((P_raw / lp) - b) / (cs * o(r) * p_calc) * r * r
% %                                    with P_raw = sum(P_raw_hr) * range_gate_hr / range_gate'
% %     beta_raw_hr   
% %            Size:       32x5760
% %            Dimensions: range_hr,time
% %            Datatype:   single
% %            Attributes:
% %                        units     = ''
% %                        long_name = 'normalized range corrected signal (high resolution)
% %                                    ((P_raw_hr / lp) - b) / (cs * o(r) * p_calc) * r * r'
% %     nn1           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        long_name = 'nn1'
% %     nn2           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        long_name = 'nn2'
% %     nn3           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        long_name = 'nn3'
% %     pbl           
% %            Size:       3x5760
% %            Dimensions: layer,time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'aerosol layer in PBL'
% %     pbs           
% %            Size:       3x5760
% %            Dimensions: layer,time
% %            Datatype:   int8
% %            Attributes:
% %                        long_name = 'quality score for aerosol layer in PBL'
% %     tcc           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int8
% %            Attributes:
% %                        long_name = 'total cloud cover'
% %     bcc           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int8
% %            Attributes:
% %                        long_name = 'base cloud cover'
% %     sci           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int8
% %            Attributes:
% %                        long_name = 'sky condition index
% %                                    (0: nothing, 1: rain, 2: fog, 3: snow,
% %                                    4: precipitation or particles on window)'
% %     vor           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'vertical optical range'
% %     voe           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'vertical optical range error'
% %     mxd           
% %            Size:       5760x1
% %            Dimensions: time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'maximum detection height'
% %     cbh           
% %            Size:       3x5760
% %            Dimensions: layer,time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'cloud base height'
% %     cbe           
% %            Size:       3x5760
% %            Dimensions: layer,time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'cloud base height variation'
% %     cdp           
% %            Size:       3x5760
% %            Dimensions: layer,time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'cloud depth'
% %     cde           
% %            Size:       3x5760
% %            Dimensions: layer,time
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'cloud depth variation'
% %     cho           
% %            Size:       1x1
% %            Dimensions: 
% %            Datatype:   int16
% %            Attributes:
% %                        units     = 'm'
% %                        long_name = 'cloud height offset'
