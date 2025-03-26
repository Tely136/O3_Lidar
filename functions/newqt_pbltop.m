%% Quality control of the pblh
%function pbltfp = qt_pbltop(Height,TimeInHour,raw,pbltop,cbh)
function pbltfp = newqt_pbltop(TimeInHour,pbltop,cbh,varargin)

% parse
p = inputParser;
% Default setting

nightTime = 2;
midnightTime=5;
SunriseTime = 12;
noonTime=17;
afternoonTime=19;
SunsetTime=21;

%PBLH climatology parameter
mean_night=1;
std_night=0.3;
mean_sunrise=0.474;
std_sunrise=0.2;
mean_noon=1.5;
std_noon=0.38;
mean_afternoon=1.67;
std_afternoon=0.35;
mean_sunset=1.55;
std_sunset=0.6;
mean_evening=1;
std_evening=0.6;
mean_midnight=0.6;
std_midnight=0.3;





addParameter(p,'night',nightTime,@isnumeric)
addParameter(p,'sunrise',SunriseTime,@isnumeric)
addParameter(p,'noon',noonTime,@isnumeric)
addParameter(p,'afternoon',afternoonTime,@isnumeric)

addParameter(p,'sunset',SunsetTime,@isnumeric)

addParameter(p,'meannight',mean_night,@isnumeric)
addParameter(p,'stdnight',std_night,@isnumeric)
addParameter(p,'meansunrise',mean_sunrise,@isnumeric)
addParameter(p,'stdsunrise',std_sunrise,@isnumeric)
addParameter(p,'meannoon',mean_noon,@isnumeric)
addParameter(p,'stdnoon',std_noon,@isnumeric)
addParameter(p,'meanafternoon',mean_afternoon,@isnumeric)
addParameter(p,'stdafternoon',std_afternoon,@isnumeric)
addParameter(p,'meansunset',mean_sunset,@isnumeric)
addParameter(p,'stdsunset',std_sunset,@isnumeric)
addParameter(p,'meanevening',mean_evening,@isnumeric)
addParameter(p,'stdevening',std_evening,@isnumeric)

lowT_midnight=0.08;
lowT_evening=0.1;
lowT_morning=0.1;
lowT_noon=0.2;

addParameter(p,'lowT_midnight',lowT_midnight,@isnumeric)
addParameter(p,'lowT_evening',lowT_evening,@isnumeric)
addParameter(p,'lowT_morning',lowT_morning,@isnumeric)
addParameter(p,'lowT_noon',lowT_noon,@isnumeric)

movnum1=39;% 10 min 
movnum2=112;% 30min
movnum3=225;% 1hour
T1=0.1;
T2=0.15;
T3=0.3;

addParameter(p,'movnum1',movnum1,@isnumeric)
addParameter(p,'movnum2',movnum2,@isnumeric)
addParameter(p,'movnum3',movnum3,@isnumeric)
addParameter(p,'T1',T1,@isnumeric)
addParameter(p,'T2',T2,@isnumeric)
addParameter(p,'T3',T3,@isnumeric)

HalfCldSpan=5;
addParameter(p,'HalfCldSpan',HalfCldSpan,@isnumeric)

parse(p,varargin{:});

n=length(pbltop);
pbltfp=false(1,n);
% pbltfp2=false(1,n);
% pbltfp3=false(1,n);
% 1: Check if there is no cloud aboud the pbl but the backscatter around the detected pbl is higher than the
% ground level, then it is likely to be the loft layer instead of pbl
% ground=15:25;% ground index range
% Thresh_loft=0.02;
% for j=1:n
%     if (pbltop(j)>1)
%     if (isnan(mean(cbh(max(j-3,1):min(j+3,n)),'omitnan'))&& (~(isnan(pbltop(j)))))% if there is no low cloud around and pbl is not nan
%        k=abs(Height-pbltop(j))<0.05;
%        pbltfp1(j)=(mean(raw(k,j))-mean(raw(ground,j)))>Thresh_loft;
% %        if ((mean(raw(k,j))-mean(raw(ground,j)))>Thresh_loft)% if around detected pbl (100m) is higher than the gound level, 
% %             pbltfp(j)=1;%
% %        end
%     end
%     end
% end
% HalfCldSpan=5;
% 
for j=1:n
    if (min(cbh(max(j-p.Results.HalfCldSpan,1):min(j+p.Results.HalfCldSpan,n)),[],'omitnan')<3)% for each detected pbl if around there exist cloud lower than 3km 
        pbltfp(j)=1;
    end
end
%% 1: Based on the seasonal average, set the threshold for different time
% period.
% [0,2) H<=1km
% [2,11) H<=0.67km
% [11,20) H<=0.17x-1.04
% [20,24) H<= 7.47-0.26x

nightInd=find(TimeInHour>p.Results.night,1,'first')-1;
midnightInd=find(TimeInHour>midnightTime,1,'first')-1;
SunriseInd=find(TimeInHour>p.Results.sunrise,1,'first')-1;
middayInd=find(TimeInHour>p.Results.noon,1,'first')-1;
AfternoonInd=find(TimeInHour>p.Results.afternoon,1,'first')-1;
SunsetInd=find(TimeInHour>p.Results.sunset,1,'first')-1;


for j=1:nightInd
 if (pbltop(j)>p.Results.meannight+1.5*p.Results.stdnight||pbltop(j)<p.Results.lowT_evening)
     % pbltop(j)=NaN;
     pbltfp(j)=1;
 end

end
for j=nightInd+1:SunriseInd
    if (pbltop(j)>p.Results.meansunrise+1.5*p.Results.stdsunrise||pbltop(j)<p.Results.lowT_midnight)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
    
end
MT_slop=(p.Results.meannoon-p.Results.meansunrise)/(p.Results.noon-p.Results.sunrise);
MT_intc=MT_slop*(-p.Results.sunrise)+p.Results.meansunrise;
for j=SunriseInd+1:middayInd
    if (pbltop(j)>(TimeInHour(j)*MT_slop+MT_intc+1.5*(p.Results.stdnoon+p.Results.stdsunrise)/2)||pbltop(j)<p.Results.lowT_morning)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end
for j=middayInd+1:SunsetInd
    if (pbltop(j)>p.Results.meanafternoon+1.5*p.Results.stdafternoon||pbltop(j)<p.Results.lowT_noon)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end
% ET_slop=(mean_midnight-mean_evening)/(24-SunsetTime);
% ET_intc=ET_slop*(-SunsetTime)+mean_evening;
for j=SunsetInd:n
    % if (pbltop(j)>7.79-0.26*TimeInHour(j))
    % if (pbltop(j)>(TimeInHour(j)*ET_slop+ET_intc+1.5*std_evening)||pbltop(j)<(TimeInHour(j)*ET_slop+ET_intc-1.5*std_evening))
    if (pbltop(j)>p.Results.meansunset+1.5*p.Results.stdsunset||pbltop(j)<p.Results.lowT_evening)
        % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end
%% in the morning transition, if the earlier pbl is higher than the later pbltop is is loft layer

%% 2: Check the continuity of the pbltop
% tmp=pbltop;
% movnum1=100;% 30min
% % stdThreshold=2;
% for j=1:nightInd
%     movnumMedian=median(tmp(max(j-round((movnum1-1)/2),1):min(j+round((movnum1-1)/2),n)),'omitnan');
%     %movnumStd=std(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),0,'omitnan');
%     if (isnan(tmp(j)))
%         continue;
%     end
%    % if(abs(tmp(j)-movnumMean)>stdThreshold*movnumStd)
%     if(abs(tmp(j)-movnumMedian)>0.1)
%         pbltfp(j)=1;
%     end
% end 


% from 0-night, 30min movmedian (112), 0.1km threshold, apply 1 time
pbltfp=checkContinuity(pbltop,pbltfp,1,nightInd,p.Results.movnum2,0.1);
% from night-sunrise,30min movmedian (112), 0.1km threshold, apply 1 time  ,1hr movmedian (112), 0.1km threshold, apply 1 time 
pbltfp=checkContinuity(pbltop,pbltfp,nightInd,SunriseInd,p.Results.movnum2,0.1);
pbltfp=checkContinuity(pbltop,pbltfp,nightInd,SunriseInd,p.Results.movnum3,0.1);

% movnum=225;% one hour
% % stdThreshold=2;
% for j=nightInd+1:SunriseInd
%     movnumMedian=median(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),'omitnan');
%     %movnumStd=std(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),0,'omitnan');
%     if (isnan(tmp(j)))
%         continue;
%     end
%    % if(abs(tmp(j)-movnumMean)>stdThreshold*movnumStd)
%     if(abs(tmp(j)-movnumMedian)>0.1)
%         pbltfp(j)=1;
%     end
% end 

% from night-sunrise,10 min movmedian (39), 0.1km threshold, apply 1 time  
% ,30min movmedian (112), 0.1km threshold, apply 1 time 
pbltfp=checkContinuity(pbltop,pbltfp,SunriseInd,AfternoonInd,p.Results.movnum1,0.1);
pbltfp=checkContinuity(pbltop,pbltfp,SunriseInd,AfternoonInd,p.Results.movnum2,0.15);

% movnum=112;% half hour
% % stdThreshold=2;
% for j=SunriseInd+1:AfternoonInd
%     movnumMedian=median(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),'omitnan');
%     %movnumStd=std(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),0,'omitnan');
%     if (isnan(tmp(j)))
%         continue;
%     end
%    % if(abs(tmp(j)-movnumMean)>stdThreshold*movnumStd)
%     if(abs(tmp(j)-movnumMedian)>0.2)
%         pbltfp(j)=1;
%     end
% end 

% from AfternoonInd,30 min movmedian (112), 0.15km threshold, apply 1 time  
% ,1hr movmedian (225), 0.3km threshold, apply 1 time 
pbltfp=checkContinuity(pbltop,pbltfp,AfternoonInd,n,p.Results.movnum2,0.15);
% pbltfp=checkContinuity(pbltop,pbltfp,AfternoonInd,n,p.Results.movnum2,0.3);

% movnum=225;
% % stdThreshold=2;
% for j=AfternoonInd:n
%     movnumMedian=median(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),'omitnan');
%     %movnumStd=std(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),0,'omitnan');
%     if (isnan(tmp(j)))
%         continue;
%     end
%    % if(abs(tmp(j)-movnumMean)>stdThreshold*movnumStd)
%     if(abs(tmp(j)-movnumMedian)>0.3)
%         pbltfp(j)=1;
%     end
% end 

for j=midnightInd+1:middayInd
    if (pbltop(j)>median(pbltop(max(j,middayInd-100):middayInd),'omitnan')+0.1)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end

function pbltfp_output=checkContinuity(pblh,pbltfp,startInd,endInd,movnum,threshold)
pbltfp_output=pbltfp;
tmp=pblh;
tmp(pbltfp)=nan;
n=length(pblh);
% stdThreshold=2;
for j=startInd:endInd
    movnumMedian=median(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),'omitnan');
    %movnumStd=std(tmp(max(j-round((movnum-1)/2),1):min(j+round((movnum-1)/2),n)),0,'omitnan');
    if (isnan(tmp(j)))
        continue;
    end
   % if(abs(tmp(j)-movnumMean)>stdThreshold*movnumStd)
    if(abs(tmp(j)-movnumMedian)>threshold)
        pbltfp_output(j)=1;
    end
end 

