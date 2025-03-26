%% Quality control of the pblh
%function pbltfp = qt_pbltop(Height,TimeInHour,raw,pbltop,cbh)
function pbltfp = newlv2qt_pbltop(TimeInHour,pbltop,varargin)

% parse
p = inputParser;
% Default setting



midnight=5;
Sunrise=12;
Midday=17;
Sunset=23;
%PBLH climatology parameter
Tnight=1.2;
Tmorningtrans=1.5;
Tnoon=2.4;
Teveningtrans=2;


addParameter(p,'midnight',midnight,@isnumeric)
addParameter(p,'sunrise',Sunrise,@isnumeric)
addParameter(p,'noon',Midday,@isnumeric)
addParameter(p,'sunset',Sunset,@isnumeric)
addParameter(p,'ThresholdNight',Tnight,@isnumeric)
addParameter(p,'ThresholdMT',Tmorningtrans,@isnumeric)
addParameter(p,'ThresholdNoon',Tnoon,@isnumeric)
addParameter(p,'ThresholdET',Teveningtrans,@isnumeric)

mov1=3;T1=0.1;
mov2=12;T2=0.15;
mov3=20;T3=0.2;
addParameter(p,'movnum1',mov1,@isnumeric)
addParameter(p,'movnum2',mov2,@isnumeric)
addParameter(p,'movnum3',mov3,@isnumeric)
addParameter(p,'T1',T1,@isnumeric)
addParameter(p,'T2',T2,@isnumeric)
addParameter(p,'T3',T3,@isnumeric)

parse(p,varargin{:})

n=length(pbltop);
pbltfp=false(1,n);
midnightInd=find(TimeInHour>p.Results.midnight,1,'first')-1;
SunriseInd=find(TimeInHour>p.Results.sunrise,1,'first')-1;
MiddayInd=find(TimeInHour>p.Results.noon,1,'first')-1;
SunsetInd=find(TimeInHour>p.Results.sunset,1,'first')-1;

% 2: Check the time before 12th hour(UTC) or 8th hour(EST) should not
% higher than 0.82km

for j=midnightInd:SunriseInd
 %   if (pbltop(j)>1)
 if (pbltop(j)>p.Results.ThresholdNight)
     % pbltop(j)=NaN;
     pbltfp(j)=1;
 end
end
for j=SunriseInd:MiddayInd
    if (pbltop(j)>p.Results.ThresholdMT)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end
for j=MiddayInd:SunsetInd
    if (pbltop(j)>p.Results.ThresholdNoon)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end
for j=SunsetInd:n
    if (pbltop(j)>p.Results.ThresholdET)
       % pbltop(j)=NaN;
        pbltfp(j)=1;
    end
end
% for j=1:midnightInd
%     if (pbltop(j)>mean(pbltop(SunsetInd:end)))
%        % pbltop(j)=NaN;
%         pbltfp(j)=1;
%     end
% end

% 3: Check the continuity of the pbltop
pbltfp=checkContinuity(pbltop,pbltfp,1,midnightInd,p.Results.movnum1,p.Results.T1);
pbltfp=checkContinuity(pbltop,pbltfp,midnightInd,SunriseInd,p.Results.movnum2,p.Results.T2);
pbltfp=checkContinuity(pbltop,pbltfp,midnightInd,SunriseInd,p.Results.movnum3,p.Results.T3);
pbltfp=checkContinuity(pbltop,pbltfp,SunriseInd,MiddayInd,p.Results.movnum1,p.Results.T1);
pbltfp=checkContinuity(pbltop,pbltfp,SunriseInd,MiddayInd,p.Results.movnum2,p.Results.T2);
pbltfp=checkContinuity(pbltop,pbltfp,MiddayInd,SunsetInd,p.Results.movnum2,p.Results.T2);
pbltfp=checkContinuity(pbltop,pbltfp,MiddayInd,SunsetInd,p.Results.movnum3,p.Results.T3);
pbltfp=checkContinuity(pbltop,pbltfp,SunsetInd,n,p.Results.movnum2,p.Results.T3);


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
