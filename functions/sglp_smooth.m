%% Vertical smoothing of data file
% input:
% prof: the prof to smooth
% h0: the alititude of the end point of the first segment to smooth
% h1: the alititude of the end point of the second segment to smooth
% h0: the alititude of the end point of the first segment to smooth
% h1: the alititude of the end point of the second segment to smooth
% prof: float, len_height x len_time  the signal profile to perform smoothing
% output:
% 
function smprof = sglp_smooth(prof,h0,h1,h2,movnum1,movnum2,movnum3)

smprof=prof;
order=1;
start_bin = h0;
temp=nan(size(prof));
temp(start_bin:end,:)=sgolayfilt(prof(start_bin:end,:),order,movnum1);
smprof(start_bin:h1,:)=temp(start_bin:h1,:);
temp(start_bin:end,:)=sgolayfilt(prof(start_bin:end,:),order,movnum2);
smprof(h1+1:h2,:)=temp(h1+1:h2,:);
temp(start_bin:end,:)=sgolayfilt(prof(start_bin:end,:),order,movnum3);
smprof(h2+1:end,:)=temp(h2+1:end,:);

smprof = movmean(smprof,movnum1,1);

end