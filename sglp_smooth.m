%% Vertical smoothing of data file
% input:
% start_bin: int, 1x1, the ordinal number of the bin that the Gate opens
% sgwin_len: int, 1x1, the window size of the SG filter smoothing
% lp_flag: bool, 1x1, the flag of wether using low pass filter after SG
% filtering, 1: yes, 0: no
% wpass: normalized passband frequency wpass in units of π rad/sample 
% for lowpass filtering  
% prof: float, len_height x len_time  the signal profile to perform smoothing
% output:
% 
function smprof = sglp_smooth(prof,start_bin,h1,h2,movnum1,movnum2,movnum3)

smprof=prof;
order=1;
% smprof(start_bin:end,:)=sgolayfilt(sgolayfilt(prof(start_bin:end,:),order,sgwin_len),order,sgwin_len);
temp=nan(size(prof));
temp(start_bin:end,:)=sgolayfilt(prof(start_bin:end,:),order,movnum1);
smprof(start_bin:h1,:)=temp(start_bin:h1,:);
temp(start_bin:end,:)=sgolayfilt(prof(start_bin:end,:),order,movnum2);
smprof(h1+1:h2,:)=temp(h1+1:h2,:);
temp(start_bin:end,:)=sgolayfilt(prof(start_bin:end,:),order,movnum3);
smprof(h2+1:end,:)=temp(h2+1:end,:);

smprof = movmean(smprof,movnum1,1);

end