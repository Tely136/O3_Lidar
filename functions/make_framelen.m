%% make_framelen
% Make a profile of frame length for SG differentiator
%% Input
%   hkm: altitude [km]
%   len0: the framelen for the SG differentiator from surface to h0 [km]
%   len1: the framelen for the SG differentiator from h1[km] to top 
%   h0: altitude [km], below h0 the framelen is len0
%   h1: altitude [km], above h1 the framelen is len1, between h0 and h1,
%   the framelen is the linearly increasing from len0 to len1
%% Output
% frame_len: frame_len profile to do the derivative 
function frame_len=make_framelen(hkm,h0,h1,len0,len1)
len_h = length(hkm);
ind0=find(hkm>=h0,1,'first');
ind1=find(hkm>=h1,1,'first');
%len0=11;len1=ceil(1500/3.75);
frame_len = len0*ones(len_h,1);
frame_len(ind0:ind1)=linspace(len0,len1,ind1-ind0+1);
frame_len(ind1:end)=len1;
% round to the nearest odd number
frame_len = 2*floor(frame_len/2)+1;