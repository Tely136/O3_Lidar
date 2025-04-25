% Calculate th frame length and vertical resolution 
%% Input
%   hkm: altitude [km]
%   dsigma_hkm: ozone differential cross section profile [m^2/molecule]
%   len0: the framelen for the SG differentiator from surface to h0 [km]
%   len1: the framelen for the SG differentiator from h1[km] to top 
%   h0: altitude [km], below h0 the framelen is len0
%   h1: altitude [km], above h1 the framelen is len1, between h0 and h1,
%   the framelen is the linearly increasing from len0 to len1
%   dR: distance of each bin  [m]
%% Output
% frame_len: frame_len profile to do the derivative 
% vers_eff: effective vertical resolution [km]

function [frame_len,vres_eff] = framelen_vrtRes(hkm,h0,h1,len0,len1,dR)
    % Using covolution with the SG filter to calculate the first order
    % derivative, change the verticle resolution 
    % h0=0.225;h1=10;
    len_h = length(hkm);
    ind0 = find(hkm>=h0,1,'first');
    ind1 = find(hkm>=h1,1,'first');
    %len0=11;len1=ceil(1500/3.75);
    frame_len = len0*ones(len_h,1);
    frame_len(ind0:ind1) = linspace(len0,len1,ind1-ind0+1);
    frame_len(ind1:end) = len1;
    % round to the nearest odd number
    frame_len = 2*floor(frame_len/2)+1;
    vres_eff = nan(size(frame_len));
    % 
    % Calculate the vertical resolution 
    for j = 1:len_h
        y3 = zeros(1,2*frame_len(j)+1);
        y3(frame_len(j)+1) = 1;
        [~,g] = sgolay(2,frame_len(j));
        [vres_eff(j), ~, ~, ~] = NDACC_ResolIR2016(dR/1000, g(:,2), y3);
    end 
end

