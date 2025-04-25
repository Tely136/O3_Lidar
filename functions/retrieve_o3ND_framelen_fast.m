%% [N_O3_SG,vres_eff]=retrieve_o3ND_framelen_fast(pon,poff,dsigma_hkm,len_array,h_array,hkm)
% Calculate the ozone number density from DIAL signal
%% Input
%   pon: signal of 287.2 nm [MHz or mV]
%   poff: signal of 299.1 nm [MHz or mV]
%   hkm: altitude [km]
%   dsigma_hkm: ozone differential cross section profile [m^2/molecule]
%   h_array: [h1,h2,h3,h4..] height of the framelength segment, from
%   surface to h1, the framelength will be len0, from h1 to h2, the frame
%   length is len1,..
%   len_array: [len0, len1, len2, len3, len4 ..], the number of bins for 
%   each framelength, the size of the len_array must be 1 greater than
%   h_array, and each element of the len_array must be odd number
%% Output
% N_O3_SG: ozone number density [molecule/m^3]
% vers_eff: effective vertical resolution [km]



function [N_O3_SG]=retrieve_o3ND_framelen_fast(pon,poff,dsigma_hkm,h_array,len_array,hkm)

[len_h,len_t]=size(pon);
ratio_P=pon./poff;
ratio_P(ratio_P<0)=nan;
% ratio_P = sgolayfilt(sgolayfilt(ratio_P,1,sm_len),1,sm_len);
Ln_ratio_P=log(ratio_P);
dR=3.75;
N_O3_SG=nan(size(ratio_P));

num_h_array=length(h_array);
ind = nan(1,num_h_array);
for j = 1:num_h_array
    ind(j)=find(hkm>=h_array(j),1,'first');
end


% Using covolution with the SG filter to calculate the first order
% derivative

diff0=nan(len_h,length(len_array));
for ii = 1:len_t
    for j =1:length(len_array)
        [~,g] = sgolay(2,len_array(j));
        diff0(:,j)= conv(Ln_ratio_P(:,ii), factorial(1)/(dR)^1 * g(:,2), 'same');
    end
    diff = diff0(:,1);
    for hh =1:num_h_array
        diff(ind(hh):end) = diff0(ind(hh):end,hh+1);% concate the diff0 of different segements
    end
    N_O3_SG(:,ii)=(1/2)*diff./dsigma_hkm;
end

