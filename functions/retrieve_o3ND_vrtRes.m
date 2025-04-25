%% [N_O3_SG,vres_eff]=retrieve_o3ND_vrtRes(pon,poff,hkm,dsigma_hkm,h0,h1,len0,len1)
% Calculate the ozone number density from DIAL signal
%% Input
%   pon: signal of 287.2 nm [MHz or mV]
%   poff: signal of 299.1 nm [MHz or mV]
%   hkm: altitude [km]
%   dsigma_hkm: ozone differential cross section profile [m^2/molecule]
%   len0: the framelen for the SG differentiator from surface to h0 [km]
%   len1: the framelen for the SG differentiator from h1[km] to top 
%   h0: altitude [km], below h0 the framelen is len0
%   h1: altitude [km], above h1 the framelen is len1, between h0 and h1,
%   the framelen is the linearly increasing from len0 to len1
%% Output
% N_O3_SG: ozone number density [molecule/m^3]
% vers_eff: effective vertical resolution [km]



function [N_O3_SG,vres_eff]=retrieve_o3ND_vrtRes(pon,poff,hkm,dsigma_hkm,h0,h1,len0,len1)

[len_h,len_t]=size(pon);
ratio_P=pon./poff;
ratio_P(ratio_P<0)=nan;
% ratio_P = sgolayfilt(sgolayfilt(ratio_P,1,sm_len),1,sm_len);
Ln_ratio_P=log(ratio_P);
dR=3.75;
N_O3_SG=nan(size(ratio_P));


% Using covolution with the SG filter to calculate the first order
% derivative, change the verticle resolution 
% h0=0.225;h1=10;
ind0=find(hkm>=h0,1,'first');
ind1=find(hkm>=h1,1,'first');
%len0=11;len1=ceil(1500/3.75);
frame_len = 11*ones(len_h,1);
frame_len(ind0:ind1)=linspace(len0,len1,ind1-ind0+1);
frame_len(ind1:end)=len1;
% round to the nearest odd number
frame_len = 2*floor(frame_len/2)+1;
vres_eff=nan(size(frame_len));
% 
% Calculate the vertical resolution 
for j =1:len_h
    y3 = zeros(1,2*frame_len(j)+1);
    y3(frame_len(j)+1)=1;
    [~,g] = sgolay(2,frame_len(j));
    [vres_eff(j), ~, ~, ~]=NDACC_ResolIR2016(dR/1000, g(:,2), y3);
end 
% figure
% plot(vres_eff,hkm,'-');
% xlabel('Effective Vertical Resolution (km)');
% ylabel('Altitude (km)');
% ylim([0,12])
% grid on;
% set(gca,'FontSize',15);


hf0 = round((frame_len(1)-1)/2);
hf1 = round((frame_len(end)-1)/2);
for i=1:len_t
    for j=1+hf0:len_h-hf1
        hf = round((frame_len(j)-1)/2);
        [~,g0] = sgolay(2,frame_len(j));
        diff0= Ln_ratio_P(j-hf:j+hf,i)'*g0(:,2);
        N_O3_SG(j,i)= -(1/2)*diff0/dsigma_hkm(j)/dR;
    end
end