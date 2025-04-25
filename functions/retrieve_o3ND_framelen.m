%% [N_O3_SG,vres_eff]=retrieve_o3ND_vrtRes(pon,poff,dsigma_hkm,frame_len_prof)
% Calculate the ozone number density from DIAL signal
%% Input
%   pon: signal of 287.2 nm [MHz or mV]
%   poff: signal of 299.1 nm [MHz or mV]
%   hkm: altitude [km]
%   dsigma_hkm: ozone differential cross section profile [m^2/molecule]
%   frame_len_prof: the framelen profile for the SG differentiator 
%% Output
% N_O3_SG: ozone number density [molecule/m^3]
% vers_eff: effective vertical resolution [km]



function [N_O3_SG]=retrieve_o3ND_framelen(pon,poff,dsigma_hkm,frame_len)

[len_h,len_t]=size(pon);
ratio_P=pon./poff;
ratio_P(ratio_P<0)=nan;
% ratio_P = sgolayfilt(sgolayfilt(ratio_P,1,sm_len),1,sm_len);
Ln_ratio_P=log(ratio_P);
dR=3.75;
N_O3_SG=nan(size(ratio_P));


hf0 = round((frame_len(1)-1)/2);
hf1 = round((frame_len(end)-1)/2);
for i=1:len_t
    for j=1+hf0:len_h-hf1
        hf = round((frame_len(j)-1)/2);
        if any(isnan(Ln_ratio_P(j-hf:j+hf,i))) % if any of the element is nan
            % skip this calculation
            continue;
        end
        [~,g0] = sgolay(2,frame_len(j));
        diff0= Ln_ratio_P(j-hf:j+hf,i)'*g0(:,2);
        N_O3_SG(j,i)= -(1/2)*diff0/dsigma_hkm(j)/dR;
    end
end