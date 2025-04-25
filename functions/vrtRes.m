% Calculate the frame length and vertical resolution 
%% Input
%   hkm: altitude [km]
%   dsigma_hkm: ozone differential cross section profile [m^2/molecule]
%   frame_len: frame_len profile to do the derivative %% Output
%% Output
% vers_eff: effective vertical resolution [km]

function vres_eff=vrtRes(hkm,frame_len,dR)
len_h = length(hkm);
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

