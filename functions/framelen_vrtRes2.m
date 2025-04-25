% Calculate th frame length and vertical resolution 
%% Input
%   hkm: altitude [km]
%   dsigma_hkm: ozone differential cross section profile [m^2/molecule]
%   h_array: [h1,h2,h3,h4..] height of the framelength segment, from
%   surface to h1, the framelength will be len1, from h1 to h2, the frame
%   length is len2,..
%   len_array: [len0, len1, len2, len3, len4 ..], the number of bins for 
%   each framelength, the size of the len_array must be 1 greater than
%   h_array, and each element of the len_array must be odd number
%% Output
% frame_len: frame_len profile to do the derivative 
% vers_eff: effective vertical resolution [km]

function [frame_len,vres_eff]=framelen_vrtRes2(hkm,h_array,len_array,dR)
num_h_array = length(h_array);
num_len_array = length(len_array);
if num_len_array~=num_h_array+1
   error('Error. \n The size of framelength array must be 1 greater than the size of the h_array ')
end

% Using covolution with the SG filter to calculate the first order
% derivative, change the verticle resolution 
% h0=0.225;h1=10;
len_h = length(hkm);
frame_len = len_array(1)*ones(len_h,1);
ind = nan(1,num_h_array);
for j = 1:num_h_array
    ind(j)=find(hkm>=h_array(j),1,'first');
    frame_len(ind(j):end) = len_array(j+1);    
end

% 
% Calculate the vertical resolution 
y3 = zeros(1,2*len_array(1)+1);
y3(len_array(1)+1)=1;
[~,g] = sgolay(2,len_array(1));
[v, ~, ~, ~]=NDACC_ResolIR2016(dR/1000, g(:,2), y3);
vres_eff= v*ones(len_h,1);
for j = 1:num_h_array
        y3 = zeros(1,2*len_array(j+1)+1);
        y3(len_array(j+1)+1)=1;
        [~,g] = sgolay(2,len_array(j+1));
        [v, ~, ~, ~]=NDACC_ResolIR2016(dR/1000, g(:,2), y3);
        vres_eff(ind(j):end) = v;
end


% figure
% plot(vres_eff,hkm,'-');
% xlabel('Effective Vertical Resolution (km)');
% ylabel('Altitude (km)');
% ylim([0,12])
% grid on;
% set(gca,'FontSize',15);