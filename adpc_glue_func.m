% function for Gluing the LICEL A/D and Photon-Counting(PC) data
% By Yonghua Wu at CCNY 
% Last updated on Apr.19, Apr.25, 2022
% Updated by Dingdong Li Jul. 27th, 2022

function [ch1_newp regR]=adpc_glue_func(ch1_avep,ch1_avepPC,PM1,PM2,zkm,zkm0)
%% Input parameters:
% ch1_avep: A/D signals (matrix) at the channel-1
% ch1_avepPC: Photon-counting (PC) data at the channel-1 
% PM1,PM2: PC data range for the linear regression between A/D and PC data,e.g. PM1=5MHz, PM2=20MHz
% zkm: Altitude in km
% zkm0: Initial altitude for the valid signals (skip the gating-off altitude)
%% Output parameters:
% ch1_newp: the glued A/D and PC signal profiles
%         (low-range:Converted A/D to PC sig.; middle or glue-range:average PC and AD; high-range: PC data)
% regR: correlation coefficient R, regression slope (a) and intercept(b) (y=ax+b)  (for diagose only)

[mm nn]=size(ch1_avep);
ch1_newp=zeros(mm,nn);
regR=[];
for i=1:nn
    p1=ch1_avep(:,i); % A/D signal data
    p1PC=ch1_avepPC(:,i); % PC signal data
    k=find(p1PC>=PM1 & p1PC<=PM2 & zkm>zkm0); % choose the PC and AD data for the regression, e.g. 5-20 MHz
    x1=p1(k); % A/D data
    y1=p1PC(k); %PC data
    ps= polyfit(x1,y1,1);
    k1=k(1);  % low-altitude index
    k2=k(length(k)); % high-altitude index
    R=corrcoef(x1,y1);
%     if R(2)<0.8
%         %check the 
%     end
    regR=[regR;R(2),ps(1),ps(2)]; 
    % get new "glued" signal profile
    ch1_newp(1:k1,i)=p1(1:k1)*ps(1)+ps(2); % when the PC>PM2 20MHz (lower altitude), Convert A/D data to PC (using AD only)
    ch1_newp(k2:mm,i)=p1PC(k2:mm); % when the PC<PM1 (5MHz), using the PC data only
    % when the PM>PM1 & PM<PM2 (middle range), use the average between the PC and converted PC from AD data
    tmp1=ps(1)*p1(k1:k2)+ps(2);  % converted PC data from the AD data
    tmp2=p1PC(k1:k2); % PC data
    ch1_newp(k1:k2,i)=(tmp1+tmp2)/2; % average them (PC and converted PC from AD)  
%     
%     figure;
%     plot(x1,y1,'bo'); hold on;
%     xlabel('AD');
%     ylabel('PC');
%     plot(x1,tmp1,'r-');  hold on;

end

return;