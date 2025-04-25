% Function: Iterative (backward) method for retrieving aerosol backscatter coefficient from the O3-DIAL at 299-nm
% Reference: Eq.(10) in Kuang et al., "Differential Absorption Lidar ...," IEEE Transactions on Geoscience and Remote Sensing, 49(1),557-571, 2011
% Last updated on June 14-15, 2022 by Yonghua Wu at CCNY

function [aero_para]=aer_iterative_func(zkm_o3,lidar_p,bsmol,o3n0,Sa,Sm,Rc,Z2,sigmaOff)
% Input: 
% zkm_o3: Ozone-lidar signal altitude in km
% lidar_p: Ozone-lidar signal intensity at 299-nm (mV or photon-count)
% bsmol: Molecular backscatter coefficientat (/1km/sr) at 299-nm
% o3n0: Ozone number density (molecule/cm3)
% Sm: Molecular extinction-to-backscatter-ratio(sr): 8*pi/3
% Sa: Aerosol extinction-to-backscatter ratio (sr)
% Rc: Reference value for the total backscatter to molecuar backsatter ratio (Rc=1 for pure molecular air or aerosol-free air)
% Z2: maximum altitude in km where the Rc is assigend
% sigmaOff: Ozone cross section at 299 nm; (cm2/molecule)
% Return:
% aero_para:[altidue/km, aer-ext/km/sr, aer-bs/km, molbs/km/sr, molext/km, o3n/cm3}

% Take the data with the altitude<Z2
k=find(zkm_o3<=Z2);
zkm2=zkm_o3(k);
lid_p=lidar_p(k);
bsmol2=bsmol(k);
o3n=o3n0(k);
molext2=bsmol2*Sm;

% sigmaOff=0.44054e-18; % cm2/molecule at 299 nm; (o3n: 1/cm3)
O3ext=(1e+5)*sigmaOff.*o3n;  % ozone extinction array (unit: 1/km)
deltaz=zkm2(2)-zkm2(1); % range-interval in km

% Assin the reference value of aerosol backscatter (extinction) at the max altitude
mn=length(lid_p);
bsa1=zeros(mn,1); % aerosol backscatter array
alpha1=zeros(mn,1); % aerosol etinction array
bsa1(mn)=(Rc-1)*bsmol2(mn); % Aerosol backscatter reference value
alpha1(mn)=bsa1(mn)*Sa; % Aerosol extinction reference value

%Backward iterative method for calculting aerosol backscatter coefficient
% % Iterative-initial calculation 
for i=(mn-1):(-1):1  % loop from max.altitude to min.altitude
    
    a1=log(lid_p(i)/lid_p(i+1));  % 1st-item in Eq.(10)
    
    b1=(o3n(i)+o3n(i+1))/2*sigmaOff*(1e+5)*deltaz;  %2nd-item in Eq.(10)
    
    c1=(molext2(i)+molext2(i+1))/2; % molecular extinction
    d1=alpha1(i+1); % Aerosol extinction
    e1=(c1+d1)*deltaz; %(3rd-item in Eq.(10))
    
    f1=(zkm2(i)/zkm2(i+1))^2;
    g1=bsmol2(i+1);
%     g1=(bsmol2(i)+bsmol2(i+1))/2; %Molecular backscatter
    h1=bsa1(i+1); %Aerosol backscatter
    m1=f1*(g1+h1); %4th item

    % all items for calculating aerosol backscatter coefficient
    bsa1(i)=exp(a1-2*b1-2*e1)*m1-bsmol2(i);
    alpha1(i)=bsa1(i)*Sa;
end
aero_para1=[zkm2,alpha1,bsa1,bsmol2,bsmol2*Sm,o3n];

% Iterative-2
bsa2=zeros(mn,1);
alpha2=zeros(mn,1);
bsa2(mn)=bsa1(mn); % aerosol reference value: using the calculated value from above (iterative-1)
alpha2(mn)=alpha1(mn);
for i=(mn-1):(-1):1
    a1=log(lid_p(i)/lid_p(i+1));  % 1st-item
    b1=(o3n(i)+o3n(i+1))/2*sigmaOff*(1e+5)*deltaz;  % 2nd-item
    
    c1=(molext2(i)+molext2(i+1))/2; 
    d1=(alpha1(i)+alpha2(i+1))/2;
    e1=(c1+d1)*deltaz; % 3rd item
    
    f1=(zkm2(i)/zkm2(i+1))^2;
    g1=(bsmol2(i)+bsmol2(i+1))/2;
%     g1=bsmol2(i+1);
    h1=(bsa1(i)+bsa2(i+1))/2;  % Now: change it with the bsa1(i) previously calculated
    m1=f1*(g1+h1);

    % all items: aerosol backscatter coefficient
    bsa2(i)=exp(a1-2*b1-2*e1)*m1-bsmol2(i);
    alpha2(i)=bsa2(i)*Sa;
end

% Calculate their relative error
bs_diff=abs(bsa2-bsa1);
bs_rel=sum(bs_diff)/sum(bsa1);

aero_para2=[zkm2,alpha2,bsa2,bsmol2,bsmol2*Sm,o3n];
aero_para=aero_para2;

% figure; 
% subplot(1,2,1); hold on;
% plot(aero_para2(:,2),aero_para2(:,1),'m-'); hold on;
% plot(aero_para2(:,5),aero_para2(:,1),'g-'); hold on;
% plot(O3ext,zkm2,'r-'); hold on;
% legend('\alpha_a','\alpha_m','O3-abs');
% ylabel('Altitude (km)');
% xlabel('Extinction (\km)');
% box on;
% subplot(1,2,2); hold on;
% plot(aero_para2(:,3),aero_para2(:,1),'m-'); hold on;
% plot(aero_para2(:,4),aero_para2(:,1),'g-'); hold on;
% legend('\beta_a','\beta_m');
% %ylabel('Altitude (km)');
% xlabel('Backscatter (\km\sr)');
% box on;

end

 