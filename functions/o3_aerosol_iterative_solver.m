%% Function: Iterative (backward) method for retrieving aerosol backscatter coefficient from the O3-DIAL at 299-nm
% Reference: Eq.(10) in Kuang et al., "Differential Absorption Lidar ...," IEEE Transactions on Geoscience and Remote Sensing, 49(1),557-571, 2011
% Last updated on June 14-15, 2022 by Yonghua Wu at CCNY, revised by D.Li

%% Input: 
% prof_299: the signal profile of 299-nm to retrieve the aerosol
% hkm: the altitude [km]
% N_O3: Ozone number density (molecule/m3)
% Sm: Molecular extinction-to-backscatter-ratio(sr): 8*pi/3
% Sa: Aerosol extinction-to-backscatter ratio (sr)
% Rc: Reference value for the total backscatter to molecuar backsatter ratio (Rc=1 for pure molecular air or aerosol-free air)
% zsurf: the lowest altitude to retrieve the aerosol profile (km)
% zref: Altitude to prescribe the reference value (km)
% am299: molecular extinction coefficient at 299 nm (m^-1)
% sigmaOff: O3 cross section at 299.1 (cm^2/molecule)
% ae: Angstrom Exponent, positive value, usually 1.5 for urban aerosol
% dsigma_hkm: ozone differential cross-section [m^2/molecule]
% frame_len: framelength to calculate the aerosol gradient
% dzm:  distance of each bin to calculate the aerosol gradient [m]
%% Return:
% abs299: aerosol backscatter coefficient at 299-nm [/sr/km]
% aext299: aerosol extinction coefficient at 299-nm [/km]
% N_absc_299: aerosol backscatter differential correction term to ozone number density [molec/m3]
% D_aext_299: aerosol extinction correction term to ozone number density [molec/m3]
% iter_num: number of iteration
% rel_NO3: relative difference of the aerosol corrected ozone number density between two adjacent iteration 

function [N_O3_corr,abs299,aext299,N_absc_299,D_aext_299,rel_NO3] = o3_aerosol_iterative_solver(prof_299,hkm,N_O3,Sm,Sa,Rc,zsurf,zref,am299,sigmaOff,ae,dsigma_hkm,frame_len,dzm)   
[len_h,len_t]=size(N_O3);
eps_O3 =0.001;% relative difference of two ozone iterations 
num_iter = 3; % number of iteration
bm_km=am299*1e3/Sm;% the molecular backscatter at 299-nm [/sr/km]
N_O3_corr = nan(size(N_O3));% Ozone number density with aerosol correction
abs299 = nan(size(N_O3));% aero backscatter 299 [/km/sr]
aext299 = nan(size(N_O3));% aero extinction 299 [/km]
N_absc_299 = nan(size(N_O3));%aero backscatter differental correction term [mole/m3]
D_aext_299 = nan(size(N_O3));%aero extinction correction term [mole/m3]
rel_NO3 = nan(num_iter,len_t);% number of iteration
abs_299_1 =nan(size(hkm)); % first estimation of aerosol backscatter

% Z1 = 4; % km
% Z2 =5; % km, K1,K2 are range window in searching for free-aerosol layer in km
for t = 1:len_t % for each profiles
    %     if cldNum(t) ==0
    %         Z1_c= 20;
    %         Z2_c=21;
    %     else
    %         Z1_c= cldBaseZ(1,t);
    %         Z2_c= cldTopZ(1,t);
    %     end
    N_O3_iter = nan(len_h,num_iter+1);% Ozone number density (molecule/m3) for iteration
            
    for i = 1: num_iter % for each time of iteration
        if i == 1
            N_O3_iter(:,1)= N_O3(:,t); % first iteration
        end
        no3_cm3 =1e-6*N_O3_iter(:,i); % O3 number density profile [molc/cm3]
        no3_cm3(hkm<zsurf)=0;% the O3 number density below zsurf set to be 0
        [aero_para] = aer_iterative_func(hkm,prof_299(:,t),bm_km,no3_cm3,Sa,Sm,Rc,zref,sigmaOff);
        %     [aero_para]=fernald_func_aercld_o3(hkm_nr,prof_299_fr_nr,bm_km,Sa,Rc,Z1,Z2,Z1_c,Z2_c,ao3_km) ;
        % if is a complex number
        img_id = (imag(aero_para)~=0);
        aero_para(img_id)=1e-5;
        % Calculate the total backscatter of 299 and 287, and the aero backscatter differental correction term
        abs_299_1(hkm<=zref) = aero_para(:,3); % aerosol backscatter coeff of 299nm [/km/sr]
        abs_287=(299.1/287.2)^ae*abs_299_1; % aerosol backscatter coeff of 287nm [/km/sr]
        am299_km= am299*1e3; % molecular extinction of 299 [/km]
        mbs_299= am299_km/Sm; % molecular backscatter coeff of 299nm [/km/sr]
        mbs_287=(299.1/287.2)^4*mbs_299; % molecular backscatter of 287nm [/km/sr]
        
        totbs_299off=abs_299_1+mbs_299; % total backscatter (aero+mole) of 299 [/km/sr]
        totbs_287on=abs_287+mbs_287; % total backscatter (aero+mole) of 287 [/km/sr]
        ratio_totbsc_onoff=totbs_287on./totbs_299off; % ratio of total backscatter (B_on/B_off)
        ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
        ln_bsc=log(ratio_totbsc_onoff); % log of ratio of total backscatter ( ln(B_on/B_off) )
        
        % Using second order sg filter to calculate the derivative of log backscatter on-off ratio
        N_absc=nan(size(ln_bsc));
        hf0 = round((frame_len(1)-1)/2);
        hf1 = round((frame_len(end)-1)/2);
        href = find(hkm<=zref,1,'last');
        for j=1+hf0:href-hf1
            hf = round((frame_len(j)-1)/2);
            [~,g0] = sgolay(2,frame_len(j));
            diff0= ln_bsc(j-hf:j+hf)'*g0(:,2);
            N_absc(j)= -(1/2)*diff0/dsigma_hkm(j)/dzm;  % backscatter item on O3 calculation
        end
        
        delta_aext=((299.1/287.2)^ae-1)*(1e-3*Sa*abs_299_1);% aerosol extinction in (/m)
        D_aext=delta_aext./double(dsigma_hkm); % extinctin item on O3 calculation
        img_id = (imag(D_aext)~=0);
        D_aext(img_id)=nan;
        img_id = (imag(N_absc)~=0);
        N_absc(img_id)=nan;
        
        % O3 correction with aerosol backscatter and extinction contribution
        N_O3_corr_temp = N_O3(:,t)- fillmissing(N_absc,'constant',0) - fillmissing(D_aext,'constant',0);
        rel_NO3_1= sum(abs(N_O3_iter(:,i)-N_O3_corr_temp),'omitnan')/sum(N_O3_iter(:,i),'omitnan');
        abs299(:,t) = abs_299_1;% aero backscatter 299 [/km/sr]
        aext299(:,t)= Sa*abs_299_1;% aero extinction 299 [/km]
        
        N_absc_299(:,t) = N_absc; %aero backscatter differental correction term [mole/m3]
        D_aext_299(:,t) = D_aext; %aero extinction correction term [mole/m3]
        rel_NO3(i,t) = rel_NO3_1;
        N_O3_corr(:,t) = N_O3_corr_temp;
        
        if rel_NO3_1 < eps_O3 % if the difference of two ozone iteration is less than eps, output the value and terminate the iteration
            break;
        else
            N_O3_iter(:,i+1) = N_O3_corr_temp; % 
        end
    end
end

