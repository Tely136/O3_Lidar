% Test a forward iterative algorithm for calculating aerosol backscatter coefficient
% Tested/Wrote by Wu on Oct.16, 2019 at CCNY
% 
% clear all;
% close all;

% load the ceilometer attenuated backscatter coefficient
load('cldata_prof.mat');  % cldata_prof=[zkm, attenuated backscatter];
rv0=cldata_prof(:,1);
attbs=cldata_prof(:,2);

 % Load the radiosonde-measured temperature and pressure profiles
 file_path3='C:\Users\Dingdong Li\Documents\MATLAB\CLDataAnalysis\test_iterative_method\test_iterative_method\';
 file_nm_radio=[file_path3 '052616_okx_12z.dat']; % '053111_OKX_12a.dat
 
 % %Function: Calculate molecular extinction at 532-nm
 mol_am_532=[];
 mol_am_532=mol_ext_sonde(file_nm_radio);  %% Z/km, ext/km, ncn3, Pmb
 hm=mol_am_532(:,1);  %% Z-km
 am=mol_am_532(:,2);  %% mol_ext
 S2=(8*pi)/3; 
 kc=(905/532)^(-4);% 905 change to 1064 for 1064nm-Lufft ceilometers
 bm=interp1(hm,am,rv0)/S2*kc;  % mol.bs at 905-nm, km-1
 k=find(rv0<hm(1));   % take care of a few points below the radiosomde lowest altitude
 bm(k)=bm(k(length(k))+2);
 am_905=bm*S2;  % mol.ext_905
 am_905(k)=am_905(k(length(k))+2);  
 % molecular two-way transmittance
 tm2=zeros(length(rv0),1)+1;      
 for mk=2:length(rv0)
    taum=trapz(rv0(1:mk),am_905(1:mk));         
    tm2(mk)=exp(-2*taum);
 end
  
 S1_aer=40; % lidar-ratio assumption
  % Any low clouds: cloud base, top and lidar-ratio
 Z1_c=20; Z2_c=21; S1_cld=18.2; 
 
 % Calculate aerosol backscatter with a forward iterative method
 [aero_para ntimes]=iterative_func2(rv0,attbs,bm,tm2,S1_aer,Z1_c,Z2_c,S1_cld); 
 
 zmax=6.0;
 figure(1); hold on;
 plot(attbs, rv0, 'r--');
 plot(aero_para(:,2),aero_para(:,1),'b--'); hold on;
 plot(bm,rv0,'k--'); hold on;
 legend('att.bs','aer.bs','mol.bs'); 
 box on; grid on;
 axis tight;
 xlim([1e-5 1e-3]);
 ylim([0 zmax]);
 xlabel('Backscatter coefficient (/km/sr)');
 ylabel('Altitude (km)');
  
 