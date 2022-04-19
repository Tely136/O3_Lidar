%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the molecular extinction at 532-nm based on the Rayleigh theory
% Input the radiosonde data
% March 07, 2006, Yonghua WU at CCNY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input radiosonde format: Pmbar, height/m, temperature/C,  DWPT/C(no use)   RELH/%(no use) 
% Output format: Height/km, extinction(km^-1) at 532-nm, density (mol/cm3), Pmb, Tk, Tau

function mol_am_532=mol_ext_sonde(fl_nm)

 mol_am_532=[];
 lamda_1=532; %% unit: nm

 sounding_data=load(fl_nm);

 PI=3.1415926; 
 ns= (7.247249e+18)* 1013.5/288.15; %% /* unit: molecule/cm3 */
    % 7.247249e+18=Ns*Ys/Ps (S: surface, Ns=2.546899e+19 /cm3, P3=1013.25mbar, Ts=288.15K % 
 WAVE_1=lamda_1/(1e+3);  %% wavelength unit: micro, um
 ms=6432.8 + 2949810/(146- 1.0/WAVE_1/WAVE_1) + 25540/(41-1/WAVE_1/WAVE_1);
 ms=ms*(1e-8)+1.0;
 sigma= 8*PI*PI*PI*(ms*ms-1)*(ms*ms-1)/3.0/WAVE_1/WAVE_1/WAVE_1/WAVE_1/ns/ns*(1e+16);
               %%% Rayleigh scattering cross section, unit: cm2 
 sigma1= sigma* (6+3*0.035)/(6-7*0.035);  %% Corrected by depolarization ratio delta=0.035
 
 pmb=sounding_data(:,1);
 tk=sounding_data(:,3)+273.15; %% Unit: K
 high=sounding_data(:,2)/1000;  %% k ncm3= (7.247249e+18)* pmb./tk; %% /* unit: /cm3 */
 ncm3= (7.247249e+18)* pmb./tk; %% /* molecular number density, unit: mol/cm3 */
 
 molex_1=ncm3*sigma1* (1e+5);   %% unit: km-1 */ 
 
 mm=size(molex_1,1); tau_1=zeros(mm,1); tau_2=zeros(mm,1);
 for i=2:mm-1  %% Molecular optical depth
    tmp1=( molex_1(i)+ molex_1(i+1))/2*(high(i+1)-high(i));
    tau_1(i)=tau_1(i-1)+tmp1;
 end   
 
 %% Output results
 mol_am_532=[high,molex_1,ncm3,pmb,tk,tau_1];  %% Output format: Z(km), ext(km^-1), n (1/cm3), Pmbar, Tk, Tau
 
 return;
 