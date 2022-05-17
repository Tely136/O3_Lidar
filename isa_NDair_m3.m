
function [NDAir_m3,D_molex,molex_1,molex_2] = isa_NDair_m3(hkm,p0,T0,d_sigma)
% p0 = 1013.25*1e2;
% T0 =288.15;
% Constant
R=8.31446; %ideal gas constant(J/mol.K)
g=9.80665;% gravitational acceleration (m/s^2)
M=0.0289652;% molar mass of dry air 0.0289652 (kg/mol)
Na=6.02214*10^23;% [/mol]

T=T0-6.5*hkm;
T(hkm>11)=216.65;
P=((1-6.5*hkm/T0).^5.2561).*p0;
P(hkm>11)=226.32*1e2*exp(-g/(R*216.65)*((hkm(hkm>11)-11)*1e3));
% figure
% plot(T,height/1000)
% standar air number density
NDAir_m3=Na/R.*P./T;
%% calculate the molecular number density and extinction
lamda_1=287.2; %% unit: nm
WAVE_1=lamda_1/(1e+3);  %% wavelength unit: micro, um
lamda_2=299.1; %% unit: nm
WAVE_2=lamda_2/(1e+3);  %% wavelength unit: micro, um
PI=3.1415926; 
ns= (7.247249e+18)* 1013.5/288.15; %% /* unit: molecule/cm3 */
    % 7.247249e+18=Ns*Ys/Ps (S: surface, Ns=2.546899e+19 /cm3, P3=1013.25mbar, Ts=288.15K % 
    
ms1=1.0+(1e-8)*(6432.8 + 2949810/(146- 1.0/WAVE_1/WAVE_1) + 25540/(41-1/WAVE_1/WAVE_1));
sigma1=(6+3*0.035)/(6-7*0.035)*8*PI*PI*PI*(ms1*ms1-1)*(ms1*ms1-1)/3.0/WAVE_1/WAVE_1/WAVE_1/WAVE_1/ns/ns*(1e+16);
               %%% Rayleigh scattering cross section, unit: cm2 
               %%% Corrected by depolarization ratio delta=0.035
molex_1=sigma1*1e-4*NDAir_m3;   %% unit: m^-1 
molex_2=(lamda_2/lamda_1)^(-4)* molex_1;
delta_molex=(1-(lamda_2/lamda_1)^(-4))* molex_1;
D_molex=delta_molex./d_sigma;


