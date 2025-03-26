%% Convert number density to mixing ratio by volume using the temperature profile
% Down Sampling using linear interpolation 
% Input: 
%        N_O3      double  1D vector Number Density of ozone (molecules/m^3)
%        Temp      double  1D Temperature profile (K)
%        Pres0     double  surface pressure (Pa), if Pressure profile is
%                  unavailable, hydrostatic equation and Prs0 is used to obtain pressure vertical profile
% Output: ppbv_O3 double  ozone mixing ratio by volume 
function [ppbv_O3] = ND2ppbv(NDO3,Temp,Pres0,h)
%% Calculate the pressure profile: using the internation barometric formula
% P=P0(1-L(h-h0)/T0)^(gM/RL)
% Constant
R=8.31446;%ideal gas constant(J/mol.K)
g=9.80665;% gravitational acceleration (m/s^2)
M=0.0289652;% molar mass of dry air 0.0289652 (kg/mol)
L=0.0065; %K/m
Na=6.02214*10^23;
T0=Temp(1);% surface temperature
h0=h(1);
Pres=Pres0.*(1-L*(h-h0)./T0).^(g*M/(R*L));%Pa

%% Calculate the mixing ratio with the temperature and pressure profile
% air number density:ND_air= n_air*Na/V = Pres*Na/R*T
% ozone mixing ratio = ND_O3 /ND_air
ND_air= Na/R.*Pres./Temp;
ppbv_O3=NDO3./ND_air*10^9;
end
