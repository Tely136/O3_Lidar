function [NDAir_m3,D_molex,molex_1,molex_2,dsigma_hkm,P_intp,T_intp] = sonde_NDair_m3_dsigma_t(sondefile,hkm,dsigma,dsigma_t,delimiter)
    sonde = readtable(sondefile,'Delimiter',delimiter,'ReadVariableNames',true);
    P = 1e2*sonde.pressure;% convert hpa to pa
    height = sonde.height;% in meter
    T = sonde.temperature + 273.15;% convert from T to K
    
    P_intp = fillmissing(interp1(height./1e3,P,hkm,'linear'),'nearest');
    T_intp = fillmissing(interp1(height./1e3,T,hkm,'linear'),'nearest');
    
    % Constant
    R = 8.31446; %ideal gas constant(J/mol.K)
    g = 9.80665;% gravitational acceleration (m/s^2)
    M = 0.0289652;% molar mass of dry air 0.0289652 (kg/mol)
    Na = 6.02214*10^23;% [/mol]
    
    % figure
    % plot(T,height/1000)
    % standar air number density
    NDAir_m3 = Na/R.*P_intp./T_intp;
    
    %% calculate the molecular number density and extinction
    lamda_1 = 287.2; %% unit: nm
    WAVE_1 = lamda_1/(1e+3);  %% wavelength unit: micro, um
    lamda_2 = 299.1; %% unit: nm
    WAVE_2 = lamda_2/(1e+3);  %% wavelength unit: micro, um
    PI = 3.1415926; 
    ns = (7.247249e+18)* 1013.5/288.15; %% /* unit: molecule/cm3 */
        % 7.247249e+18=Ns*Ys/Ps (S: surface, Ns=2.546899e+19 /cm3, P3=1013.25mbar, Ts=288.15K % 
        
    ms1 = 1.0+(1e-8)*(6432.8 + 2949810/(146- 1.0/WAVE_1/WAVE_1) + 25540/(41-1/WAVE_1/WAVE_1));
    sigma1 = (6+3*0.035)/(6-7*0.035)*8*PI*PI*PI*(ms1*ms1-1)*(ms1*ms1-1)/3.0/WAVE_1/WAVE_1/WAVE_1/WAVE_1/ns/ns*(1e+16);
                   %%% Rayleigh scattering cross section, unit: cm2 
                   %%% Corrected by depolarization ratio delta=0.035
    molex_1 = sigma1*1e-4*NDAir_m3;   %% unit: m^-1 
    molex_2 = (lamda_2/lamda_1)^(-4)* molex_1;
    delta_molex = (1-(lamda_2/lamda_1)^(-4))* molex_1;
    
    %% interp T
    dsigma_hkm = interp1(dsigma_t,dsigma,T_intp);
    D_molex = delta_molex./double(dsigma_hkm);
end

