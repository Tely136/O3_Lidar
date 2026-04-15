% Uncertainty in ozone number density due to uncertainty in air number density profile
% This is the case when Na is derived from temperature and pressure profiles
% and there is full correlation between temperature and pressure profiles
%   Equation 89 
function u = uNO3_Na2_FC(dS_O3,dS_M,dS_O2,q_O2,Na,pa,Ta,u_pa,u_Ta)
    u = abs((dS_M + q_O2* dS_O2) ./ dS_O3) .* Na .* abs((u_pa ./ pa) - (u_Ta ./ Ta));
end