% Uncertainty in ozone number density due to uncertainty in air number density profile
% This is the case when Na is derived from temperature and pressure profiles
% and there is no correlation between temperature and pressure profiles
%   Equation 86 
function u = uNO3_Na2_NC(dS_O3,dS_M,dS_O2,q_O2,Na,pa,Ta,u_pa,u_Ta)
    u = abs((dS_M + q_O2* dS_O2) ./ dS_O3) .* Na .* sqrt((u_pa.^2 ./ pa.^2) + (u_Ta.^2 ./ Ta.^2));
end