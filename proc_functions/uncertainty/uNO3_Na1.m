% Uncertainty in ozone number density due to uncertainty in air number
% density profile
%  This is the case when Na is not derived from temperature and pressure
%   Equation 81 
function u = uNO3_Na1(dS_O3,dS_M,dS_O2,q_O2,u_Na)
    u = abs((dS_M + q_O2* dS_O2) ./ dS_O3) .* u_Na;
end