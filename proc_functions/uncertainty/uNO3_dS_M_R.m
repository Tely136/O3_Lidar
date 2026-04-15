% Random uncertainty in Ozone number density owing to Rayleigh extinction cross section differential
%   Equation 57
function u = uN_dS_M_R(N,dS_M,udS_M_on,usS_M_off)
    u = (2 * N ./ abs(dS_M)) * sqrt(udS_M_on^2 + usS_M_off^2);
end