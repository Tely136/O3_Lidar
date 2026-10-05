% Random uncertainty in Ozone number density owing to ozone cross section differential
%   Equation 47
function u = uNO3_dS_O3_R(N,dS_O3,udS_O3_on,usS_O3_off)
    u = (2 * N ./ abs(dS_O3)) * sqrt(udS_O3_on^2 + usS_O3_off^2);
end