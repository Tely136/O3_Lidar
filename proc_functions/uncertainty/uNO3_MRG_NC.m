% Uncertainty in ozone number density from combined ozone profiles based on
% weighting coefficients w
% This is the case when there is no correlation in vertical dimension,
% which is the case for uncertainty due to detection noise
%   Equation 95
function u = uNO3_MRG_NC(UNO3_up,uNO3_lo,w)
    u = sqrt((w.* uNO3_lo).^2 + ((1-w).*UNO3_up).^2);
end