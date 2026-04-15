% Uncertainty in ozone number density from combined ozone profiles based on
% weighting coefficients w
% This is the case when there is full correlation in vertical dimension,
% which is the case for uncertainty due to detection saturation and
% background correction
%   Equation 97
function u = uNO3_MRG_FC(UNO3_up,uNO3_lo,w)
    u = abs(w.* uNO3_lo + (1-w).*UNO3_up);
end