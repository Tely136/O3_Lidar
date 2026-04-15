%  Uncertainty in photon counts due to background noise extraction
% Assumes background has linear form
%   B(k) = b0 + b1 * z(k);
function u = uP_BKG(ub0,ub1,rb0b1,z)
    u = sqrt(ub0^2 + ub1^2*z.^2 + 2*z*ub0*ub1*rb0b1);
end