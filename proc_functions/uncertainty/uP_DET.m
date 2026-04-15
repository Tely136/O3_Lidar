% Uncertainty in photon counts due to detection noise
% Need to confirm effect of signal averaging
%  Equations 30 and 31 in paper
function u = uP_DET(P,R)
    u = (P./R).^2 .* sqrt(R);
end