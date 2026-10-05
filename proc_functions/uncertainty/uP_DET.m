% Uncertainty in photon counts due to detection noise
%  Equations 30 and 31 in paper
function u = uP_DET(P,R)
% R is raw total signal counts across all shots in lidar signal
% P is dead-time and background corrected total signal count

    u = (P./R).^2 .* sqrt(R);
end