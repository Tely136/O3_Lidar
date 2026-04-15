% Uncertainty in photon counts due to saturation correction
% Equations 34 and 35 in paper
function u = uP_SAT(P,dz,L,uTau)
    c = physconst("LightSpeed");
    
    u = (c / (2*dz*L)) * P.^2 * uTau;
end