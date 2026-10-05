% Uncertainty owing to detection noise
% Need to confirm how this is impacted by signal averaging
    % Equations 28 and 29 in paper
function u = uDET(R)
% R should be total counts
% laser shots * photon/sec * sec/bin
    u = sqrt(R);
end