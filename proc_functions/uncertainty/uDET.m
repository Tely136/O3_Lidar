function u = uDET(R)
    % Eqs. 28-29: R is accumulated photon counts.

    u = sqrt(R);
    u(~isfinite(R) | R < 0) = NaN;
end
