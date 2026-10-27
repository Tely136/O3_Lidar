function u = uNO3_dS_O3_R(N,dS_O3,udS_O3_on,udS_O3_off)
    % Eq. 47: ozone density; single-path cross sections in m^2.

    u = N ./ abs(dS_O3) .* hypot(udS_O3_on,udS_O3_off);
end
