function u = uNO3_dS_M_R(Na,dS_O3,udS_M_on,udS_M_off)
    % Eq. 57: AIR density and single-path OZONE cross-section difference.

    u = Na ./ abs(dS_O3) .* hypot(udS_M_on,udS_M_off);
end
