function u = uNO3_dS_M_S(Na,dS_O3,udS_M_on,udS_M_off)
    % Eq. 55: correlated Rayleigh cross sections; single-path m^2.

    u = Na .* abs((udS_M_on - udS_M_off) ./ dS_O3);
end
