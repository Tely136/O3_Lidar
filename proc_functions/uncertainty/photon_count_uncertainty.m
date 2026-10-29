function [total,components] = photon_count_uncertainty(R,P,N,Na,uB,rb,id_on,id_off,options)
    % R and P are raw/corrected accumulated counts for each output profile.
    % N is this range's ozone density; Na is air density; L is total shots.

    if isscalar(options.uTau)
        uTau = [options.uTau options.uTau];
    else
        uTau = options.uTau;
    end

    if isscalar(options.L)
        L_on = options.L;
        L_off = options.L;
    else
        L_on = options.L(:,:,id_on);
        L_off = options.L(:,:,id_off);
    end

    % Detection noise.
    detection_on = uP_DET(P(:,:,id_on),R(:,:,id_on));
    detection_off = uP_DET(P(:,:,id_off),R(:,:,id_off));

    % Dead-time uncertainty.
    saturation_on = uP_SAT(P(:,:,id_on),options.dz,L_on,uTau(1));
    saturation_off = uP_SAT(P(:,:,id_off),options.dz,L_off,uTau(2));

    % Background fit uncertainty.
    z = (1:size(P,1))' .* options.dz;
    background_on = uP_BKG(uB(2,:,id_on),uB(1,:,id_on),rb(:,id_on)',z);
    background_off = uP_BKG(uB(2,:,id_off),uB(1,:,id_off),rb(:,id_off)',z);

    [~,components.DET,components.BKG,components.SAT] = uNO3(P(:,:,id_on),P(:,:,id_off),detection_on,detection_off,background_on,background_off,saturation_on,saturation_off,options);

    % Ozone cross sections: independent wavelength errors (Eq. 47).
    u_sigma_o3_on = options.u_sigma_o3_rel .* options.sigma_o3_on;
    u_sigma_o3_off = options.u_sigma_o3_rel .* options.sigma_o3_off;
    components.O3 = uNO3_dS_O3_R(N,options.dsigma,u_sigma_o3_on,u_sigma_o3_off);

    % Rayleigh cross sections: correlated analytical model (Eq. 55).
    u_sigma_m_on = options.u_sigma_m_rel .* options.sigma_m_on;
    u_sigma_m_off = options.u_sigma_m_rel .* options.sigma_m_off;
    components.RAYLEIGH = uNO3_dS_M_S(Na,options.dsigma,u_sigma_m_on,u_sigma_m_off);

    total = sqrt(components.DET.^2 + components.BKG.^2 + components.SAT.^2 + components.O3.^2 + components.RAYLEIGH.^2);
end
