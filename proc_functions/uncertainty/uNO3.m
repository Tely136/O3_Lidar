function [uNO3_total,components] = uNO3(R,P,N,Na,uB,rb,id_on,id_off,options)
    % R and P are raw/corrected accumulated counts for each output profile.
    % N is this range's ozone density; Na is air density; L is total shots.

    % R: Summed raw photon counts
    % P: Summed dead-time and background corrected photon counts
    % N: Ozone number density
    % Na: Air number density
    % L: Total laser shots
    % uB: 
    % rb:
    % id_on:
    % id_off:
    % options: 

    z = (1:size(P,1))' .* options.dz;
    
    uTau = options.uTau + zeros(1,2);

    u_sigma_o3_on =  options.u_sigma_o3_rel .* options.sigma_o3_on;
    u_sigma_o3_off = options.u_sigma_o3_rel .* options.sigma_o3_off;

    u_sigma_m_on = options.u_sigma_m_rel .* options.sigma_m_on;
    u_sigma_m_off = options.u_sigma_m_rel .* options.sigma_m_off;

    % uB = options.uB;
    % rb = options.rb;
    % id_on = options.id_on;
    % id_off = options.id_off;
    L = options.L;

    P_on =  P(:,:,id_on);
    P_off = P(:,:,id_off);
    R_on =  R(:,:,id_on);
    R_off = R(:,:,id_off);

    uB_on = uB(:,:,id_on);
    uB_off = uB(:,:,id_off);

    rb_on = rb(:,id_on);
    rb_off = rb(:,id_off);

    % Uncertainty in photon counts due to detection noise
    uP_DET_on =  uP_DET(P_on,  R_on);
    uP_DET_off = uP_DET(P_off, R_off);

    % Uncertainty in photon counts due to saturation/dead-time correction
    uP_SAT_on =  uP_SAT(P_on,  options.dz, L, uTau(1));
    uP_SAT_off = uP_SAT(P_off, options.dz, L, uTau(2));

    % Uncertainty in photon counts due to background correction
    uP_BKG_on =  uP_BKG(uB_on(2,:),  uB_on(1,:),  rb_on',  z);
    uP_BKG_off = uP_BKG(uB_off(2,:), uB_off(1,:), rb_off', z);

    % Uncertainty in ozone number density due to detection noise
    uNO3_DET = uNO3_P_X_NC(P_on,P_off,uP_DET_on,uP_DET_off,options.dsigma,options.dz,options.C,'independent');

    % Uncertainty in ozone number density due to saturation/dead-time correction
    uNO3_SAT = uNO3_P_X_NC(P_on,P_off,uP_SAT_on,uP_SAT_off,options.dsigma,options.dz,options.C,options.channel_correlation);

    % Uncertainty in ozone number density due to background correction
    uNO3_BKG = uNO3_P_X_NC(P_on,P_off,uP_BKG_on,uP_BKG_off,options.dsigma,options.dz,options.C,options.channel_correlation);

    % TODO: check this
    % Uncertainty in ozone number density due to differential ozone cross section random uncertainty (Eq. 47).
    O3 = uNO3_dS_O3_R(N,options.dsigma,u_sigma_o3_on,u_sigma_o3_off);

    % TODO: check this
    % Rayleigh cross sections: correlated analytical model (Eq. 55).
    RAYLEIGH = uNO3_dS_M_S(Na,options.dsigma,u_sigma_m_on,u_sigma_m_off);

    % Total uncertainty in ozone number density due to detection, saturation, background, ozone cross-section, and molecular cross-section
    uNO3_total = sqrt(uNO3_DET.^2 + uNO3_BKG.^2 + uNO3_SAT.^2 + O3.^2 + RAYLEIGH.^2);

    % Add components to struct
    components = struct( ...
        uNO3_DET=uNO3_DET, ...
        uNO3_BKG=uNO3_BKG, ...
        uNO3_SAT=uNO3_SAT, ...
        O3=O3, ...
        RAYLEIGH=RAYLEIGH ...
    );
end
