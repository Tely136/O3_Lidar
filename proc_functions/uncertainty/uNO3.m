% Combined total uncertainty in ozone number density
function [u, uNO3_P_DET_NC, uNO3_P_BKG_NC, uNO3_P_SAT_NC] = uNO3(P_on,P_off,R_on,R_off,dz,L,uTau,dsigma_o3,C,uB_on,uB_off,rb_on,rb_off,z)
    arguments
        P_on        % corrected total signal counts for on channel
        P_off       % corrected total signal counts for off channel
        R_on        % raw signal counts for on channel
        R_off       % raw signal counts for off channel
        dz          % bin size in meters
        L           % number of laser shots in raw lidar signal
        uTau       % uncertainty in dead-time constant
        dsigma_o3   % differential ozone absorption cross section
        C           % range-dependent filter coefficients
        uB_on        % background correction coefficient uncertainty for on_channel
        uB_off       % background correction coefficient uncertainty for on_channel
        rb_on          % correlation coefficient of background correction coefficients
        rb_off      % correlation coefficient of background correction coefficients
        z           % altitude grid
    end
    
    % Uncertainty in photon counts due to detection noise
    uP_DET_on = uP_DET(P_on,R_on);
    uP_DET_off = uP_DET(P_off,R_off);
    
    % Uncertainty in photon counts due to background subtraction
    uP_BKG_on = uP_BKG(uB_on(1),uB_on(2),rb_on,z)'; % fill in with real values for coefficients and their correlation
    uP_BKG_off = uP_BKG(uB_off(1),uB_off(2),rb_off,z)';
    
    % Uncertainty in photon counts due to saturation correction
    uP_SAT_on = uP_SAT(P_on,dz,L,uTau); % get value for uTau on and off
    uP_SAT_off = uP_SAT(P_off,dz,L,uTau);
    
    
    % uncertainty in ozone number density due to detection noise
    uNO3_P_DET_NC = uNO3_P_X_NC(P_on,P_off,uP_DET_on,uP_DET_off,dsigma_o3,dz,C);
    
    % Uncertainty in ozone number density from background subtraction
    uNO3_P_BKG_NC = uNO3_P_X_NC(P_on,P_off,uP_BKG_on,uP_BKG_off,dsigma_o3,dz,C);
    
    % Uncertainty in ozone number density from saturation correction
    uNO3_P_SAT_NC = uNO3_P_X_NC(P_on,P_off,uP_SAT_on,uP_SAT_off,dsigma_o3,dz,C);
    
    % Total combined uncertainty
    u2 = uNO3_P_DET_NC.^2 + uNO3_P_BKG_NC.^2 + uNO3_P_SAT_NC.^2; % plus other terms later
    u = sqrt(u2);
end