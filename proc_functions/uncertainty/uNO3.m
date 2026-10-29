function [u,u_DET,u_BKG,u_SAT] = uNO3(P_on,P_off,detection_on,detection_off,background_on,background_off,saturation_on,saturation_off,options)
    % Propagate supplied signal uncertainties to the ozone number density.
    % All signals/uncertainties use the same accumulated-count normalization.

    u_DET = uNO3_P_X_NC(P_on,P_off,detection_on,detection_off,options.dsigma,options.dz,options.C,'none');

    if ~isfield(options,'channelCorrelation')
        options.channelCorrelation = 'independent';
    end

    if strcmp(options.channelCorrelation,'independent')
        u_BKG = uNO3_P_X_NC(P_on,P_off,background_on,background_off,options.dsigma,options.dz,options.C,'full');
        u_SAT = uNO3_P_X_NC(P_on,P_off,saturation_on,saturation_off,options.dsigma,options.dz,options.C,'full');
    else
        u_BKG = uNO3_P_X_FC(P_on,P_off,background_on,background_off,options.dsigma,options.dz,options.C);
        u_SAT = uNO3_P_X_FC(P_on,P_off,saturation_on,saturation_off,options.dsigma,options.dz,options.C);
    end

    u = sqrt(u_DET.^2 + u_BKG.^2 + u_SAT.^2);
end
