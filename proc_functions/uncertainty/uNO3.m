function [u,u_DET,u_BKG,u_SAT] = uNO3(P_on,P_off,R_on,R_off,dz,L,uTau,dsigma_o3,C,uB_on,uB_off,rb_on,rb_off,z,options)
    % Photon-count subtotal: detection, background and saturation.
    % Counts, single-path cross sections; uB order is [slope; intercept].
    % options can supply signal uncertainties propagated before averaging.

    if nargin < 15
        options = struct();
    end

    if ~isfield(options,'nAverage')
        options.nAverage = 1;
    end

    if ~isfield(options,'channelCorrelation')
        options.channelCorrelation = 'independent';
    end

    n_profiles = size(P_on,2);

    if isvector(P_on)
        P_on = P_on(:);
        P_off = P_off(:);
        R_on = R_on(:);
        R_off = R_off(:);
        n_profiles = 1;
    end

    if ~isequal(size(uB_on),[2 n_profiles]) || ~isequal(size(uB_off),[2 n_profiles])
        error('O3:BackgroundSize','uB must be 2 x profiles, [slope; intercept].');
    end

    background_on = uP_BKG(uB_on(2,:),uB_on(1,:),reshape(rb_on,1,[]),z);
    background_off = uP_BKG(uB_off(2,:),uB_off(1,:),reshape(rb_off,1,[]),z);

    if isfield(options,'uDET_on') && isfield(options,'uDET_off')
        detection_on = options.uDET_on;
        detection_off = options.uDET_off;
    else
        detection_on = uP_DET(P_on,R_on,options.nAverage);
        detection_off = uP_DET(P_off,R_off,options.nAverage);
    end

    if isscalar(uTau)
        tau_on = uTau;
        tau_off = uTau;
    elseif numel(uTau) == 2
        tau_on = uTau(1);
        tau_off = uTau(2);
    else
        error('O3:DeadTimeSize','uTau scalar or [on off] seconds required.');
    end

    if isfield(options,'uSAT_on') && isfield(options,'uSAT_off')
        saturation_on = options.uSAT_on;
        saturation_off = options.uSAT_off;
    else
        saturation_on = uP_SAT(P_on,dz,L,tau_on);
        saturation_off = uP_SAT(P_off,dz,L,tau_off);
    end

    u_DET = uNO3_P_X_NC(P_on,P_off,detection_on,detection_off,dsigma_o3,dz,C,'none');
    mode = validatestring(options.channelCorrelation,{'independent','full'});

    if strcmp(mode,'independent')
        u_BKG = uNO3_P_X_NC(P_on,P_off,background_on,background_off,dsigma_o3,dz,C,'full');
        u_SAT = uNO3_P_X_NC(P_on,P_off,saturation_on,saturation_off,dsigma_o3,dz,C,'full');
    else
        u_BKG = uNO3_P_X_FC(P_on,P_off,background_on,background_off,dsigma_o3,dz,C);
        u_SAT = uNO3_P_X_FC(P_on,P_off,saturation_on,saturation_off,dsigma_o3,dz,C);
    end

    u = sqrt(u_DET.^2 + u_BKG.^2 + u_SAT.^2);
end
