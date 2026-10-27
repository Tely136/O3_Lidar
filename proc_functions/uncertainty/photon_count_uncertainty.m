function [total,components] = photon_count_uncertainty(R,P,B,uB,rb,id_on,id_off,C,dsigma,dz,nAverage,L,fd,uTau,channelCorrelation)
    % R is raw acquisition counts; P and B are averaged signal/background.
    % uB is [slope SE; intercept SE]; uTau is seconds; fd is MHz.

    if nargin < 15
        channelCorrelation = 'independent';
    end

    validateattributes(nAverage,{'numeric'},{'scalar','integer','positive'});

    if isscalar(uTau)
        uTau = [uTau uTau];
    end

    if numel(uTau) ~= 2
        error('O3:DeadTimeSize','uTau scalar or [on off] required.');
    end

    n_profiles = size(P,2);
    n_bins = size(P,1);

    if size(R,2) < nAverage * n_profiles
        error('O3:ExposureSize','Not enough raw acquisitions.');
    end

    channels = [id_on id_off];
    detection = NaN(n_bins,n_profiles,2);
    saturation = detection;
    raw_mean = detection;

    for j = 1:2
        channel = channels(j);

        for i = 1:n_profiles
            acquisitions = (i - 1) * nAverage + (1:nAverage);
            raw = R(:,acquisitions,channel);
            rate = scale_binary_pc(raw,dz,L);
            corrected_total = dt_correction(raw,rate,fd);
            corrected = corrected_total - B(:,i,channel);
            detection_per_acquisition = uP_DET(corrected,raw);
            detection(:,i,j) = sqrt(sum(detection_per_acquisition.^2,2)) ./ nAverage;
            saturation(:,i,j) = mean(uP_SAT(corrected,dz,L,uTau(j)),2);
            raw_mean(:,i,j) = mean(raw,2);
            invalid = any(~isfinite(rate) | rate >= fd | rate < 0,2);
            detection(invalid,i,j) = NaN;
            saturation(invalid,i,j) = NaN;
        end
    end

    options = struct('channelCorrelation',channelCorrelation, 'uDET_on',detection(:,:,1),'uDET_off',detection(:,:,2), 'uSAT_on',saturation(:,:,1),'uSAT_off',saturation(:,:,2));
    z = (1:n_bins)' .* dz;
    [total,components.DET,components.BKG,components.SAT] = uNO3(P(:,:,id_on),P(:,:,id_off),raw_mean(:,:,1),raw_mean(:,:,2), ...
        dz,L,uTau,dsigma,C,uB(:,:,id_on),uB(:,:,id_off),rb(:,id_on),rb(:,id_off),z,options);
end
