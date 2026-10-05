function [B,b,uB,rb] = bg_correction(z,data,bg_start)

    sz = size(data);
    n_bins = sz(1);

    z = z(:);
    if numel(z) ~= n_bins
        error('z must have the same length as the first dimension of data.');
    end
    if bg_start < 1 || bg_start > n_bins
        error('bg_start must be within the first dimension of data.');
    end

    data_2d = reshape(data,n_bins,[]);
    n_profile = size(data_2d,2);

    B = NaN(size(data_2d));
    b = NaN(2,n_profile);
    uB = NaN(2,n_profile);
    rb = NaN(1,n_profile);

    zBg = z(bg_start:end);

    for i = 1:n_profile

        yBg = data_2d(bg_start:end,i);

        valid = isfinite(zBg) & isfinite(yBg);
        xv = zBg(valid);
        yv = yBg(valid);

        if numel(xv) < 2 || numel(unique(xv)) < 2
            continue
        end

        mdl = fitlm(xv,yv,'linear');

        beta = mdl.Coefficients.Estimate;   % [intercept; slope]
        se = mdl.Coefficients.SE;

        b(:,i) = [beta(2); beta(1)];    % [slope; intercept]
        uB(:,i) = [se(2); se(1)];

        B(:,i) = beta(1) + beta(2) * z;

        C = mdl.CoefficientCovariance;
        rb(i) = C(1,2) / sqrt(C(1,1) * C(2,2));
    end

    B = reshape(B,size(data));
    b = reshape(b,[2,sz(2:end)]);
    uB = reshape(uB,[2,sz(2:end)]);
    rb = reshape(rb,sz(2:end));

end
