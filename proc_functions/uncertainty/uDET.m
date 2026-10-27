function u = uDET(R,nAverage)
    % Eqs. 28-29: raw counts; nAverage independent acquisitions.

    if nargin < 2
        nAverage = 1;
    end

    validateattributes(nAverage,{'numeric'},{'real','finite','positive'});
    u = sqrt(R ./ nAverage);
    u(~isfinite(R) | R < 0) = NaN;
end
