function u = uNO3_P_X_FC(Pon,Poff,u_Pon_X,u_Poff_X,dsigma,dz,cp)
    % Eqs. 38/45: correlated ON/OFF channels and altitude bins.
    % dsigma is single-path; preserve the signed sum as printed.

    [Pon,Poff,u_on,u_off,dsigma] = uncertainty_inputs(Pon,Poff,u_Pon_X,u_Poff_X,dsigma,dz,cp);
    u = NaN(size(Pon));

    for profile = 1:size(Pon,2)
        for k = 1:size(Pon,1)
            filter = cp{k};
            coefficients = filter(:,2);
            half_width = (numel(coefficients) - 1) / 2;

            if half_width ~= floor(half_width)
                error('O3:FilterLength','Filters must have odd length.');
            end

            if k - half_width < 1 || k + half_width > size(Pon,1)
                continue;
            end

            bins = k - half_width:k + half_width;
            on = Pon(bins,profile);
            off = Poff(bins,profile);
            uncertainty_on = u_on(bins,profile);
            uncertainty_off = u_off(bins,profile);

            if any(~isfinite(on) | ~isfinite(off) | on <= 0 | off <= 0 | ~isfinite(uncertainty_on) | ~isfinite(uncertainty_off))
                continue;
            end

            if isfinite(dsigma(k,profile)) && dsigma(k,profile) ~= 0
                u(k,profile) = sum(coefficients .* abs(uncertainty_on ./ on - uncertainty_off ./ off)) ./ (2 .* abs(dsigma(k,profile)) .* dz);
            end
        end
    end
end
