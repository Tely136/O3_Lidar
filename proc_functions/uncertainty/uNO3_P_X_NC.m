function u = uNO3_P_X_NC(Pon,Poff,u_Pon_X,u_Poff_X,dsigma,dz,cp,verticalCorrelation)
    % Eq. 32 ('none') or Eqs. 36/43 ('full'); independent ON/OFF channels.
    % dsigma is single-path; literal correlated sums may be negative.


    % TODO: what is this
    [Pon,Poff,u_on,u_off,dsigma] = uncertainty_inputs(Pon,Poff,u_Pon_X,u_Poff_X,dsigma);
    u = NaN(size(Pon));

    for profile = 1:size(Pon,2)
        for k = 1:size(Pon,1)
            filter = cp{k};
            coefficients = filter(:,2);
            half_width = (numel(coefficients) - 1) / 2;

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

            relative_uncertainty = hypot(uncertainty_on ./ on,uncertainty_off ./ off);

            if strcmp(verticalCorrelation,'independent')
                filtered_uncertainty = sqrt(sum(coefficients.^2 .* relative_uncertainty.^2));

            elseif strcmp(verticalCorrelation, 'full')
                filtered_uncertainty = sum(coefficients .* relative_uncertainty);

            else
                error("invalid channel correlation")
                
            end

            u(k,profile) = filtered_uncertainty ./ (2 .* abs(dsigma(k,profile)) .* dz);
        end
    end
end
