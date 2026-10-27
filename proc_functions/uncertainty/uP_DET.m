function u = uP_DET(P,R,nAverage)
    % Eqs. 30-31: P is corrected counts; R is raw counts.
    % Propagate per acquisition before averaging.

    if nargin < 3
        nAverage = 1;
    end

    u = (P ./ R).^2 .* uDET(R,nAverage);
    u(~isfinite(P) | ~isfinite(R) | R <= 0) = NaN;
end
