function u = uP_DET(P,R)
    % Eqs. 30-31: corrected P and raw R have the same accumulated shots.

    u = (P ./ R).^2 .* uDET(R);
    u(~isfinite(P) | ~isfinite(R) | R <= 0) = NaN;
end
