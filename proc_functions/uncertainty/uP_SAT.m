function u = uP_SAT(P,dz,L,uTau)
    % Eqs. 34-35: P in counts, dz in meters, L in shots, uTau in seconds.

    validateattributes(dz,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(L,{'numeric'},{'real','finite','positive'});
    validateattributes(uTau,{'numeric'},{'real','nonnegative'});
    u = 299792458 ./ (2 .* dz .* L) .* P.^2 .* uTau;
    u(~isfinite(P)) = NaN;
end
