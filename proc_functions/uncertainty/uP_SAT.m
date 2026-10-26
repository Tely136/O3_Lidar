function u = uP_SAT(P,dz,L,uTau)
% Literal Eqs.34-35. P corrected COUNTS; dz meters; L shots; uTau SECONDS.
% Fixed dead time is correlated in altitude/time: average this linearly.
    validateattributes(dz,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(L,{'numeric'},{'real','finite','positive'});
    validateattributes(uTau,{'numeric'},{'real','nonnegative'});
    u=299792458 ./ (2.*dz.*L) .* P.^2 .* uTau;
    u(~isfinite(P))=NaN;
end
