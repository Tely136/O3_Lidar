function u = uP_DET(P,R,nAverage)
% Literal Eqs.30-31: corrected P and raw R in COUNTS, same normalization.
% Compute per acquisition BEFORE nonlinear averaging, then average variances.
% nAverage=1 by default. See README for the paper's P/background ambiguity.
    if nargin<3, nAverage=1; end
    u=(P./R).^2 .* uDET(R,nAverage);
    u(~isfinite(P) | ~isfinite(R) | R<=0)=NaN;
end
