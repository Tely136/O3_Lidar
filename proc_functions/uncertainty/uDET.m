function u = uDET(R,nAverage)
% Eqs.28-29: raw TOTAL counts, never MHz. For the mean of independent,
% equal-exposure acquisitions var(Rmean)=Rmean/nAverage. Default=1.
    if nargin<2, nAverage=1; end
    validateattributes(nAverage,{'numeric'},{'real','finite','positive'});
    u=sqrt(R./nAverage);
    u(~isfinite(R) | R<0)=NaN;
end
