function [Pon,Poff,uon,uoff,dsigma] = uncertainty_inputs(Pon,Poff,uon,uoff,dsigma,dz,C)
    validateattributes(Pon,{'numeric'},{'2d','real'});
    validateattributes(Poff,{'numeric'},{'size',size(Pon),'real'});
    validateattributes(dz,{'numeric'},{'scalar','real','finite','positive'});
    if numel(C)~=size(Pon,1), error('O3:FilterSize','One filter required per altitude.'); end
    uon=expand_input(uon,size(Pon)); uoff=expand_input(uoff,size(Pon));
    dsigma=expand_input(dsigma,size(Pon));
    if any(uon(:)<0) || any(uoff(:)<0)
        error('O3:NegativeUncertainty','Signal SD inputs must be nonnegative.');
    end
end
function x = expand_input(x,sz)
    if isvector(x) && numel(x)==sz(1), x=x(:); end
    try
        x=x+zeros(sz);
    catch
        error('O3:UncertaintySize','Incompatible dimensions.');
    end
    if ~isequal(size(x),sz), error('O3:UncertaintySize','Incompatible dimensions.'); end
end
