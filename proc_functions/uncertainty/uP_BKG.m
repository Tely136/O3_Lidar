function u = uP_BKG(ub0,ub1,rb0b1,z)
% Eqs.40-42: B=b0+b1*z. ub0=INTERCEPT SE; ub1=SLOPE SE.
% bg_correction returns [slope; intercept], so reverse its SE order here.
% z altitude column; parameters scalars or row vectors of profiles.
    z=z(:);
    validateattributes(ub0,{'numeric'},{'real','nonnegative'});
    validateattributes(ub1,{'numeric'},{'real','nonnegative'});
    validateattributes(rb0b1,{'numeric'},{'real','>=',-1,'<=',1});
    v=ub0.^2+ub1.^2.*z.^2+2.*z.*ub0.*ub1.*rb0b1;
    u=sqrt(max(v,0));
    u(~isfinite(v))=NaN;
end
