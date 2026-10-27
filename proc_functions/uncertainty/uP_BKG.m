function u = uP_BKG(ub0,ub1,rb0b1,z)
    % Eqs. 41-42: ub0 is intercept SE; ub1 is slope SE.

    z = z(:);
    validateattributes(ub0,{'numeric'},{'real','nonnegative'});
    validateattributes(ub1,{'numeric'},{'real','nonnegative'});
    validateattributes(rb0b1,{'numeric'},{'real','>=', - 1,'<=',1});
    variance = ub0.^2 + ub1.^2 .* z.^2 + 2 .* z .* ub0 .* ub1 .* rb0b1;
    u = sqrt(max(variance,0));
    u(~isfinite(variance)) = NaN;
end
