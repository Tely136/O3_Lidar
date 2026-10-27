function [c,derivative] = resolution_coefficients(C)
    if ~isvector(C)
        C = C(:,2);
    end

    validateattributes(C,{'numeric'},{'vector','real','finite','nonempty'});
    c = C(:);
    N = (numel(c) - 1) / 2;

    if N ~= floor(N)
        error('O3:FilterLength','Filters must have odd length.');
    end

    tolerance = 1e-10 * max(abs(c));
    derivative = max(abs(c + flipud(c))) <= tolerance;

    if derivative
        scale = (-N:N) * c;
    elseif max(abs(c - flipud(c))) <= tolerance
        scale = sum(c);
    else
        error('O3:FilterSymmetry','Filters must be symmetric or antisymmetric.');
    end

    if abs(scale) <= tolerance
        error('O3:FilterNormalization','Filter normalization must be nonzero.');
    end

    c = c ./ scale;
end
