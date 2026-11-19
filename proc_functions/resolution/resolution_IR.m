function [dzIR,response] = resolution_IR(C,dz)
    % Eqs. 41-46: delta response for smoothing, step response for derivatives.
    % C is a coefficient vector or filter_gen cell array (column 2).

    validateattributes(dz,{'numeric'},{'scalar','real','finite','positive'});

    if iscell(C)
        dzIR = NaN(size(C));
        response = cell(size(C));

        for k = 1:numel(C)
            [dzIR(k),response{k}] = resolution_IR(C{k},dz);
        end

        return;
    end

    [c,derivative] = resolution_coefficients(C);
    response = flipud(c);

    if derivative
        response = cumsum(response);
        response(end) = 0;
    end

    response = [0;response;0];
    half_max = max(response) / 2;
    crossings = find((response(1:end-1) >= half_max) ~= (response(2:end) >= half_max));
    x = crossings + (half_max - response(crossings)) ./ (response(crossings + 1) - response(crossings));
    dzIR = dz * max(1,x(end) - x(1));
end
