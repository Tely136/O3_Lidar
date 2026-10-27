function [dzFC,gain,f] = resolution_FC(C,dz)
    % Eqs. 49-54: first crossing of normalized gain at 0.5.
    % C is a coefficient vector or filter_gen cell array (column 2).

    validateattributes(dz,{'numeric'},{'scalar','real','finite','positive'});

    if iscell(C)
        dzFC = NaN(size(C));
        gain = cell(size(C));
        f = cell(size(C));

        for k = 1:numel(C)
            [dzFC(k),gain{k},f{k}] = resolution_FC(C{k},dz);
        end

        return;
    end

    [c,derivative] = resolution_coefficients(C);
    N = (numel(c) - 1) / 2;
    f = linspace(0,0.5,max(4097,64 * numel(c) + 1))';
    gain = zeros(size(f));

    for n = 1:N
        if derivative
            gain = gain + 2 * c(N + 1 + n) .* sin(2 * pi * n .* f);
        else
            gain = gain + 2 * c(N + 1 + n) .* cos(2 * pi * n .* f);
        end
    end

    if derivative
        gain(2:end) = gain(2:end) ./ (2 * pi .* f(2:end));
        gain(1) = 1;
    else
        gain = gain + c(N + 1);
    end

    k = find(gain <= 0.5,1);

    if isempty(k)
        fc = 0.5;
    else
        fc = f(k-1) + (0.5 - gain(k-1)) * (f(k) - f(k-1)) / (gain(k) - gain(k-1));
    end

    dzFC = dz / (2 * fc);
end
