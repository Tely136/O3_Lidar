function [total,merged] = merge_pc_uncertainty(low,high,w)
    % Eqs. 95/97: merge components before combining the subtotal.
    % wLow must match the retrieval's effective LOW-range weights.

    signal_size = size(low.DET);
    w = w + zeros(signal_size);

    if ~isequal(size(w),signal_size) || any(w(:) < 0 | w(:) > 1)
        error('O3:MergeWeights','Weights must match profiles and lie in [0,1].');
    end

    fields = {'DET','BKG','SAT','O3','RAYLEIGH'};

    for j = 1:numel(fields)
        name = fields{j};
        low_component = low.(name);
        high_component = high.(name);

        if ~isequal(size(low_component),signal_size) || ~isequal(size(high_component),signal_size)
            error('O3:MergeSize','Component profiles must have identical sizes.');
        end

        if strcmp(name,'DET')
            uncertainty = uNO3_MRG_NC(high_component,low_component,w);
            
        elseif strcmp(name,'O3')
            % One shared fractional error; preserve the ozone-density sign.
            uncertainty = w .* low_component + (1-w) .* high_component;

        elseif strcmp(name,'RAYLEIGH')
            % The same atmospheric correction applies to both ranges.
            uncertainty = low_component;

        else
            uncertainty = uNO3_MRG_FC(high_component,low_component,w);

        end

        uncertainty(w == 1) = low_component(w == 1);
        uncertainty(w == 0) = high_component(w == 0);
        merged.(name) = uncertainty;
    end

    total = sqrt(merged.DET.^2 + merged.BKG.^2 + merged.SAT.^2 + merged.O3.^2 + merged.RAYLEIGH.^2);
end
