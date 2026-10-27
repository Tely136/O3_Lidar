function [total,merged] = merge_pc_uncertainty(low,high,wLow)
    % Eqs. 95/97: merge components before combining the subtotal.
    % wLow must match the retrieval's effective LOW-range weights.

    signal_size = size(low.DET);
    wLow = wLow + zeros(signal_size);

    if ~isequal(size(wLow),signal_size) || any(wLow(:) < 0 | wLow(:) > 1)
        error('O3:MergeWeights','Weights must match profiles and lie in [0,1].');
    end

    fields = {'DET','BKG','SAT'};

    for j = 1:numel(fields)
        name = fields{j};
        low_component = low.(name);
        high_component = high.(name);

        if ~isequal(size(low_component),signal_size) || ~isequal(size(high_component),signal_size)
            error('O3:MergeSize','Component profiles must have identical sizes.');
        end

        if strcmp(name,'DET')
            uncertainty = uNO3_MRG_NC(high_component,low_component,wLow);
        else
            uncertainty = uNO3_MRG_FC(high_component,low_component,wLow);
        end

        uncertainty(wLow == 1) = low_component(wLow == 1);
        uncertainty(wLow == 0) = high_component(wLow == 0);
        merged.(name) = uncertainty;
    end

    total = sqrt(merged.DET.^2 + merged.BKG.^2 + merged.SAT.^2);
end
