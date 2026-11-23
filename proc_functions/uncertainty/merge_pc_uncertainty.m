function [total,merged] = merge_pc_uncertainty(low,high,w)
    % w is the low-range weight returned by merge_profiles.

    w = w + zeros(size(low.uNO3_DET));

    merged.uNO3_DET = uNO3_MRG_NC(high.uNO3_DET,low.uNO3_DET,w);
    merged.uNO3_BKG = uNO3_MRG_FC(high.uNO3_BKG,low.uNO3_BKG,w);
    merged.uNO3_SAT = uNO3_MRG_FC(high.uNO3_SAT,low.uNO3_SAT,w);

    % Cross-section errors are shared between ranges.
    merged.O3 = w .* low.O3 + (1-w) .* high.O3;
    merged.RAYLEIGH = low.RAYLEIGH;

    fields = fieldnames(merged);
    for i = 1:length(fields)
        name = fields{i};
        lo = low.(name) + zeros(size(w));
        hi = high.(name) + zeros(size(w));
        u = merged.(name) + zeros(size(w));
        u(w == 1) = lo(w == 1);
        u(w == 0) = hi(w == 0);
        merged.(name) = u;
    end

    total = sqrt(merged.uNO3_DET.^2 + merged.uNO3_BKG.^2 + merged.uNO3_SAT.^2 + merged.O3.^2 + merged.RAYLEIGH.^2);
end
