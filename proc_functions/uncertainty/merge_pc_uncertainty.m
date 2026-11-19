function [total,merged] = merge_pc_uncertainty(low,high,w)
    % w is the low-range weight returned by merge_profiles.

    w = w + zeros(size(low.DET));

    merged.DET = uNO3_MRG_NC(high.DET,low.DET,w);
    merged.BKG = uNO3_MRG_FC(high.BKG,low.BKG,w);
    merged.SAT = uNO3_MRG_FC(high.SAT,low.SAT,w);

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

    total = sqrt(merged.DET.^2 + merged.BKG.^2 + merged.SAT.^2 + merged.O3.^2 + merged.RAYLEIGH.^2);
end
