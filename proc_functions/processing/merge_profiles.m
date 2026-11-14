function [merged,w] = merge_profiles(prof_low,prof_high,start_merge,end_merge)
    [n_bin,n_prof] = size(prof_low);
    w = NaN(n_bin,n_prof);

    if isscalar(start_merge)
        start_merge = repmat(start_merge,1,n_prof);
    end
    if isscalar(end_merge)
        end_merge = repmat(end_merge,1,n_prof);
    end

    for i = 1:n_prof
        if isnan(start_merge(i)) || isnan(end_merge(i))
            continue
        end

        merge_ind = start_merge(i):end_merge(i);
        w(1:start_merge(i)-1,i) = 1;
        w(end_merge(i)+1:end,i) = 0;
        w(merge_ind,i) = linspace(1,0,length(merge_ind));
    end

    % Return the weights actually used when one range is missing.
    w(isnan(prof_low) & ~isnan(prof_high) & w < 1) = 0;
    w(~isnan(prof_low) & isnan(prof_high) & w > 0) = 1;

    merged = w .* prof_low + (1-w) .* prof_high;
    merged(w == 1) = prof_low(w == 1);
    merged(w == 0) = prof_high(w == 0);
end
