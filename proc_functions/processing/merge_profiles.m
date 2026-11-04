function [merged,w] = merge_profiles(prof_low,prof_high,start_merge,end_merge)
    
    [n_bin, n_prof] = size(prof_low);
    
    merged = NaN(size(prof_low));
    w = NaN(n_bin,n_prof);
    for i = 1:n_prof
        low_temp = prof_low(:,i);
        high_temp = prof_high(:,i);

        if isscalar(start_merge) || isscalar(end_merge)
            merge_ind = start_merge:end_merge;
        else
            merge_ind = start_merge(i):end_merge(i);
        end
        
        w(1:merge_ind(1)-1, i) = 1;
        w(end_merge+1:end, i) = 0;
        w(merge_ind,i) = (merge_ind(end) - merge_ind) ./ (merge_ind(end) - merge_ind(1));
        merged(:,i) = mean([low_temp high_temp], 2, "omitmissing", "Weights", [w(:,i) 1-w(:,i)]);
    end
end