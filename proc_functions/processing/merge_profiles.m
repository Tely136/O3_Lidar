function [merged,w] = merge_profiles(prof_low,prof_high,start_merge,end_merge)
    
    n_prof = size(prof_low,2);
    
    merged = NaN(size(prof_low));
    for i = 1:n_prof
        low_temp = prof_low(:,i);
        high_temp = prof_high(:,i);

        if isscalar(start_merge) || isscalar(end_merge)
            merge_ind = start_merge:end_merge;
        else
            merge_ind = start_merge(i):end_merge(i);
        end
        
        w = (merge_ind(end) - merge_ind) ./ (merge_ind(end) - merge_ind(1));
        avgs = mean([low_temp(merge_ind) high_temp(merge_ind)], 2, "Weights", [w' 1-w']);
    
        merged(:,i) = low_temp; % initialize with low profile

        if isscalar(end_merge)
            merged(end_merge:end,i) = high_temp(end_merge:end); % add high profile above merge range
            merged(merge_ind,i) = avgs; % add weighted average to merged range
        else
            merged(end_merge(i):end,i) = high_temp(end_merge(i):end); % add high profile above merge range
            merged(merge_ind,i) = avgs; % add weighted average to merged range
        end
    end

    % figure;
    % plot(low_temp, 'DisplayName','low')
    % hold on
    % plot(high_temp, 'DisplayName','high')
    % plot(merged(:,i), 'DisplayName','merge')
    % hold off
    % legend
end