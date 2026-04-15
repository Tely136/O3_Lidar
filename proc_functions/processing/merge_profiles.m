function [merged,w] = merge_profiles(prof_low,prof_high,k,start_merge,end_merge)
    merge_ind = find(k > start_merge & k < end_merge);
    h_merge = k(merge_ind);
    
    w = (h_merge - h_merge(1)) ./ (h_merge(end) - h_merge(1));
    
    n_prof = size(prof_low,2);
    avgs2 = NaN(length(w),n_prof);
    for i = 1:length(merge_ind)
        temp = cat(1,prof_high(merge_ind(i),:),prof_low(merge_ind(i),:));
        avgs2(i,:) = mean(temp,1,"omitmissing","Weights",[w(i) 1-w(i)]);
    end

    % avgs = w.*prof_high(merge_ind,:) + (1-w).*prof_low(merge_ind,:);
    merged = prof_low;
    merged(merge_ind,:) = avgs2;
    merged(k >= end_merge,:) = prof_high(k >= end_merge,:);
end