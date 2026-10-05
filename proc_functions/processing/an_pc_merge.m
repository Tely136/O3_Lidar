function [pc_merge,p,r2] = an_pc_merge(an,pc,window_length,min_toggle,max_toggle)
    n_prof = size(an,2);

    p = NaN(n_prof,2);
    r2 = NaN(n_prof,1);

    pc_merge = pc(:,:,1);
    for i=1:n_prof
        an_tmp = an(:,i);
        pc_tmp = pc(:,i);
    
        valid_mask = pc_tmp >= min_toggle & pc_tmp <= max_toggle & ~isnan(an_tmp) & ~isnan(pc_tmp);
        valid_idx = find(valid_mask);
    
        n_valid = length(valid_idx);
        n_run = n_valid - window_length + 1;
    
        p_run = NaN(n_run,2);
        r2_run = NaN(n_run,1);
        for j=1:n_run
            windowIdx = valid_idx(j:j+window_length-1);
            an_window = an_tmp(windowIdx);
            pc_window = pc_tmp(windowIdx);
    
            [p_temp,s_temp] = polyfit(an_window,pc_window,1);
            p_run(j,:) = p_temp;
            r2_run(j) = s_temp.rsquared;
        end
    
        [r2_best,id_best] = max(r2_run);

        if r2_best >= 0.8
            p(i,:) = p_run(id_best,:);
            r2(i) = r2_best;
        
            pc_merge(pc_tmp>max_toggle,i) = an_tmp(pc_tmp>max_toggle)*p(i,1) + p(i,2);
        else
            if isempty(r2_best)
                warning("No valid fits for profile: %d",i);
            else
                warning("Best r2 is %0.2f",r2_best);
            end
        end
    end
end