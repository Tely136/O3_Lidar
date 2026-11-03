function range = linear_pc_range(an,pc,window_length)
    [n_bin,n_prof] = size(an);

    n_run = n_bin - window_length + 1;

    range = NaN(2,n_prof);
    ids = NaN(n_run,2);
    r2 = NaN(n_run);
    for i=1:n_prof
        an_tmp = an(:,i);
        pc_tmp = pc(:,i);


        for j=1:n_run
            windowIdx = j:j+window_length-1;
            an_window = an_tmp(windowIdx);
            pc_window = pc_tmp(windowIdx);

            [~,s_temp] = polyfit(an_window,pc_window,1);
            r2(j) = s_temp.rsquared;

            ids(j,:) = [windowIdx(1) windowIdx(end)]; 
        end

        valid_id = find(r2 > 0.9);
        lowest = valid_id(1);
        highest = valid_id(end);

        range(1,i) = ids(lowest,1); % take first index from lower range
        range(2,i) = ids(highest,1); % take first index from upper range

        % if i==1
        %     break
        % end
    end


    % figure;
    % scatter(an_tmp, pc_tmp)
    % 
    % figure;
    % plot(r2)
    % xline(range(1))
    % xline(range(2))
    % yline(0.9)
end