function [range,fit] = linear_fr_nr_range(nr,fr,window_length,search_range,step,min_r2,max_residual)
    n_prof = size(nr,2);
    starts = search_range(1):step:search_range(2)-window_length+1;
    n_run = length(starts);

    range = NaN(2,n_prof);
    fit.centres = round(starts + (window_length-1)/2);
    fit.r2 = NaN(n_run,n_prof);
    fit.residual = NaN(n_run,n_prof);
    fit.slope = NaN(n_run,n_prof);
    fit.valid = false(n_run,n_prof);
    fit.min_r2 = min_r2;
    fit.max_residual = max_residual;
    fit.search_range = search_range;

    for i = 1:n_prof
        for j = 1:n_run
            ids = starts(j):starts(j)+window_length-1;
            x = nr(ids,i);
            y = fr(ids,i);
            if any(~isfinite(x) | ~isfinite(y) | x <= 0 | y <= 0) || max(x) == min(x)
                continue
            end

            [p,s] = polyfit(x,y,1);
            fit.r2(j,i) = s.rsquared;
            fit.residual(j,i) = s.normr / sqrt(length(y)) / mean(y);
            fit.slope(j,i) = p(1);
        end

        fit.valid(:,i) = fit.r2(:,i) >= min_r2 & fit.residual(:,i) <= max_residual & fit.slope(:,i) > 0;
        edges = diff([false;fit.valid(:,i);false]);
        first = find(edges == 1);
        last = find(edges == -1)-1;
        [length_best,id] = max(last-first+1);

        if ~isempty(id) && length_best >= 3
            range(:,i) = [fit.centres(first(id));fit.centres(last(id))];
        end
    end
end
