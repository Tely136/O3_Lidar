function [range,fit] = linear_pc_range(an,pc,window_length,search_range,step,min_r2,max_residual,reference_rates,max_slope_difference)
    n_prof = size(an,2);
    starts = search_range(1):step:search_range(2)-window_length+1;
    n_run = length(starts);

    range = NaN(2,n_prof);
    fit.centres = round(starts + (window_length-1)/2);
    fit.r2 = NaN(n_run,n_prof);
    fit.residual = NaN(n_run,n_prof);
    fit.slope = NaN(n_run,n_prof);
    fit.slope_difference = NaN(n_run,n_prof);
    fit.reference_range = NaN(2,n_prof);
    fit.reference_slope = NaN(1,n_prof);
    fit.reference_intercept = NaN(1,n_prof);
    fit.reference_r2 = NaN(1,n_prof);
    fit.reference_residual = NaN(1,n_prof);
    fit.valid = false(n_run,n_prof);
    fit.min_r2 = min_r2;
    fit.max_residual = max_residual;
    fit.reference_rates = reference_rates;
    fit.max_slope_difference = max_slope_difference;
    fit.search_range = search_range;

    for i = 1:n_prof
        % Fit the reference using corrected, background-subtracted PC rates.
        reference = pc(:,i) >= reference_rates(1) & pc(:,i) <= reference_rates(2) & isfinite(an(:,i));
        reference(1:search_range(1)-1) = false;
        [first,last] = longest_run(reference);
        if isnan(first)
            continue
        end

        fit.reference_range(:,i) = [first;last];
        ids = first:last;
        x_ref = an(ids,i);
        y_ref = pc(ids,i);
        if length(x_ref) >= window_length && max(x_ref) > min(x_ref)
            [p_ref,s_ref] = polyfit(x_ref,y_ref,1);
            fit.reference_slope(i) = p_ref(1);
            fit.reference_intercept(i) = p_ref(2);
            fit.reference_r2(i) = s_ref.rsquared;
            fit.reference_residual(i) = s_ref.normr / sqrt(length(y_ref)) / mean(y_ref);
        end

        for j = 1:n_run
            ids = starts(j):starts(j)+window_length-1;
            x = an(ids,i);
            y = pc(ids,i);

            if any(~isfinite(x) | ~isfinite(y) | y <= 0) || max(x) == min(x)
                continue
            end

            [p,s] = polyfit(x,y,1);
            fit.r2(j,i) = s.rsquared;
            fit.residual(j,i) = s.normr / sqrt(length(y)) / mean(y);
            fit.slope(j,i) = p(1);
        end

        reference_valid = fit.reference_slope(i) > 0 & fit.reference_r2(i) >= min_r2 & fit.reference_residual(i) <= max_residual;
        fit.slope_difference(:,i) = abs(fit.slope(:,i) ./ fit.reference_slope(i)-1);
        fit.valid(:,i) = reference_valid & fit.r2(:,i) >= min_r2 & fit.residual(:,i) <= max_residual & fit.slope_difference(:,i) <= max_slope_difference;
        [first,last] = longest_run(fit.valid(:,i));

        % Require at least three consecutive passing windows.
        if last-first+1 >= 3
            range(:,i) = [fit.centres(first);fit.centres(last)];
        end
    end
end

function [first,last] = longest_run(mask)
    edges = diff([false;mask;false]);
    starts = find(edges == 1);
    ends = find(edges == -1)-1;
    [~,id] = max(ends-starts+1);

    first = NaN;
    last = NaN;
    if ~isempty(id)
        first = starts(id);
        last = ends(id);
    end
end
