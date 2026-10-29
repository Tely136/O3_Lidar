function [an_bin_out,pc_P_out,times_avg,ids_avg,an_B,pc_B,pc_uB,pc_rb,pc_R_sum,n_shots_total] = preprocess(an_R,pc_R,times,dz,n_avg,n_shots,fd)
    n_bins = size(an_R,1);
    n_prof = size(an_R,2);
    n_channels = size(pc_R,3);
    n_prof_avg = floor(n_prof / n_avg);

    z = (1:n_bins) .* dz;
    times = times(:,2);

    if isscalar(n_shots)
        n_shots = n_shots .* ones(n_prof,n_channels);
    end

    an_bin_out = NaN(n_bins,n_prof_avg,n_channels);
    pc_R_sum = NaN(n_bins,n_prof_avg,n_channels);
    n_shots_total = NaN(1,n_prof_avg,n_channels);
    times_avg = NaT(n_prof_avg,1);
    ids_avg = NaN(n_prof,1);

    % Accumulate complete groups of files; the remainder is left unused.
    for i = 1:n_prof_avg
        ids = (i-1)*n_avg+1:(i-1)*n_avg+n_avg;
        ids_avg(ids) = i;

        an_bin_out(:,i,:) = sum(an_R(:,ids,:),2);
        pc_R_sum(:,i,:) = sum(pc_R(:,ids,:),2);
        n_shots_total(1,i,:) = sum(n_shots(ids,:),1);
        times_avg(i) = mean(times(ids),1);
    end

    % Correct the accumulated photon counts using their total exposure.
    pc_r = scale_binary_pc(pc_R_sum,dz,n_shots_total);
    pc_P_out = dt_correction(pc_R_sum,pc_r,fd);
    pc_P_out(~isfinite(pc_r) | pc_r >= fd | pc_r < 0) = NaN;

    an_B = bg_correction(z,an_bin_out,7000);
    [pc_B,~,pc_uB,pc_rb] = bg_correction(z,pc_P_out,7000);
end
