function [an_bin_out,pc_P_out,times_avg,an_B,pc_B,pc_uB,pc_rb] = preprocess(an_R,pc_R,times,dz,n_avg,n_shots,fd)
    n_bins = size(an_R,1);
    n_prof = size(an_R,2);
    
    z = (1:n_bins) .* dz;

    % Use midpoint time to represent time
    % TODO: make use of start and end times
    times = times(:,2);
    
    % time-average data, this doesn't handle discontinuous data
    n_prof_avg = floor(n_prof / n_avg);
    mod(n_prof,n_avg); % can check if remainder is >= n_avg/2 and add another profile
    
    pc_r = scale_binary_pc(pc_R,dz,n_shots);
    pc_P = dt_correction(pc_R,pc_r,fd);
   
    
    an_bin_out = NaN(n_bins,n_prof_avg,4);
    pc_P_out = NaN(n_bins,n_prof_avg,4);
    
    times_avg = NaT(n_prof_avg,1);
    for i = 1:n_prof_avg
        ids = (i-1)*n_avg+1:(i-1)*n_avg+n_avg;
    
        an_bin_out(:,i,:) = mean(an_R(:,ids,:),2);
        pc_P_out(:,i,:) = mean(pc_P(:,ids,:),2);
    
        times_avg(i) = mean(times(ids),1);
    end

    an_B = bg_correction(z,an_bin_out,7000);
    [pc_B,~,pc_uB,pc_rb] = bg_correction(z,pc_P_out,7000);
end