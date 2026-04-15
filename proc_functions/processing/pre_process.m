function data = pre_process(raw_data,times,avg_time,td,bg_bins,PC)
    n_bin = size(raw_data,1);

    % Time average data    
    d = raw_data';
    TT = timetable(times,d);
    TT2 = retime(TT,'regular','mean','TimeStep',avg_time);
    data = TT2.d';
    
    times_avg = TT2.times;
    n_profiles = length(times_avg);
    
    % Dead-time correct PC data
    if PC == true % only for PC data
        data = data./(1-data/td);
    end
    
    % Background Subtraction
    bg_ind = n_bin-bg_bins:n_bin;
    solar_background = NaN(n_profiles);
    
    for i = 1:n_profiles
        bg = mean(data(bg_ind,i));
        solar_background(i) = bg;
        
        data(:,i) = data(:,i) - bg; 
        data(data(:,i)<0,i) = NaN;
    end
end