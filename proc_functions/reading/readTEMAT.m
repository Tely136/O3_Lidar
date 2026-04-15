function result = readTEMAT(fname)
    temp_f = load(fname);
    
    % TODO: fix to make it not need harcoded ids
    % TODO: also make this general to not DIAL, loop over all active
    % datasets and output data and configuration
    FR_on = cat(4,temp_f.TR(1).data(:,1),temp_f.TR(1).data(:,2));
    FR_off = cat(4,temp_f.TR(2).data(:,1),temp_f.TR(2).data(:,2));
    NR_on = cat(4,temp_f.TR(3).data(:,1),temp_f.TR(3).data(:,2));
    NR_off = cat(4,temp_f.TR(4).data(:,1),temp_f.TR(4).data(:,2));

    % TODO: remove harcoded ids and make it take all active datasets
    ids = 1:4;

    config = struct;
    for i = 1:length(ids)
        config.shots(i) = temp_f.acquired_shots;
        config.bins(i) = temp_f.TR(i).bins;
        config.adcbits(i) = temp_f.TR(i).adcbits;
        config.range(i) = temp_f.TR(i).range;
        config.binwidth(i) = temp_f.TR(i).binWidth;
    end
    
    result.time = [temp_f.start_time temp_f.start_time + diff([temp_f.start_time temp_f.end_time])/2 temp_f.end_time];
    result.data = cat(3,FR_on,FR_off,NR_on,NR_off);
    result.config = config;
end