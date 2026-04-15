function result = readLicelASCII(fname)
    fid = fopen(fname, 'rb');
    
    header = strings(0);
    
    line = 1;
    while true
        line_chars = fgetl(fid);
        header(line) = line_chars;
        if isempty(line_chars)
            break
        end
        line = line+1;
    end
    
    l3_vals = sscanf(header(3), '%d %d %d %d %d');
    n_datasets = l3_vals(5);
    n_bins = NaN(n_datasets,1);
    an_pc = NaN(n_datasets,1);
    
    for i = 1:n_datasets
        dataset_info = textscan(header(3+i), '%d %d %d %d %d %s %f %s %d %d %d %d %d %d %d %f %s');
        n_bins(i) = dataset_info{4};
        an_pc(i) = dataset_info{2};
    end
    
    max_bins = max(n_bins);
    data = NaN(max_bins,n_datasets);
    
    for i = 1:max_bins
        temp_data = fgetl(fid);
        temp_data = double(split(string(temp_data)));
        if any(size(temp_data(1:end-1)) == 0) 
            continue
        end
        data(i,:) = temp_data(1:end-1);
        fgetl(fid);
    end

    fclose(fid);

    % TODO: add the actual data to this config
    config = struct;
    config.shots = 1;
    config.bins = max_bins;
    config.adcbits = 1;
    config.range = 1;
    config.binwidth = 1;

    temp_str = strsplit(header(2));

    start_date = temp_str(3);
    start_t = temp_str(4);

    end_date = temp_str(5);
    end_t = temp_str(6);

    start_time = datetime(strjoin([start_date start_t]),'InputFormat','dd/MM/yyyy HH:mm:ss');
    end_time = datetime(strjoin([end_date end_t]),'InputFormat','dd/MM/yyyy HH:mm:ss');

    result.time = [start_time start_time + diff([start_time end_time])/2 end_time];
    result.data = data;
    result.config = config;
end
