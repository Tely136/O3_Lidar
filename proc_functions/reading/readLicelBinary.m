function result = readLicelBinary(fname)
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
    n = 16;
    
    l3_vals = sscanf(header(3), '%d %d %d %d %d');
    n_datasets = l3_vals(5);
    tr = NaN(n,1);
    n_bins = NaN(n,1);
    an_pc = NaN(n,1);
    anrange = NaN(n,1);
    disc = NaN(n,1);
    adcbits = NaN(n,1);
    binwidth = NaN(n,1);
    shots = NaN(n,1);

    for i = 1:n_datasets
        dataset_info = textscan(header(3+i), '%d %d %d %d %d %d %f %s %d %d %d %d %d %d %f %s');
        dscr = char(dataset_info{16});
        tr(i) = str2double(dscr(3));
        if tr(i)>3% need to fix this later
            tr(i)=3;
        end
        tr_idx = tr(i) + 1;

        an_pc(i) = dataset_info{2};
        n_bins(tr_idx) = dataset_info{4};
        binwidth(tr_idx) = dataset_info{7};
        shots(tr_idx) = dataset_info{14};

        temp = dataset_info{15};
        if dataset_info{2} == 0
        adcbits(tr_idx) = dataset_info{13};

            switch temp
                case 0.5
                    anrange(tr_idx) = 0;
                case 0.1
                    anrange(tr_idx) = 1;
                case 0.02
                    anrange(tr_idx) = 2;
            end
        elseif dataset_info{2} == 1
            disc(tr_idx) = temp;
        end

    end
    n_tr = length(unique(tr(~isnan(tr))));

    n_bins(n_tr+1:end) = [];
    % an_pc(n_tr+1:end) = [];
    anrange(n_tr+1:end) = [];
    adcbits(n_tr+1:end) = [];
    binwidth(n_tr+1:end) = [];
    shots(n_tr+1:end) = [];
    
    max_bins = max(n_bins);
    data = NaN(max_bins,1,n_tr,2);
    for i = 1:n_datasets
        data(1:n_bins(tr(i)+1),1,tr(i)+1,an_pc(i)+1) = fread(fid, n_bins(tr(i)+1), 'uint32');
        fgetl(fid);
    end
    fclose(fid);

    config = struct;
    config.shots = shots;
    config.bins = n_bins;
    config.adcbits = adcbits;
    config.range = anrange;
    config.binwidth = binwidth;

    temp_str = strsplit(header(2));

    start_date = temp_str(3);
    start_t = temp_str(4);

    end_date = temp_str(5);
    end_t = temp_str(6);

    try % this try catch needs to be removed
        start_time = datetime(strjoin([start_date start_t]),'InputFormat','dd/MM/yyyy HH:mm:ss');
        end_time = datetime(strjoin([end_date end_t]),'InputFormat','dd/MM/yyyy HH:mm:ss');
    catch
        start_date = temp_str(4);
        start_t = temp_str(5);
    
        end_date = temp_str(6);
        end_t = temp_str(7);

        start_time = datetime(strjoin([start_date start_t]),'InputFormat','dd/MM/yyyy HH:mm:ss');
        end_time = datetime(strjoin([end_date end_t]),'InputFormat','dd/MM/yyyy HH:mm:ss');

    end

    result.time = [start_time start_time + diff([start_time end_time])/2 end_time];
    result.data = data;
    result.config = config;
end
