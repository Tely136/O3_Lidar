function result = tempo_o3_colocate(files,params)
    n_files = length(files);

    result = struct;
    result.o3_prof = NaN(24,n_files);
    result.alt = NaN(24,n_files);
    result.alt_bounds = NaN(24,2,n_files);
    result.time = NaT(1,n_files,'TimeZone','UTC');

    for i = 1:n_files
        temp_data = read_tempo_o3(files(i));

        temp_o3_prof = temp_data.o3_prof;

        temp_alt = temp_data.o3_prof_alt;
        temp_alt_bounds = temp_data.o3_prof_alt_bound;

        temp_lat = temp_data.lat;
        temp_lon = temp_data.lon;

        [n_row,n_col] = size(temp_lat);
        % n_level = size(temp_alt,1);
    
        distance_m = distance(params.lat,params.lon,temp_lat,temp_lon,wgs84Ellipsoid);

        idx_valid = find(distance_m <= params.max_distance_m & temp_data.cld_frac <= params.max_cld_frac & temp_data.sza <= params.max_sza & temp_data.vza <= params.max_vza);
        
        if isempty(idx_valid)
            % warning("no valid profiles for %s", files(i))
            continue
        end

        [~,temp_idx_closest] = min(distance_m(idx_valid),[],"all");
        idx_closest = idx_valid(temp_idx_closest);

        % [rows_valid,cols_valid] = ind2sub([n_row,n_col],idx_valid);
        [row_closest,col_closest] = ind2sub([n_row,n_col],idx_closest);

        result.o3_prof(:,i) = temp_o3_prof(:,row_closest,col_closest);
        result.alt(:,i) = temp_alt(:,row_closest,col_closest);
        result.alt_bounds(:,:,i) = temp_alt_bounds(:,:,row_closest,col_closest)';
        result.time(i) = datetime(temp_data.time(col_closest),'ConvertFrom','epochtime','Epoch','1980-01-06','TimeZone','UTC');
    end

    nan_cols = all(isnan(result.o3_prof),1);

    result.o3_prof(:,nan_cols) = [];
    result.alt(:,nan_cols) = [];
    result.alt_bounds(:,:,nan_cols) = [];
    result.time(nan_cols) = [];
end