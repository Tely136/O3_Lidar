function [an_data,pc_data,times,configs] = read_lidar_datafiles(folder_path,file_format,prefix,start_time,end_time)
    arguments
        folder_path string
        file_format string
        prefix string = "ol"
        start_time datetime = datetime(2000,1,1);
        end_time datetime = datetime(2100,1,1);
    end

    % create FDS of ozone lidar data files
    switch file_format
        case "mat"
            fds = fileDatastore(folder_path, 'ReadFcn', @readTEMAT, 'FileExtensions','.mat');

        case "bin"
            dir_info = dir(folder_path);
            files = dir_info(~[dir_info.isdir]);
            filenames = string({files.name});

            filename_dates = parse_filename(filenames);
            
            ol_ids = startsWith(filenames,prefix) & filename_dates >= start_time & filename_dates <= end_time;
            
            ol_files = files(ol_ids);
            ol_paths = fullfile(string({ol_files.folder}), string({ol_files.name}));
            fds = fileDatastore(ol_paths,'ReadFcn',@readLicelBinary);

        case "txt"
            fds = fileDatastore(folder_path, 'ReadFcn', @readLicelASCII, 'FileExtensions','.txt');

        case "list"
            filename_list = folder_path;
            fds = fileDatastore(filename_list,'ReadFcn',@readLicelBinary);
    end
    
    % preview data within fds
    preview_data = preview(fds);
    
    % get filenames and number of files
    fullFileNames = fds.Files;
    num_files = length(fullFileNames);
    n_channels = size(preview_data.data,3);
    % n_datasets = size(preview_data.data,4);

    % initialize raw data array
    exp = 10;
    an_data = NaN(2^exp,num_files,n_channels);
    pc_data = NaN(2^exp,num_files,n_channels);

    % initialize array for times
    times = NaT(num_files,3);
    
    % initialize struct array for config
    configs(num_files) = preview_data.config;
    
    % loop over files, add data and time to arrays
    for i = 1:num_files
        Data = read(fds);
    
        times(i,:)       = Data.time;
        configs(i)       = Data.config;
        for j = 1:n_channels
            temp = Data.data(:,:,j,:);

            an_data(1:Data.config.bins(j),i,j) = temp(:,:,:,1);
            pc_data(1:Data.config.bins(j),i,j) = temp(:,:,:,2);
        end
    end
end