function [data,times,configs] = read_lidar_datafiles(folder_path,file_format,prefix)
    arguments
        folder_path string
        file_format string
        prefix string = "ol"
    end

    % create FDS of ozone lidar datafiles
    switch file_format
        case "mat"
            fds = fileDatastore(folder_path, 'ReadFcn', @readTEMAT, 'FileExtensions','.mat');

        case "bin"
            dir_info = dir(folder_path);
            filenames = string({dir_info.name});
            ol_ids = ~[dir_info().isdir] & startsWith(filenames,prefix);
            
            ol_files = dir_info(ol_ids);
            ol_paths = fullfile(string({ol_files.folder}), string({ol_files.name}));
            fds = fileDatastore(ol_paths,'ReadFcn',@readLicelBinary);

        case "txt"
            fds = fileDatastore(folder_path, 'ReadFcn', @readLicelASCII, 'FileExtensions','.txt');

    end
    
    % preview data within fds
    preview_data = preview(fds);
    
    % get filenames and number of files
    fullFileNames = fds.Files;
    num_files = length(fullFileNames);
    n_channels = size(preview_data.data,3);
    n_datasets = size(preview_data.data,4);

    % initialize raw data array
    data = NaN(2^14,num_files,n_channels,n_datasets);
    
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
            data(1:Data.config.bins(j),i,j,:) = Data.data(:,:,j,:);
        end
    end
end