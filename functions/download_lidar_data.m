function download_lidar_data(date, download_path)
    current_dir = pwd;
    cd(download_path);
    
    folderNames = {'txt', 'TXT'};
    success = false;
    
    for i = 1:numel(folderNames)
        site = ['https://datadb.noaacrest.org/public/ORSL/Archived_Lidar/O3DIAL/', ...
            year(date), '/', string(date, 'uuuMMdd'), '/', folderNames{i}];
        command = strjoin(['"', current_dir, '\wget" -r -q -np -nH --cut-dirs=5 -R "index.html*" ', ...
            site, ' --no-check-certificate'], '');
        [s, ~] = system(command);
        if s == 0
            success = true;
            break;
        end
    end
    
    if ~success
        disp('Error downloading lidar data')
    end
    
    cd(current_dir);

end

