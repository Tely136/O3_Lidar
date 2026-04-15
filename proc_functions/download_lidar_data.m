function download_lidar_data(date, download_path)
    current_dir = pwd;
    cd(download_path);
    
    success = false;
    
    site = ['https://datadb.noaacrest.org/public/ORSL/Archived_Lidar/O3DIAL/', ...
        year(date), '/', string(date, 'uuuMMdd')];
        command = strjoin(['"', current_dir, '\wget" -r -q -np -nH --cut-dirs=5 -R "index.html*" ', ...
            site, ' --no-check-certificate'], '');

    [s, ~] = system(command);
        
    if s == 0
        success = true;
    end
  
    
    if ~success
        disp('Error downloading lidar data')
    end
    
    cd(current_dir);
end

