function write_o3_nasa(fname,folder,config)
    fileID = fopen(fullfile(folder,fname), 'w');
    
    header_lines(fileID,fname,config)
    o3_lines(fileID,config)
    fwrite(fileID,newline);
    
    o3_data(fileID,config);
    fwrite(fileID,newline);
    
    fclose(fileID);
end


function header_lines(fid,fname,config)
    start_time = string(datetime(config.start_time,"Format","dd/MM/uuuu HH:mm:ss"));
    stop_time = string(datetime(config.end_time,"Format","dd/MM/uuuu HH:mm:ss"));
    location = sprintf("%-8s",'CCNY');


    height = 10;
    longitude = -77.1;
    latitude = 44.6;
    zenith = 0;
    azimuth = 0;

    l1_shots = config.max_shots;
    l1_rep = config.repRate;
    l2_shots = 0;
    l2_rep = 0;
    l3_shots = 0;
    l3_rep = 0;
    n_ds = 8;

    line1 = ' ' + fname;
    line2 = strjoin(["" location start_time stop_time sprintf('%04d %011.6f %011.6f %04.1f %04.1f',height,longitude,latitude,zenith,azimuth)]);
    line3 = sprintf(' %07d %04d %07d %04d %02d %07d %04d %07d %04d',l1_shots,l1_rep,l2_shots,l2_rep,n_ds,l3_shots,l3_rep,0,0);

    fprintf(fid, '%s\n%s\n%s\n%s\n', line1, line2, line3);
end


function o3_data(fid,config)
    fields = ["FR_on" "FR_off" "NR_on" "NR_off"];
    
    for i = 1:4
        ds_name = fields(i);

        fwrite(fid,config.(ds_name).data(:,1),'uint32');
        fwrite(fid,newline);

        fwrite(fid,config.(ds_name).data(:,2),'uint32');
        fwrite(fid,newline);
    end

end


function o3_lines(fid,config)
    fields = ["FR_on" "FR_off" "NR_on" "NR_off"];
    bins = config.bins;
    shots = config.acquired_shots; % change acquire data to put the actual acquired shots in config, not result, use that here

    for i = 1:4 
        ds_name = fields(i);
        ds = config.(ds_name);

        switch ds.range
            case 0
                an_range = 0.5;
            case 1
                an_range = 0.1;
            case 2
                an_range = 0.02;
        end
    
        if ds_name == "FR_on" || ds_name == "FR_off"
            wavelength = config.lam_on;
        elseif ds_name == "NR_on" || ds_name == "NR_off"
            wavelength = config.lam_off;
        end

        an_line = ds_line(ds.TR,ds.enable,0,1,bins,0,ds.pmtv,ds.binWidth,wavelength,'o',ds.binShift,0,ds.adcBits,shots,an_range,'BT');
        pc_line = ds_line(ds.TR,ds.enable,1,1,bins,0,ds.pmtv,ds.binWidth,wavelength,'o',ds.binShift,0,ds.adcBits,shots,ds.discr,'BC');
        
        fwrite(fid,an_line);
        fwrite(fid,pc_line);
    end
end


function line = ds_line(tr,active,type,laser,bins,laser_pol,hv,bin_width,lambda,pol,bin_shift,bin_shift_dec,adc,shots,range,discr)
    % active is 1 or 0
    active_str = sprintf("%d",active);

    % type is 1,2,3,4 or 5
    type_str = sprintf("%d",type);

    % laser is 1,2,3, or 4
    laser_str = sprintf("%d",laser);

    % bins is a 5 digit integer
    bins_str = sprintf("%05d",bins);

    % laser_pol is 1,2,3, or 4
    laser_pol_str = sprintf("%d",laser_pol);

    % hv is a 4 digit integer
    hv_str = sprintf("%04d",hv);

    % bin_width is 4 digit, 2 after decimal place
    bin_width_str = sprintf("%04.2f",bin_width);

    % lambda is 5 digits and a period, pol is 1 character added right after
    lambda_pol_str = sprintf("%05d.%c",lambda,pol);

    % bin_shift is 2 digits
    bin_shift_str = sprintf("%02d",bin_shift);

    % bin_sift_dec is 3 decimal points of the bin_shift
    bin_shift_dec_str = sprintf("%03d",bin_shift_dec);

    % adc is 2 digits of analog, 0 otherwise
    if type == 0
        adc_str = sprintf("%02d",adc);
    else
        adc_str = sprintf("%02d",0);
    end

    % shots is 6 digits
    shot_str = sprintf("%06d",shots);

    % range is 1 period 3 digits if analog, 1 period 4 digits if photon counting
    if ismember(type,[0 2 4 5])
        range_str = sprintf("%04.3f",range);
    elseif ismember(type,[1 3])
        range_str = sprintf("%05.3f",range);
    end

    % discr is 2-3 characters to describe the dataset followed by thr TR number
    disc_str = sprintf("%s%c",discr,dec2hex(tr));

    line = strjoin(["" active_str type_str laser_str bins_str laser_pol_str...
                    hv_str bin_width_str lambda_pol_str "0 0" bin_shift_str...
                    bin_shift_dec_str adc_str shot_str range_str disc_str newline]);
end
