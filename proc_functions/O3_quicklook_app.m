function O3_quicklook_app(app)
    % disp('test')

    % TODO: change this to O3_quicklook_wrapper or something
    % add indication in app that there are files present and when it is
    % processing
    % modify code so it take single time averaged profile at a time
    % make it a seperate function called by this one

    files = dir(app.DirectorytoProcessEditField.Value);
    num_files = length(find([files.isdir] == 0));
    app.FilesEditField.Value = num_files;

    % TODO: add check that ozone plot is the active screen
    if num_files <= app.n_files || ~(app.TabGroup.SelectedTab == app.AcquisitionTab && app.TabGroup2.SelectedTab == app.OzoneTimeHeightQuicklookTab)
        return; % Exit the function if there are not enough files
    end

    % update UI to show that data is being processed and how many data
    % files are in the working directory
    app.EditField_3.Value = 'Processing';
    drawnow;

    n_channels = 8;
    min_avg = app.minutestoaverageEditField.Value;
    bg_bins = 100;
    td = 270;
    dz = 3.75;
    min_toggle = 3;
    max_toggle = 40;
    
    % On and off wavelengths
    lam_on = 287.2; % nm
    lam_off = 299.1; % nm
    
    fr_fl_M1 = 11; fr_fl_M2 = 321; fr_fl_h1 = 1.2; fr_fl_h2 = 10;
    nr_fl_M1 = 11; nr_fl_M2 = 53; nr_fl_h1 = .2; nr_fl_h2 = 1.2;
    
    o3_params_fr = params_struct(lam_on,lam_off,n_channels,min_avg,bg_bins,td,dz,min_toggle,max_toggle,fr_fl_M1,fr_fl_M2,fr_fl_h1,fr_fl_h2);
    o3_params_nr = params_struct(lam_on,lam_off,n_channels,min_avg,bg_bins,td,dz,min_toggle,max_toggle,nr_fl_M1,nr_fl_M2,nr_fl_h1,nr_fl_h2);


    [raw_data,times,configs] = read_lidar_datafiles(app.DirectorytoProcessEditField.Value,app.FileFormatButtonGroup.SelectedObject.Text);

    
    % Use midpoint time to represent time
    % TODO: make use of start and end times
    times = times(:,2);
    
    % Convert each dataset to physical units
    num_bins = size(raw_data,1);
    num_files = size(raw_data,2);
    num_devices = size(raw_data,3);
    scaled_data = NaN(size(raw_data));
    h = NaN(num_bins,num_files,num_devices);
    
    for i = 1:num_files
        temp_config = configs(i);
        for j = 1:num_devices
            % TODO: remove hardcoded 40 MHz
            scaled_data(:,i,j,1) = scale_binary_analog(raw_data(:,i,j,1),temp_config.range(j),temp_config.adcbits(j),temp_config.shots(j));
            scaled_data(:,i,j,2) = scale_binary_pc(raw_data(:,i,j,2),40,temp_config.shots(j));
    
            h(:,i,j) = (1:num_bins).*temp_config.binwidth(j);
        end
    end
    hkm = h./1000;
    % TODO: potentially handle different bin width for different profiles,
    % right now hkm is a single profile
    hkm = hkm(:,1,1);
    
    out_data_fr = struct;
    out_data_fr.an_on = scaled_data(:,:,1,1);
    out_data_fr.pc_on = scaled_data(:,:,1,2);
    out_data_fr.hkm = hkm;
    
    out_data_fr.an_off = scaled_data(:,:,2,1);
    out_data_fr.pc_off = scaled_data(:,:,2,2);
    
    
    out_data_nr = struct;
    out_data_nr.an_on = scaled_data(:,:,3,1);
    out_data_nr.pc_on = scaled_data(:,:,3,2);
    out_data_nr.hkm = hkm;
    
    out_data_nr.an_off = scaled_data(:,:,4,1);
    out_data_nr.pc_off = scaled_data(:,:,4,2);
    
    output_fr = O3_quicklook(out_data_fr,times,configs,o3_params_fr);
    output_nr = O3_quicklook(out_data_nr,times,configs,o3_params_nr);
    times_avg = output_fr.on.times_avg;
    
    % merge nr and fr
    start_merge = app.startmergeheightkmEditField.Value;
    end_merge = app.endmergeheightkmEditField.Value;
    % [No3_merged,w_fr_nr] = merge_profiles(output_nr.No3,output_fr.No3,hkm,start_merge,end_merge);
    qo3_merged = merge_profiles(output_nr.qo3,output_fr.qo3,hkm,start_merge,end_merge);

    %% 
    % % O3 Time-height Plot
    ND_plot_ax(app.UIAxes2,times_avg,hkm,qo3_merged,[0 120],app.YlimitSlider_2.Value,'');
    % 
    % % Signal Plot
    % % app.IDSpinner.Limits = [1, length(times_avg)];
    % % id = app.IDSpinner.Value;
    % % plot(app.UIAxes4,data_glued(:,id,1),hkm);
    % % % hold(app.UIAxes4,"on");
    % % % 
    % % % hold(app.UIAxes4,"off");
    % % xlim(app.UIAxes4,app.XlimitSlider_2.Value);
    % % ylim(app.UIAxes4,app.YlimitSlider_3.Value);
    % 
    % % O3 Profile Plot
    % app.IDSpinner_2.Limits = [1, length(times_avg)];
    % id = app.IDSpinner_2.Value;
    % loglog(app.UIAxes5,qo3_merged(:,id),hkm);
    % % hold(app.UIAxes5,"on");
    % % 
    % % hold(app.UIAxes5,"off");
    % 
    % xlim(app.UIAxes5,app.XlimitSlider_3.Value);
    % ylim(app.UIAxes5,app.YlimitSlider_4.Value);
    % 
    app.EditField_3.Value = 'Idle';
    app.n_files = num_files;
end