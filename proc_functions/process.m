function out = process(input_path,output_path)
    repo_dir = fileparts(fileparts(mfilename('fullpath')));
    model_dir = fullfile(fileparts(repo_dir),'signal_analysis');
    addpath(genpath(fullfile(repo_dir,'proc_functions')));
    addpath(genpath(fullfile(model_dir,'functions')));

    %% Signal preparation
    % TODO: move parameters into config file
    data_format = "bin";
    prefix = "ol";
    n_avg = 30;
    fd = 270;                           % inverse dead time (MHz)
    tau_d = 1 / (fd * 1e6);             % dead-time in seconds

    start_merge = 1.5;
    end_merge = 2;

    min_toggle = .1;
    max_toggle = 60;

    % Relative standard uncertainties
    u_sigma_o3_rel = 0.02;
    u_sigma_m_rel = 0.02;
    uTau_rel = 0.10;                    % provisional standard uncertainty until calibrated

    uTau = uTau_rel * tau_d;

    % Channel IDs
    on_fr_id = 1;
    off_fr_id = 2;
    on_nr_id = 3;
    off_nr_id = 4;

    % Bin width
    dz = 3.75;

    % Window sizes for an-pc merge
    window_length_fr = 500;
    window_length_nr = 100;

    % Read raw binary analog and photon counting data
    [an_R,pc_R,times,configs] = read_lidar_datafiles(input_path,data_format,prefix);
    n_bins = size(an_R,1);
    n_prof = size(an_R,2);

    % Shots per file and channel from the Licel headers.
    n_shots = reshape([configs.shots],size(pc_R,3),[])';

    n_prof_avg = floor(n_prof / n_avg);

    % Alt vector
    z = (1:n_bins) .* dz;
    z_km = z.*1e-3;

    % change below comment
    % Sum raw files, correct dead time, then fit the accumulated background.
    [an_bin_sum, pc_P_sum, times_avg, ids_avg, an_B_sum, pc_B_sum, pc_uB, pc_rb, pc_R_sum, n_shots_total] = preprocess(an_R,pc_R,times,dz,n_avg,n_shots,fd);
    an_bin_bg = an_bin_sum - an_B_sum;
    pc_P_bg = pc_P_sum - pc_B_sum;

    an_bin_bg(an_bin_bg <= 0) = NaN;
    pc_P_bg(pc_P_bg <= 0) = NaN;

    % Convert binary analog signal and background to mV
    an_mV_avg = scale_binary_analog(an_bin_bg,0,16,n_shots_total); % TODO loop over devices and scale each one according to config values
    an_mV_bg_avg = scale_binary_analog(an_B_sum,0,16,n_shots_total);

    % Convert binary photon counting signal and background MHz
    pc_p_avg = scale_binary_pc(pc_P_bg,dz,n_shots_total);
    pc_p_bg_avg = scale_binary_pc(pc_B_sum,dz,n_shots_total);

    % Grid
    z_grid = repmat((1:n_bins)',1,n_prof_avg) .* dz;

    % Extract analog signal and background per TR
    an_on_fr = an_mV_avg(:,:,on_fr_id);
    an_off_fr = an_mV_avg(:,:,off_fr_id);
    an_on_nr = an_mV_avg(:,:,on_nr_id);
    an_off_nr = an_mV_avg(:,:,off_nr_id);

    an_on_fr_bg = an_mV_bg_avg(:,:,on_fr_id);
    an_off_fr_bg = an_mV_bg_avg(:,:,off_fr_id);
    an_on_nr_bg = an_mV_bg_avg(:,:,on_nr_id);
    an_off_nr_bg = an_mV_bg_avg(:,:,off_nr_id);

    % Extract photon counting signal and background per TR
    pc_on_fr = pc_p_avg(:,:,on_fr_id);
    pc_off_fr = pc_p_avg(:,:,off_fr_id);
    pc_on_nr = pc_p_avg(:,:,on_nr_id);
    pc_off_nr = pc_p_avg(:,:,off_nr_id);

    pc_on_fr_bg = pc_p_bg_avg(:,:,on_fr_id);
    pc_off_fr_bg = pc_p_bg_avg(:,:,off_fr_id);
    pc_on_nr_bg = pc_p_bg_avg(:,:,on_nr_id);
    pc_off_nr_bg = pc_p_bg_avg(:,:,off_nr_id);
    
    % Merge analog and photon counting signals for far range channels
    [mg_on_fr,coeff_on_fr,r2_on_fr]    = an_pc_merge(an_on_fr, pc_on_fr, window_length_fr, min_toggle, max_toggle);
    [mg_off_fr,coeff_off_fr,r2_off_fr] = an_pc_merge(an_off_fr, pc_off_fr, window_length_fr, min_toggle, max_toggle);

    % Merge analog and photon counting signals for near range channels
    [mg_on_nr,coeff_on_nr,r2_on_nr]    = an_pc_merge(an_on_nr, pc_on_nr, window_length_nr, min_toggle, max_toggle);
    [mg_off_nr,coeff_off_nr,r2_off_nr] = an_pc_merge(an_off_nr, pc_off_nr, window_length_nr, min_toggle, max_toggle);

    % TODO: Automatically find best FR-NR merge for PC processed ozone
    % [range_on_fr,pc_fit_on_fr] = linear_pc_range(an_on_fr,pc_on_fr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);
    % [range_off_fr,pc_fit_off_fr] = linear_pc_range(an_off_fr,pc_off_fr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);
    % [range_on_nr,pc_fit_on_nr] = linear_pc_range(an_on_nr,pc_on_nr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);
    % [range_off_nr,pc_fit_off_nr] = linear_pc_range(an_off_nr,pc_off_nr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);

    [~,start_merge_id] = min(abs(z_km-start_merge));
    [~,end_merge_id] = min(abs(z_km-end_merge));

    %% O3 Processing
    % On and off wavelengths
    lam_on = 287.2 * 1e-9; % m
    lam_off = 299.1 * 1e-9; % m

    % FR and NR filter widths
    fr_fl_M1 = 11; fr_fl_M2 = 351;
    % fr_fl_M1 = 75; fr_fl_M2 = 75; % constant 
    fr_fl_h1 = 1.2; fr_fl_h2 = 10;
    
    nr_fl_M1 = 11; nr_fl_M2 = 121; 
    % nr_fl_M1 = 45; nr_fl_M2 = 45; % constant
    nr_fl_h1 = .2;  nr_fl_h2 = 1.2;

    % FR and NR filters
    C_fr = filter_gen(fr_fl_M1,fr_fl_M2,fr_fl_h1,fr_fl_h2,z_km);
    C_nr = filter_gen(nr_fl_M1,nr_fl_M2,nr_fl_h1,nr_fl_h2,z_km);

    % Air number density and temperature
    Na = airNumberDensity(z_grid);
    Ta = atmosisa(z_grid);

    % Differential ozone absorption cross section as function of temperature
    BDM = load(fullfile(model_dir,'cross_sections','ozone','bdm','bdm_287_299.mat'));
    BDM_t = BDM.BDM_t;
    sigma_o3_on = BDM.BDM_c287 .* 1e-4;
    sigma_o3_off = BDM.BDM_c299 .* 1e-4;

    sigma_o3_on_t = NaN(n_bins,n_prof_avg);
    sigma_o3_off_t = NaN(n_bins,n_prof_avg);
    for i=1:n_prof_avg
        sigma_o3_on_t(:,i) = fillmissing(interp1(BDM_t, sigma_o3_on, Ta(:,i)), "nearest");
        sigma_o3_off_t(:,i) = fillmissing(interp1(BDM_t, sigma_o3_off, Ta(:,i)), "nearest");
    end

    dsigma_o3 = sigma_o3_on_t - sigma_o3_off_t;

    % Rayleigh extinction correction
    sigma_m_on = rayleigh_cross_section(lam_on);
    sigma_m_off = rayleigh_cross_section(lam_off);
    
    dsigma_m = sigma_m_on - sigma_m_off;

    N_o3_dm = rayleigh_correction(Na,dsigma_m,dsigma_o3);

    % Ozone number density calculated from analog signal
    N_o3_fr_an = dial(an_on_fr,an_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_an = dial(an_on_nr,an_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    % Ozone number density calculated from photon counting signal
    N_o3_fr_pc = dial(pc_on_fr,pc_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_pc = dial(pc_on_nr,pc_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    % Ozone number density calculated from analog photon counting merged signal
    N_o3_fr_mg = dial(mg_on_fr,mg_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_mg = dial(mg_on_nr,mg_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    % Merge near range and far range for analog calculated profile, photon counting calculated profile, and analog-pc calculated profile
    [N_o3_an_merged,w_o3_an] = merge_profiles(N_o3_nr_an,N_o3_fr_an,start_merge_id,end_merge_id);
    [N_o3_pc_merged,w_o3_pc] = merge_profiles(N_o3_nr_pc,N_o3_fr_pc,start_merge_id,end_merge_id);
    [N_o3_mg_merged,w_o3_mg] = merge_profiles(N_o3_nr_mg,N_o3_fr_mg,start_merge_id,end_merge_id);

    % Convert ozone number density to mixing ratio (ppb)
    q_o3_an_fr = (N_o3_fr_an ./ Na) * 1e9;
    q_o3_an_nr = (N_o3_nr_an ./ Na) * 1e9;

    q_o3_pc_fr = (N_o3_fr_pc ./ Na) * 1e9;
    q_o3_pc_nr = (N_o3_nr_pc ./ Na) * 1e9;

    q_o3_mg_fr = (N_o3_fr_mg ./ Na) * 1e9;
    q_o3_mg_nr = (N_o3_nr_mg ./ Na) * 1e9;

    q_o3_an_merged = (N_o3_an_merged ./ Na) * 1e9;
    q_o3_pc_merged = (N_o3_pc_merged ./ Na) * 1e9;
    q_o3_mg_merged = (N_o3_mg_merged ./ Na) * 1e9;

    % Photon-count uncertainty budget
    uncertainty_options.dz = dz;                                % bin spacing
    uncertainty_options.L = n_shots_total;                      % number of total laser shots per profile
    uncertainty_options.uTau = uTau;                            % dead-time uncertainty
    uncertainty_options.dsigma = dsigma_o3;                     % differential ozone cross-section
    uncertainty_options.channel_correlation = 'independent';    % correlation between on-off channels ('independent' or 'full')
    uncertainty_options.sigma_o3_on = sigma_o3_on_t;            % ozone cross absorption section at lam_on
    uncertainty_options.sigma_o3_off = sigma_o3_off_t;          % ozone absorption cross section at lam_off
    uncertainty_options.sigma_m_on = sigma_m_on;                % molecular extinction cross section at lam_on
    uncertainty_options.sigma_m_off = sigma_m_off;              % molecular extinction cross section at lam_off
    uncertainty_options.u_sigma_o3_rel = u_sigma_o3_rel;        % percent uncertainty in ozone absorption cross section
    uncertainty_options.u_sigma_m_rel = u_sigma_m_rel;          % percent uncertainty in molecular extinction cross section

    % Calculate uncertainty in ozone number density for FR
    uncertainty_options.C = C_fr;
    [uNO3_fr,pc_uncertainty_fr] = uNO3(pc_R_sum,pc_P_bg,N_o3_fr_pc,Na,pc_uB,pc_rb,on_fr_id,off_fr_id,uncertainty_options);

    % Calculate uncertainty in ozone number density for NR 
    uncertainty_options.C = C_nr;
    [uNO3_nr,pc_uncertainty_nr] = uNO3(pc_R_sum,pc_P_bg,N_o3_nr_pc,Na,pc_uB,pc_rb,on_nr_id,off_nr_id,uncertainty_options);

    % Merge FR and NR uncertainty through weighted average
    [uNO3_pc_merged,pc_uncertainty_merged] = merge_pc_uncertainty(pc_uncertainty_nr,pc_uncertainty_fr,w_o3_pc);

    % Convert uncertainty in ozone number density to ozone mixing ratio
    uqO3_fr        = uNO3_fr ./ Na .* 1e9;
    uqO3_nr        = uNO3_nr ./ Na .* 1e9;
    uqO3_pc_merged = uNO3_pc_merged ./ Na .* 1e9;

    %% Output
    out = struct();

    % Raw analog signals
    out.an_R_on_fr  = an_R(:,:,on_fr_id);
    out.an_R_off_fr = an_R(:,:,off_fr_id);
    out.an_R_on_nr  = an_R(:,:,on_nr_id);
    out.an_R_off_nr = an_R(:,:,off_nr_id);

    % Raw PC signals
    out.pc_R_on_fr  = pc_R(:,:,on_fr_id);
    out.pc_R_off_fr = pc_R(:,:,off_fr_id);
    out.pc_R_on_nr  = pc_R(:,:,on_nr_id);
    out.pc_R_off_nr = pc_R(:,:,off_nr_id);

    % Analog signal (mV)
    out.an_on_fr  = an_on_fr;
    out.an_off_fr = an_off_fr;
    out.an_on_nr  = an_on_nr;
    out.an_off_nr = an_off_nr;

    % Analog signal background (mV)
    out.an_on_fr_bg  = an_on_fr_bg;
    out.an_off_fr_bg = an_off_fr_bg;
    out.an_on_nr_bg  = an_on_nr_bg;
    out.an_off_nr_bg = an_off_nr_bg;

    % PC signal (MHz)
    out.pc_on_fr  = pc_on_fr;
    out.pc_off_fr = pc_off_fr;
    out.pc_on_nr  = pc_on_nr;
    out.pc_off_nr = pc_off_nr;

    % PC signal background (MHz)
    out.pc_on_fr_bg  = pc_on_fr_bg;
    out.pc_off_fr_bg = pc_off_fr_bg;
    out.pc_on_nr_bg  = pc_on_nr_bg;
    out.pc_off_nr_bg = pc_off_nr_bg;

    % AN-PC Merge signal (MHz)
    out.mg_on_fr  = mg_on_fr;
    out.mg_off_fr = mg_off_fr;
    out.mg_on_nr  = mg_on_nr;
    out.mg_off_nr = mg_off_nr;

    % Ozone mixing ratio (ppb)
    out.q_o3_an_fr = q_o3_an_fr;
    out.q_o3_an_nr = q_o3_an_nr;
    out.q_o3_pc_fr = q_o3_pc_fr;
    out.q_o3_pc_nr = q_o3_pc_nr;
    out.q_o3_mg_fr = q_o3_mg_fr;
    out.q_o3_mg_nr = q_o3_mg_nr;

    % Ozone number density (m^-3)
    out.N_o3_fr_an = N_o3_fr_an;
    out.N_o3_nr_an = N_o3_nr_an;
    out.N_o3_fr_pc = N_o3_fr_pc;
    out.N_o3_nr_pc = N_o3_nr_pc;
    out.N_o3_fr_mg = N_o3_fr_mg;
    out.N_o3_nr_mg = N_o3_nr_mg;

    % Ozone mixing ratio NR-FR Merged
    out.q_o3_an_merged = q_o3_an_merged;
    out.q_o3_pc_merged = q_o3_pc_merged;
    out.q_o3_mg_merged = q_o3_mg_merged;

    % Ozone number density NR-FR Merged
    out.N_o3_an_merged = N_o3_an_merged;
    out.N_o3_pc_merged = N_o3_pc_merged;
    out.N_o3_mg_merged = N_o3_mg_merged;

    % FR-NR merge weighting coefficients
    out.w_o3_an = w_o3_an;
    out.w_o3_pc = w_o3_pc;
    out.w_o3_mg = w_o3_mg;

    % Rayleigh extinction correction
    out.N_o3_dm = N_o3_dm;

    % AN-PC Merge coefficients
    out.coeff_on_fr  = coeff_on_fr;
    out.coeff_off_fr = coeff_off_fr;
    out.coeff_on_nr  = coeff_on_nr;
    out.coeff_off_nr = coeff_off_nr;

    % AN-PC Merge r^2
    out.r2_on_fr  = r2_on_fr;
    out.r2_off_fr = r2_off_fr;
    out.r2_on_nr  = r2_on_nr;
    out.r2_off_nr = r2_off_nr;

    % Altitude, raw time, averaged time
    out.z_vec = z;
    out.times = times;
    out.times_avg = times_avg;

    out.ids_avg = ids_avg;
    out.n_shots_total = n_shots_total;
    out.pc_R_sum = pc_R_sum;
    out.pc_P_sum = pc_P_bg;

    % Licel configs
    out.configs = configs;

    % SG Filter coefficients
    out.C_fr = C_fr;
    out.C_nr = C_nr;

    % Air number density and temperature
    out.Na = Na;
    out.Ta = Ta;

    % Ozone extinction cross section
    out.sigma_o3_on = sigma_o3_on_t;
    out.sigma_o3_off = sigma_o3_off_t;
    out.dsigma_o3 = dsigma_o3;
    out.sigma_m_on = sigma_m_on;
    out.sigma_m_off = sigma_m_off;
    out.dsigma_m = dsigma_m;

    % Uncertainty from photon-counting signals
    out.uNO3_fr = uNO3_fr;                                  % FR uncertainty in ozone number density
    out.uNO3_nr = uNO3_nr;                                  % NR uncertainty in ozone number density
    out.uNO3_pc_merged = uNO3_pc_merged;                    % FR-NR merged uncertainty in ozone number density

    out.uqO3_fr = uqO3_fr;                                  % FR uncertainty in ozone mixing ratio
    out.uqO3_nr = uqO3_nr;                                  % NR uncertainty in ozone mixing ratio
    out.uqO3_pc_merged = uqO3_pc_merged;                    % FR-NR merged uncertainty in ozone mixing ratio

    out.pc_uncertainty_fr = pc_uncertainty_fr;              % ?
    out.pc_uncertainty_nr = pc_uncertainty_nr;              % ?
    out.pc_uncertainty_merged = pc_uncertainty_merged;      % ?

    save(output_path, '-struct', 'out')
end

