function out = process(input_path,output_path,aerosol)
    if nargin < 3
        aerosol = [];
    end

    % Resolve dependencies from this file, independent of MATLAB's cwd.
    repo_dir = fileparts(fileparts(mfilename('fullpath')));
    model_dir = fullfile(fileparts(repo_dir),'signal_analysis');
    addpath(genpath(fullfile(repo_dir,'proc_functions')));
    addpath(genpath(fullfile(model_dir,'functions')));

    %% Signal preparation
    % parameters to be moved into config file
    data_format = "bin";
    prefix = "ol";
    n_avg = 30;
    fd = 270; % inverse dead time (MHz)
    tau_d = 1 / (fd * 1e6); % seconds
    uTau_rel = 0.10; % provisional standard uncertainty until calibrated
    uTau = uTau_rel * tau_d;
    start_merge = 1.5;
    end_merge = 2;
    min_toggle = .1;
    max_toggle = 60;

    % Channel IDs
    on_fr_id = 1;
    off_fr_id = 2;
    on_nr_id = 3;
    off_nr_id = 4;

    % Relative standard uncertainties; working estimates from Table 3.
    u_sigma_o3_rel = 0.02;
    u_sigma_m_rel = 0.02;

    % Bin width
    dz = 3.75;

    % Window sizes for an-pc merge
    window_length_fr = 500;
    window_length_nr = 100;

    % Analog-PC agreement search for the FR-NR PC merge.
    pc_search_km = [1 5];
    pc_window_length = 201;
    pc_window_step = 20;
    pc_min_r2 = 0.95;
    pc_max_residual = 0.10;
    pc_max_slope_difference = 0.15;

    % Direct FR-NR PC linearity check.
    fr_nr_min_r2 = 0.95;
    fr_nr_max_residual = 0.10;

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

    % Find continuous analog-PC agreement intervals.
    [~,pc_search_start] = min(abs(z_km-pc_search_km(1)));
    [~,pc_search_end] = min(abs(z_km-pc_search_km(2)));
    pc_search_range = [pc_search_start pc_search_end];

    [range_on_fr,pc_fit_on_fr] = linear_pc_range(an_on_fr,pc_on_fr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);
    [range_off_fr,pc_fit_off_fr] = linear_pc_range(an_off_fr,pc_off_fr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);
    [range_on_nr,pc_fit_on_nr] = linear_pc_range(an_on_nr,pc_on_nr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);
    [range_off_nr,pc_fit_off_nr] = linear_pc_range(an_off_nr,pc_off_nr,pc_window_length,pc_search_range,pc_window_step,pc_min_r2,pc_max_residual,[min_toggle max_toggle],pc_max_slope_difference);

    start_pc_merge_anpc = max([range_on_fr(1,:);range_off_fr(1,:);range_on_nr(1,:);range_off_nr(1,:)],[],1,'includemissing');
    end_pc_merge_anpc = min([range_on_fr(2,:);range_off_fr(2,:);range_on_nr(2,:);range_off_nr(2,:)],[],1,'includemissing');
    no_overlap = start_pc_merge_anpc >= end_pc_merge_anpc;
    start_pc_merge_anpc(no_overlap) = NaN;
    end_pc_merge_anpc(no_overlap) = NaN;

    % Find ON and OFF linear ranges independently, then intersect them.
    [range_fr_nr_on,fit_fr_nr_on] = linear_fr_nr_range(pc_on_nr,pc_on_fr,pc_window_length,pc_search_range,pc_window_step,fr_nr_min_r2,fr_nr_max_residual);
    [range_fr_nr_off,fit_fr_nr_off] = linear_fr_nr_range(pc_off_nr,pc_off_fr,pc_window_length,pc_search_range,pc_window_step,fr_nr_min_r2,fr_nr_max_residual);

    start_pc_merge_linear = max([range_fr_nr_on(1,:);range_fr_nr_off(1,:)],[],1,'includemissing');
    end_pc_merge_linear = min([range_fr_nr_on(2,:);range_fr_nr_off(2,:)],[],1,'includemissing');
    no_overlap = start_pc_merge_linear >= end_pc_merge_linear;
    start_pc_merge_linear(no_overlap) = NaN;
    end_pc_merge_linear(no_overlap) = NaN;

    start_pc_merge = start_pc_merge_linear;
    end_pc_merge = end_pc_merge_linear;

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
        sigma_o3_on_t(:,i) = interp1(BDM_t, sigma_o3_on, Ta(:,i));
        sigma_o3_off_t(:,i) = interp1(BDM_t, sigma_o3_off, Ta(:,i));
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

    if ~isempty(aerosol)
        [D_fr,D_nr,beta_a] = aerosol_correction(pc_off_fr,N_o3_fr_pc,z(:),C_fr,C_nr,...
            Na.*sigma_m_off,sigma_o3_off_t,dsigma_o3,lam_on,lam_off,aerosol);

        fr = isfinite(D_fr);
        nr = isfinite(D_nr);
        N_o3_fr_an(fr) = N_o3_fr_an(fr) + D_fr(fr);
        N_o3_fr_pc(fr) = N_o3_fr_pc(fr) + D_fr(fr);
        N_o3_fr_mg(fr) = N_o3_fr_mg(fr) + D_fr(fr);
        N_o3_nr_an(nr) = N_o3_nr_an(nr) + D_nr(nr);
        N_o3_nr_pc(nr) = N_o3_nr_pc(nr) + D_nr(nr);
        N_o3_nr_mg(nr) = N_o3_nr_mg(nr) + D_nr(nr);
    end

    % Merge near range and far range for analog calculated profile, photon counting calculated profile, and analog-pc calculated profile
    [N_o3_an_merged,w_o3_an] = merge_profiles(N_o3_nr_an,N_o3_fr_an,start_merge_id,end_merge_id);
    [N_o3_pc_merged,w_o3_pc] = merge_profiles(N_o3_nr_pc,N_o3_fr_pc,start_pc_merge,end_pc_merge);
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
    uncertainty_options.dz = dz;
    uncertainty_options.L = n_shots_total;
    uncertainty_options.uTau = uTau;
    uncertainty_options.dsigma = dsigma_o3;
    uncertainty_options.channelCorrelation = 'independent';
    uncertainty_options.sigma_o3_on = sigma_o3_on_t;
    uncertainty_options.sigma_o3_off = sigma_o3_off_t;
    uncertainty_options.sigma_m_on = sigma_m_on;
    uncertainty_options.sigma_m_off = sigma_m_off;
    uncertainty_options.u_sigma_o3_rel = u_sigma_o3_rel;
    uncertainty_options.u_sigma_m_rel = u_sigma_m_rel;

    uncertainty_options.C = C_fr;
    [uNO3_fr,pc_uncertainty_fr] = photon_count_uncertainty(pc_R_sum,pc_P_bg,N_o3_fr_pc,Na,pc_uB,pc_rb,on_fr_id,off_fr_id,uncertainty_options);

    uncertainty_options.C = C_nr;
    [uNO3_nr,pc_uncertainty_nr] = photon_count_uncertainty(pc_R_sum,pc_P_bg,N_o3_nr_pc,Na,pc_uB,pc_rb,on_nr_id,off_nr_id,uncertainty_options);

    [uNO3_pc_merged,pc_uncertainty_merged] = merge_pc_uncertainty(pc_uncertainty_nr,pc_uncertainty_fr,w_o3_pc);

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

    % PC merge intervals and analog-PC fit diagnostics.
    out.start_pc_merge = start_pc_merge;
    out.end_pc_merge = end_pc_merge;
    out.start_pc_merge_anpc = start_pc_merge_anpc;
    out.end_pc_merge_anpc = end_pc_merge_anpc;
    out.start_pc_merge_linear = start_pc_merge_linear;
    out.end_pc_merge_linear = end_pc_merge_linear;
    out.range_fr_nr_on = range_fr_nr_on;
    out.range_fr_nr_off = range_fr_nr_off;
    out.fit_fr_nr_on = fit_fr_nr_on;
    out.fit_fr_nr_off = fit_fr_nr_off;
    out.range_on_fr = range_on_fr;
    out.range_off_fr = range_off_fr;
    out.range_on_nr = range_on_nr;
    out.range_off_nr = range_off_nr;
    out.pc_fit_on_fr = pc_fit_on_fr;
    out.pc_fit_off_fr = pc_fit_off_fr;
    out.pc_fit_on_nr = pc_fit_on_nr;
    out.pc_fit_off_nr = pc_fit_off_nr;

    % Rayleigh extinction correction
    out.N_o3_dm = N_o3_dm;

    if ~isempty(aerosol)
        out.aerosol = struct('settings',aerosol,'beta_off',beta_a,...
            'delta_N_o3_fr',D_fr,'delta_N_o3_nr',D_nr);
    end

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

    % PC uncertainty
    out.uNO3_fr = uNO3_fr;
    out.uNO3_nr = uNO3_nr;
    out.uNO3_pc_merged = uNO3_pc_merged;
    out.uqO3_fr = uqO3_fr;
    out.uqO3_nr = uqO3_nr;
    out.uqO3_pc_merged = uqO3_pc_merged;
    out.pc_uncertainty_fr = pc_uncertainty_fr;
    out.pc_uncertainty_nr = pc_uncertainty_nr;
    out.pc_uncertainty_merged = pc_uncertainty_merged;
    out.pc_uncertainty_metadata = struct('components',{{'DET','BKG','SAT','O3','RAYLEIGH'}},...
        'complete_atmospheric_budget',false,'tau_seconds',tau_d,'uTau_seconds',uTau,...
        'uTau_relative',uTau_rel,'dead_time_uncertainty_source','provisional engineering estimate',...
        'channel_correlation','independent','formula_mode','literal_paper',...
        'averaging_method','raw_counts_accumulation','signal_normalization','total_counts',...
        'cross_section_convention','single_path',...
        'ozone_cross_section_correlation','independent','rayleigh_cross_section_correlation','full',...
        'u_sigma_o3_relative',u_sigma_o3_rel,'u_sigma_m_relative',u_sigma_m_rel,...
        'cross_section_uncertainty_source','Leblanc 2016 Table 3 working estimates');

    save(output_path, '-struct', 'out')
end

