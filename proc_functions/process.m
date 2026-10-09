function out = process(input_path,output_path)
    % DIAL Processing Functions
    addpath(genpath("..\O3_Lidar\proc_functions"));
    
    % Signal Analysis and Modeling Functions
    addpath(genpath("..\signal_analysis\functions"));

    % parameters to be moved into config file
    data_format = "bin";
    prefix = "ol";
    n_avg = 30;
    n_shots = 1200;
    fd = 270; % fix the naming convention around this
    start_merge = 1.5;
    end_merge = 2;
    min_toggle = 3;
    max_toggle = 60;

    % Bin width
    dz = 3.75;

    window_length_fr = 500;
    window_length_nr = 100;

    [an_R,pc_R,times,configs] = read_lidar_datafiles(input_path,data_format,prefix);
    n_bins = size(an_R,1);
    n_prof = size(an_R,2);

    n_prof_avg = floor(n_prof / n_avg);


    % Alt vector
    z = (1:n_bins) .* dz;
    z_km = z.*1e-3;

    [an_bin_avg, pc_P_avg,times_avg,an_B_avg,pc_B_avg] = preprocess(an_R,pc_R,times,dz,n_avg,n_shots,fd);
    an_bin_bg = an_bin_avg - an_B_avg;
    pc_P_bg = pc_P_avg - pc_B_avg;

    an_bin_bg(an_bin_bg <= 0) = NaN;
    pc_P_bg(pc_P_bg <= 0) = NaN;

    an_mV_avg = scale_binary_analog(an_bin_bg,0,16,1200); % TODO loop over devices and scale each one according to config values
    an_mV_bg_avg = scale_binary_analog(an_B_avg,0,16,1200);

    pc_p_avg = scale_binary_pc(pc_P_bg,3.75,1200);
    pc_p_bg_avg = scale_binary_pc(pc_B_avg,3.75,1200);

    % Grid
    z_grid = repmat((1:n_bins)',1,n_prof_avg) .* dz;

    % Analog
    an_on_fr = an_mV_avg(:,:,1);
    an_off_fr = an_mV_avg(:,:,2);
    an_on_nr = an_mV_avg(:,:,3);
    an_off_nr = an_mV_avg(:,:,4);

    % PC
    pc_on_fr = pc_p_avg(:,:,1);
    pc_off_fr = pc_p_avg(:,:,2);
    pc_on_nr = pc_p_avg(:,:,3);
    pc_off_nr = pc_p_avg(:,:,4);
    
    [mg_on_fr,coeff_on_fr,r2_on_fr] = an_pc_merge(an_on_fr, pc_on_fr, window_length_fr, min_toggle, max_toggle);
    [mg_off_fr,coeff_off_fr,r2_off_fr] = an_pc_merge(an_off_fr, pc_off_fr, window_length_fr, min_toggle, max_toggle);

    [mg_on_nr,coeff_on_nr,r2_on_nr] = an_pc_merge(an_on_nr, pc_on_nr, window_length_nr, min_toggle, max_toggle);
    [mg_off_nr,coeff_off_nr,r2_off_nr] = an_pc_merge(an_off_nr, pc_off_nr, window_length_nr, min_toggle, max_toggle);

    % On and off wavelengths
    lam_on = 287.2 * 1e-9; % m
    lam_off = 299.1 * 1e-9; % m

    % FR and NR filter widths
    fr_fl_M1 = 11; fr_fl_M2 = 351; fr_fl_h1 = 1.2; fr_fl_h2 = 10;
    nr_fl_M1 = 11; nr_fl_M2 = 121; nr_fl_h1 = .2; nr_fl_h2 = 1.2;

    % FR and NR filters
    C_fr = filter_gen(fr_fl_M1,fr_fl_M2,fr_fl_h1,fr_fl_h2,z_km);
    C_nr = filter_gen(nr_fl_M1,nr_fl_M2,nr_fl_h1,nr_fl_h2,z_km);

    % Air number density and temperature
    Na = airNumberDensity(z_grid);
    Ta = atmosisa(z_grid);

    % Differential ozone absorption cross section as function of temperature
    % dsigma_o3_temp = load('o3_cross_287_299.mat');
    % dsigma_o3 = NaN(n_bins,n_prof_avg);
    % for i=1:n_prof_avg
    %     dsigma_o3(:,i) = interp1(dsigma_o3_temp.dsigma_temp,dsigma_o3_temp.dsigma,Ta(:,i));
    % end

    BDM = load("cross_sections\ozone\bdm\bdm_287_299.mat");
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

    % Total Rayleigh extinction
    % alpha_m_on = sigma_m_on .* Na;
    % alpha_m_off = sigma_m_off .* Na;
    % 
    % Bs_m_on = alpha_m_on .* rayleighPhaseFunction(pi) ./ (4*pi);
    % Bs_m_off = alpha_m_off .* rayleighPhaseFunction(pi) ./ (4*pi);
    % 
    % N_o3_Bsm_fr = dial(Bs_m_on,Bs_m_off,C_fr,dz,dsigma_o3);
    % N_o3_Bsm_nr = dial(Bs_m_on,Bs_m_off,C_nr,dz,dsigma_o3);


    N_o3_fr_an = dial(an_on_fr,an_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_an = dial(an_on_nr,an_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    N_o3_fr_pc = dial(pc_on_fr,pc_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_pc = dial(pc_on_nr,pc_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    N_o3_fr_mg = dial(mg_on_fr,mg_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_mg = dial(mg_on_nr,mg_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    [N_o3_an_merged,w_o3_an] = merge_profiles(N_o3_nr_an,N_o3_fr_an,z_km,start_merge,end_merge);
    [N_o3_pc_merged,w_o3_pc] = merge_profiles(N_o3_nr_pc,N_o3_fr_pc,z_km,start_merge,end_merge);
    [N_o3_mg_merged,w_o3_mg] = merge_profiles(N_o3_nr_mg,N_o3_fr_mg,z_km,start_merge,end_merge);

    q_o3_an_fr = (N_o3_fr_an ./ Na) * 1e9;
    q_o3_an_nr = (N_o3_nr_an ./ Na) * 1e9;

    q_o3_pc_fr = (N_o3_fr_pc ./ Na) * 1e9;
    q_o3_pc_nr = (N_o3_nr_pc ./ Na) * 1e9;

    q_o3_mg_fr = (N_o3_fr_mg ./ Na) * 1e9;
    q_o3_mg_nr = (N_o3_nr_mg ./ Na) * 1e9;

    q_o3_an_merged = (N_o3_an_merged ./ Na) * 1e9;
    q_o3_pc_merged = (N_o3_pc_merged ./ Na) * 1e9;
    q_o3_mg_merged = (N_o3_mg_merged ./ Na) * 1e9;

    % Uncertainty
    % uTau = 0;
    % uNO3_fr = uncertainty(1,2,pc_P,R_avg,uB,rb,C_fr,dsigma_o3,uTau,dz,1200);
    % uNO3_nr = uncertainty(3,4,pc_P,R_avg,uB,rb,C_fr,dsigma_o3,uTau,dz,1200);
    % 
    % uqO3_fr = uNO3_fr ./ Na .* 1e9;
    % uqO3_nr = uNO3_nr ./ Na .* 1e9;

    out = struct();

    % Analog signal
    out.an_on_fr = an_on_fr;
    out.an_off_fr = an_off_fr;
    out.an_on_nr = an_on_nr;
    out.an_off_nr = an_off_nr; % add background for each channel as well

    % PC signal
    out.pc_on_fr = pc_on_fr;
    out.pc_off_fr = pc_off_fr;
    out.pc_on_nr = pc_on_nr;
    out.pc_off_nr = pc_off_nr; % add background here too

    % AN-PC Merge signal
    out.mg_on_fr = mg_on_fr;
    out.mg_off_fr = mg_off_fr;
    out.mg_on_nr = mg_on_nr;
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

    % AN-PC Merge coefficients
    out.coeff_on_fr = coeff_on_fr;
    out.coeff_off_fr = coeff_off_fr;
    out.coeff_on_nr = coeff_on_nr;
    out.coeff_off_nr = coeff_off_nr;

    % AN-PC Merge r^2
    out.r2_on_fr = r2_on_fr;
    out.r2_off_fr = r2_off_fr;
    out.r2_on_nr = r2_on_nr;
    out.r2_off_nr = r2_off_nr;

    % Altitude, raw time, averaged time
    out.z_vec = z;
    out.times = times;
    out.times_avg = times_avg;

    % Licel configs
    out.configs = configs;

    % SG Filter coefficients
    out.C_fr = C_fr;
    out.C_nr = C_nr;

    % Air number density and temperature
    out.Na = Na;
    out.Ta = Ta;

    % Ozone extinction correction
    out.sigma_o3_on = sigma_o3_on_t;
    out.sigma_o3_off = sigma_o3_off_t;
    out.dsigma_o3 = dsigma_o3;

    % out.uNO3_fr = uNO3_fr;
    % out.uNO3_nr = uNO3_nr;
    % 
    % out.uqO3_fr = uqO3_fr;
    % out.uqO3_nr = uqO3_nr;

    save(output_path, '-struct', 'out')
end

function u = uncertainty(id_on,id_off,P,R,uB,rb,C,dsigma_o3,uTau,dz,L)
    P_on = P(:,:,id_on);
    P_off = P(:,:,id_off);

    R_on = R(:,:,id_on);
    R_off = R(:,:,id_off);

    uB_on = uB(:,:,id_on);
    uB_off = uB(:,:,id_off);

    rb_on = rb(:,id_on);
    rb_off = rb(:,id_off);

    n_bins = size(P_on,1);
    n_prof = size(P_on,2);

    z = (1:n_bins).*dz;

    u = NaN(n_bins,n_prof);
    for i = 1:n_prof
        u(:,i) = uNO3(P_on(:,i),P_off(:,i),R_on(:,i),R_off(:,i),dz,L,uTau,dsigma_o3(:,i),C,uB_on(:,i),uB_off(:,i),rb_on(i),rb_off(i),z);
    end
end

