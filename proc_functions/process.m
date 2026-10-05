function out = process(input_path,output_path)
    % DIAL Processing Functions
    addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\O3_Lidar\proc_functions"));
    
    % Signal Analysis and Modeling Functions
    addpath(genpath("C:\Users\tely1\OneDrive\Desktop\EEJ9903\signal_analysis\functions"));

    % parameters to be moved into config file
    data_format = "bin";
    prefix = "ol";
    n_avg = 30;
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

    % Alt vector
    z = (1:n_bins) .* dz;
    z_km = z.*1e-3;

    % Use midpoint time to represent time
    % TODO: make use of start and end times
    times = times(:,2);

    % time-average data, this doesn't handle discontinuous data
    n_prof_avg = floor(n_prof / n_avg);
    mod(n_prof,n_avg); % can check if remainder is >= n_avg/2 and add another profile

    pc_r = scale_binary_pc(pc_R,3.75,1201);
    pc_P_temp = dt_correction(pc_R,pc_r,fd);

    an_B = bg_correction(z,an_R,7000);
    [pc_B,~,uB,rb] = bg_correction(z,pc_P_temp,7000);

    pc_P = pc_P_temp - pc_B; 
    an_bin = an_R - an_B;



    an_bin_avg = NaN(n_bins,n_prof_avg,4);
    pc_P_avg = NaN(n_bins,n_prof_avg,4);

    times_avg = NaT(n_prof_avg,1);
    for i = 1:n_prof_avg
        ids = (i-1)*n_avg+1:(i-1)*n_avg+n_avg;

        an_bin_avg(:,i,:) = mean(an_bin(:,ids,:),2);
        pc_P_avg(:,i,:) = mean(pc_P(:,ids,:),2);

        times_avg(i) = mean(times(ids),1);
    end

    pc_P_avg(pc_P_avg<=0) = NaN;
    an_bin_avg(an_bin_avg<=0) = NaN;

    pc_p_avg = scale_binary_pc(pc_P_avg,3.75,1201);
    an_mV_avg = scale_binary_analog(an_bin_avg,0,16,1201); % TODO loop over devices and scale each one according to config values


    % Grid
    z_grid = repmat((1:n_bins)',1,n_prof_avg) .* dz;
    
    [mg_on_fr,coeff_on_fr,r2_on_fr] = an_pc_merge(an_mV_avg(:,:,1),pc_p_avg(:,:,1),window_length_fr,min_toggle,max_toggle);
    [mg_off_fr,coeff_off_fr,r2_off_fr] = an_pc_merge(an_mV_avg(:,:,2),pc_p_avg(:,:,2),window_length_fr,min_toggle,max_toggle);
    [mg_on_nr,coeff_on_nr,r2_on_nr] = an_pc_merge(an_mV_avg(:,:,3),pc_p_avg(:,:,3),window_length_nr,min_toggle,max_toggle);
    [mg_off_nr,coeff_off_nr,r2_off_nr] = an_pc_merge(an_mV_avg(:,:,4),pc_p_avg(:,:,4),window_length_nr,min_toggle,max_toggle);

    an_on_fr = an_mV_avg(:,:,1);
    an_off_fr = an_mV_avg(:,:,2);
    an_on_nr = an_mV_avg(:,:,3);
    an_off_nr = an_mV_avg(:,:,4);

    pc_on_fr = pc_p_avg(:,:,1);
    pc_off_fr = pc_p_avg(:,:,2);
    pc_on_nr = pc_p_avg(:,:,3);
    pc_off_nr = pc_p_avg(:,:,4);

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
    dsigma_o3_temp = load('o3_cross_287_299.mat');
    dsigma_o3 = NaN(n_bins,n_prof_avg);
    for i=1:n_prof_avg
        dsigma_o3(:,i) = interp1(dsigma_o3_temp.dsigma_temp,dsigma_o3_temp.dsigma,Ta(:,i));
    end

    % Differential Rayleigh extinction correction
    dsigma_m = rayleigh_cross_section(lam_on) - rayleigh_cross_section(lam_off);
    N_o3_dm = rayleigh_correction(Na,dsigma_m,dsigma_o3);

    N_o3_fr_an = dial(an_on_fr,an_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_an = dial(an_on_nr,an_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    N_o3_fr_pc = dial(pc_on_fr,pc_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_pc = dial(pc_on_nr,pc_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;

    N_o3_fr_mg = dial(mg_on_fr,mg_off_fr,C_fr,dz,dsigma_o3) - N_o3_dm;
    N_o3_nr_mg = dial(mg_on_nr,mg_off_nr,C_nr,dz,dsigma_o3) - N_o3_dm;


    q_o3_an_fr = (N_o3_fr_an ./ Na) * 1e9;
    q_o3_an_nr = (N_o3_nr_an ./ Na) * 1e9;

    q_o3_pc_fr = (N_o3_fr_pc ./ Na) * 1e9;
    q_o3_pc_nr = (N_o3_nr_pc ./ Na) * 1e9;

    q_o3_mg_fr = (N_o3_fr_mg ./ Na) * 1e9;
    q_o3_mg_nr = (N_o3_nr_mg ./ Na) * 1e9;

    [q_o3_an_merged,w_o3_an] = merge_profiles(q_o3_an_nr,q_o3_an_fr,z_km,start_merge,end_merge);
    [q_o3_pc_merged,w_o3_pc] = merge_profiles(q_o3_pc_nr,q_o3_pc_fr,z_km,start_merge,end_merge);
    [q_o3_mg_merged,w_o3_mg] = merge_profiles(q_o3_mg_nr,q_o3_mg_fr,z_km,start_merge,end_merge);


    % Uncertainty
    uTau = 0;
    % uNO3_fr = uncertainty(1,2,pc_P,R_avg,uB,rb,C_fr,dsigma_o3,uTau,dz,1200);
    % uNO3_nr = uncertainty(3,4,pc_P,R_avg,uB,rb,C_fr,dsigma_o3,uTau,dz,1200);
    % 
    % uqO3_fr = uNO3_fr ./ Na .* 1e9;
    % uqO3_nr = uNO3_nr ./ Na .* 1e9;

    out = struct();

    out.q_o3_an_fr = q_o3_an_fr;
    out.q_o3_an_nr = q_o3_an_nr;
    out.q_o3_pc_fr = q_o3_pc_fr;
    out.q_o3_pc_nr = q_o3_pc_nr;
    out.q_o3_mg_fr = q_o3_mg_fr;
    out.q_o3_mg_nr = q_o3_mg_nr;

    out.N_o3_fr_an = N_o3_fr_an;
    out.N_o3_nr_an = N_o3_nr_an;
    out.N_o3_fr_pc = N_o3_fr_pc;
    out.N_o3_nr_pc = N_o3_nr_pc;
    out.N_o3_fr_mg = N_o3_fr_mg;
    out.N_o3_nr_mg = N_o3_nr_mg;

    out.q_o3_an_merged = q_o3_an_merged;
    out.w_o3_an = w_o3_an;
    out.q_o3_pc_merged = q_o3_pc_merged;
    out.w_o3_pc = w_o3_pc;
    out.q_o3_mg_merged = q_o3_mg_merged;
    out.w_o3_mg = w_o3_mg;

    out.coeff_on_fr = coeff_on_fr;
    out.coeff_off_fr = coeff_off_fr;
    out.coeff_on_nr = coeff_on_nr;
    out.coeff_off_nr = coeff_off_nr;

    out.r2_on_fr = r2_on_fr;
    out.r2_off_fr = r2_off_fr;
    out.r2_on_nr = r2_on_nr;
    out.r2_off_nr = r2_off_nr;

    out.z_vec = z;
    out.times = times;
    out.times_avg = times_avg;

    out.configs = configs;
    out.C_fr = C_fr;
    out.C_nr = C_nr;

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

