function output = O3_quicklook(data,times,params)
    % TODO: put wavelengths in configs and read it from there (maybe)
    hkm = data.hkm;

    pre_output_on = preprocess(data.an_on,data.pc_on,times,params);
    pre_output_off = preprocess(data.an_off,data.pc_off,times,params);

    % Determine filter widths and load filter coefficients
    M1 = params.M1;
    M2 = params.M2;
    h1 = params.h1;
    h2 = params.h2;
    
    fl = gen_framelength(M1,M2,h1,h2,hkm);
    
    % TODO: don't load filters here, make altitude dependent filter list
    % previously and input to processsing function
    sg = load('sg_filters.mat');
    sg_diff = sg.sg_diff;

    C = cell(1,length(fl));
    for k = 1:length(fl)
        m = fl(k);
        N = (m-1)/2;

        C{k} = sg_diff{N};
    end
    
    % Calculate S
    r = log(pre_output_off.glued./pre_output_on.glued);
    S = smooth_derivative(r,C,params.dz);
    n_avg = size(S,2);
    
    % Calculate differential ozone absorption cross section
    [t,~,ND] = standard_atmosphere(hkm.*1e3); % ND in m^-3
    
    % Differential ozone absorption cross section as function of temperature
    % TODO: make O3 cross section input better
    dsigma_o3_temp = load('o3_cross_287_299.mat');
    dsigma_o3 = interp1(dsigma_o3_temp.dsigma_temp,dsigma_o3_temp.dsigma,t);
    % TODO: maybe don't repmat here
    dsigma_o3_2d = repmat(dsigma_o3,1,n_avg);
    
    % Calculate differential Rayleigh extinction cross section
    sigma_m_on = rayleigh_cross_section(params.lam_on);   
    sigma_m_off = rayleigh_cross_section(params.lam_off);   
    dsigma_m = sigma_m_on - sigma_m_off;

    pre_output_on.sigma_m = sigma_m_on;
    pre_output_off.sigma_m = sigma_m_off;

    % Calculate NO3 and qO3
    ND_2d = repmat(ND,1,n_avg);

    D_No3 = dsigma_m.*ND_2d./dsigma_o3_2d;
    D_qo3 = dsigma_m./dsigma_o3_2d;

    No3 = S./dsigma_o3_2d./2;
    qo3 = (S./dsigma_o3_2d./ND_2d./2) .*10^9;

    output = struct;
    output.on = pre_output_on;
    output.off = pre_output_off;

    output.No3 = No3;
    output.qo3 = qo3 - D_qo3;
    output.DNo3 = D_No3;
    output.Dqo3 = D_qo3.*10^9;

    output.dso3 = dsigma_o3;
    output.ND = ND;
    output.hkm = hkm;
    output.C = C;
end

function [output] = preprocess(an,pc,times,params)
    min_avg = params.min_avg;
    bg_bins = params.bg_bins;
    td = params.td;
    min_toggle = params.min_toggle;
    max_toggle = params.max_toggle;

    % Remove data before gate
    [an,pc_avg] = rm_gate(an,pc);

    % Dead-time correction of PC data
    pc_dt = correct_deadtime(pc,td);

    % Remove background
    [an_bg_rem,an_bg] = remove_background(an,bg_bins,"simple");
    [~,pc_bg] = remove_background(pc,bg_bins,"simple");  
    pc_dt_bg_rem = pc_dt - pc_bg;

    an_bg_rem(an_bg_rem<0 | pc_dt_bg_rem<0) = NaN;
    pc_dt_bg_rem(an_bg_rem<0 | pc_dt_bg_rem<0) = NaN;

    % Glue analog and PC data for each TR
    [glued, coeffs, corrcoefs, valid] = glue_an_pc(an_bg_rem,pc_dt_bg_rem,min_toggle,max_toggle); 

    % Time average
    [glued,times_avg,times_counts] = retime_avg(glued,times,min_avg);

    output = struct;

    output.an = an;
    output.pc = pc;

    output.an_avg = an;
    output.pc_avg = pc_avg;

    output.pc_dt = pc_dt;

    output.an_bg_rem = an_bg_rem;
    output.pc_dt_bg_rem = pc_dt_bg_rem;

    output.an_bg = an_bg;
    output.pc_bg = pc_bg;

    output.glued = glued;

    output.glued_coeffs = coeffs;
    output.glued_valid = valid;
    output.glued_corrcoefs = corrcoefs;

    output.times_avg = times_avg;
    output.times_counts = times_counts;
end

function [data_avg,times_avg,times_avg_counts] = retime_avg(data,times,min_avg)
    tt = timetable(times,data');
    % TODO: look into way of using retime similar to how it was in Dingdong's
    % code
    tt_mean = retime(tt,'regular','mean','TimeStep',minutes(min_avg));
    tt_count = retime(tt,'regular','count','TimeStep',minutes(min_avg));
    data_avg = tt_mean.Var1';
    times_avg = tt_mean.times;
    times_avg_counts = tt_count.Var1(:,1);
end

function [data_bg_rem,B] = remove_background(data,bg_bins,mode)
    B = bkg(data,bg_bins,mode);

    [r,~] = size(data);

    data_bg_rem = data - repmat(B,r,1);
end

function data = correct_deadtime(data,rd)
    data = data./(1-data/rd);
end

function sigma_m = rayleigh_cross_section(lambda)
    % input lambda on and off in m and air number density

    % refractive index of air 
    n_air = @(lam) 1 + (1e-8) * (6432.8 + 2949810/(146-lam.^-2) + 25540/(41-lam.^-2));
    n = n_air(lambda * 1e6); % lambda in um
    
    % Rayleigh scattering cross-section
    rho = 0.035;
    Ns_cm3 = 2.5469e19;

    s = @ (lam, n) ...
        (24 * pi^3) ./ (lam.^4 .* Ns_cm3.^2) .* ...
        ((n.^2 - 1)./(n.^2 + 2)).^2 * ((6 + 3*rho)/(6 - 7*rho)) .* 1e-4;
    
    sigma_m = s(lambda * 1e2, n); % lambda in cm and Nd in cm^-3
end

function [an_out,pc_out] = rm_gate(an,pc)
    [~,an_gate_id] = max(an,[],1);
    for i = 1:size(an,2)
        an(1:an_gate_id(i)-1,i) = NaN;
        pc(1:an_gate_id(i)-1,i) = NaN;
    end
    an_out = an;
    pc_out = pc;
end

function B = bkg(data,bins,mode)
    switch mode
        case "simple"
            B = mean(data(bins,:),1);

        case "linear"
            % will be added
    end
end