function pc = scale_binary_pc(data_in, binwidth, n_acq)
    % Function to convert binary photon counting data to physical units

    mhz = 1e-6 .* physconst("LightSpeed") ./ (2*binwidth);
    pc = double(data_in) .* mhz ./ n_acq;
end