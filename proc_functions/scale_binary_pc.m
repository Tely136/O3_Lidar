function pc = scale_binary_pc(data_in, mhz, n_acq)
    % Function to convert binary photon counting data to physical units
    pc = double(data_in) .* mhz ./ n_acq;
end