function data_out = scale_data(mode,data_in,an_range,adc_bits,mhz,n_acq)
    %% Function to convert input binary data to physical units

    % Analog
    if mode == 1
        data_out = scale_binary_analog(data_in,an_range,adc_bits,n_acq);
    % Photon counting
    elseif mode == 2
        data_out = scale_binary_pc(data_in, mhz, n_acq);
    end
end