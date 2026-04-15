function an = scale_binary_analog(data_in,an_range,adc_bits,n_acq)
    %% Function to convert binary analog data to physical units
    
    switch an_range
        case 0
            mv = 500;
        case 1
            mv = 100;
        case 2
            mv = 20;
    end
    an = double(data_in) .* mv ./(2^adc_bits-1) ./ n_acq;
end