function ds = read_glass(file)
    time = hdfread(file,"DATETIME");
    time_start = hdfread(file,"DATETIME.START");
    time_stop = hdfread(file,"DATETIME.STOP");
    
    ds.time = mj2k(time);
    ds.time_start = mj2k(time_start);
    ds.time_stop = mj2k(time_stop);
    
    ds.alt = double(hdfread(file,"ALTITUDE"));
    ds.alt_inst = hdfread(file,"ALTITUDE.INSTRUMENT");
    
    ds.NO3     = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL");
    ds.u_NO3_R = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD");
    ds.u_NO3_S = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD");
    ds.u_NO3_C = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD");
    
    ds.qO3     = hdfread(file,"O3.MIXING.RATIO.VOLUME_DERIVED");
    ds.u_qO3_R = hdfread(file,"O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD");
    ds.u_qO3_S = hdfread(file,"O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD");
    ds.u_qO3_C = hdfread(file,"O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD");
    
    ds.resolution_fwhm = hdfread(file,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM");
    
    ds.p = hdfread(file,"PRESSURE_INDEPENDENT"); % hPa
    ds.t = hdfread(file,"TEMPERATURE_INDEPENDENT"); % k
end


function out = mj2k(in)
    out = datetime(round(in .* 24 .* 60 .* 60),'ConvertFrom','epochtime','Epoch','2000-01-01', 'TimeZone','UTC');
end