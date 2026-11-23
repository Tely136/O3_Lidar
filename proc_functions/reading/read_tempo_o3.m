function ds = read_tempo_o3(file)
    ds.time = ncread(file, "/geolocation/time");

    ds.lat = ncread(file, "/geolocation/latitude");
    ds.lon = ncread(file, "/geolocation/longitude");

    ds.raa = ncread(file, "/geolocation/relative_azimuth_angle");
    ds.sza = ncread(file, "/geolocation/solar_zenith_angle");
    ds.vza = ncread(file, "/geolocation/viewing_zenith_angle");

    ds.o3_prof = ncread(file, "/product/ozone_profile");
    ds.o3_prof_error = ncread(file, "/product/ozone_profile_error");
    ds.o3_prof_prec = ncread(file, "/product/ozone_profile_precision");

    ds.o3_strat = ncread(file, "/product/stratosphere_ozone_column");
    ds.o3_strat_error = ncread(file, "/product/stratosphere_ozone_error");
    ds.o3_strat_prec = ncread(file, "/product/stratosphere_ozone_precision");

    ds.o3_trop = ncread(file, "/product/troposphere_ozone_column");
    ds.o3_trop_error = ncread(file, "/product/troposphere_ozone_error");
    ds.o3_trop_prec = ncread(file, "/product/troposphere_ozone_precision");

    ds.o3_col = ncread(file, "/product/total_ozone_column");
    ds.o3_col_error = ncread(file, "/product/total_ozone_error");
    ds.o3_col_prec = ncread(file, "/product/total_ozone_precision");

    ds.cld_flag = ncread(file, "/support_data/cloud_flag");
    ds.cld_frac = ncread(file, "/support_data/eff_cloud_fraction");
    ds.gnd_pix_quality = ncread(file, "/support_data/ground_pixel_quality_flag");

    ds.o3_prof_apriori = ncread(file, "/support_data/ozone_apriori_profile");
    ds.o3_prof_apriori_error = ncread(file, "/support_data/ozone_apriori_profile_error");
    ds.o3_avg_kernel = ncread(file, "/support_data/ozone_averaging_kernel");

    ds.o3_prof_alt = ncread(file, "/support_data/ozone_profile_altitude");
    ds.o3_prof_alt_bound = ncread(file, "/support_data/ozone_profile_altitude_bounds");
    ds.o3_prof_temperature = ncread(file, "/support_data/ozone_profile_temperature");
end