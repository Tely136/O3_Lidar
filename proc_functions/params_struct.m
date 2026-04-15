function params = params_struct(lam_on,lam_off,n_channels,min_avg,bg_bins,td,dz,min_toggle,max_toggle,M1,M2,h1,h2)
    params = struct;
    params.lam_on = lam_on;
    params.lam_off = lam_off;
    params.n_channels = n_channels;
    params.min_avg = min_avg;
    params.bg_bins = bg_bins;
    params.td = td;
    params.dz = dz;
    params.min_toggle = min_toggle;
    params.max_toggle = max_toggle;
    params.M1 = M1;
    params.M2 = M2;
    params.h1 = h1;
    params.h2 = h2;
end