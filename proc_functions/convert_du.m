function o3_col = convert_du(NO3,prof_alt,alt,alt_bins)

    n_alt = length(alt);
    o3_col = NaN(n_alt,1);

    for i = 1:n_alt
        lower_alt = alt_bins(i,2);
        upper_alt = alt_bins(i,1);

        idx = prof_alt >= lower_alt & prof_alt < upper_alt;
        o3_col(i) = trapz(prof_alt(idx) .* 1e3,NO3(idx));
    end

    o3_col = o3_col ./ (2.68678e20);
end