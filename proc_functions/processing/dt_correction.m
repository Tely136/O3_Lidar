function P = dt_correction(R,r,fd)
    P = R./(1-r/fd);
end