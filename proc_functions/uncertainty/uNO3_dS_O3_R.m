function u = uNO3_dS_O3_R(N,dS_O3,udS_O3_on,udS_O3_off)
% Eq.47 elastic backscatter. N is OZONE m^-3; cross sections single-path m^2.
% Paper DeltaSigma=2*dS_O3, cancelling the factor 2 in Eq.47.
% Uncertainty inputs are absolute cross-section SDs, not relative fractions.
% Preserve printed N: a negative retrieved N gives a signed formula output.
    u=N./abs(dS_O3).*hypot(udS_O3_on,udS_O3_off);
end
