function u = uNO3_dS_M_R(Na,dS_O3,udS_M_on,udS_M_off)
% Eq.57 elastic backscatter. Na is AIR density; denominator is OZONE
% differential, NOT Rayleigh differential. All cross sections single-path m^2.
% Paper DeltaSigma_O3=2*dS_O3, cancelling its factor 2. Update old callers.
    u=Na./abs(dS_O3).*hypot(udS_M_on,udS_M_off);
end
