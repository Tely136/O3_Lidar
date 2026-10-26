function u = uNO3_Na2_FC(dS_O3,dS_M,dS_O2,q_O2,Na,pa,Ta,u_pa,u_Ta)
% Eq.89: positively fully correlated pressure/temperature (rho=+1).
% NO2/SO2 specified as NUMBER densities; cross sections single-path m^2.
    u=abs((dS_M+q_O2.*dS_O2)./dS_O3).*Na.*abs(u_pa./pa-u_Ta./Ta);
end
