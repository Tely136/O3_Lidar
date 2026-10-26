function u = uNO3_Na1(dS_O3,dS_M,dS_O2,q_O2,u_Na)
% Eq.81: NO2/SO2 specified as NUMBER densities. All differentials must
% use the same convention (repo: single-path m^2).
    u=abs((dS_M+q_O2.*dS_O2)./dS_O3).*u_Na;
end
