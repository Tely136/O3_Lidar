function u = uNO3_Na2_NC(dS_O3,dS_M,dS_O2,q_O2,Na,pa,Ta,u_pa,u_Ta)
% Eq.86: independent temperature/pressure; NO2/SO2 NUMBER densities.
% Pressure and its SD use identical units; temperature and its SD in kelvin.
    u=abs((dS_M+q_O2.*dS_O2)./dS_O3).*Na.*hypot(u_pa./pa,u_Ta./Ta);
end
