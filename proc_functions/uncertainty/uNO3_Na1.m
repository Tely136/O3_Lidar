function u = uNO3_Na1(dS_O3,dS_M,dS_O2,q_O2,u_Na)
    % Eq. 81: interfering gases specified as number densities.

    u = abs((dS_M + q_O2 .* dS_O2) ./ dS_O3) .* u_Na;
end
