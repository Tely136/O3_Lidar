function sigma_m = rayleigh_cross_section(lambda)
    % input lambda on and off in m and air number density

    % refractive index of air 
    n_air = @(lam) 1 + (1e-8) * (6432.8 + 2949810/(146-lam.^-2) + 25540/(41-lam.^-2));
    n = n_air(lambda * 1e6); % lambda in um
    
    % Rayleigh scattering cross-section
    rho = 0.035;
    Ns_cm3 = 2.5469e19;

    s = @ (lam, n) ...
        (24 * pi^3) ./ (lam.^4 .* Ns_cm3.^2) .* ...
        ((n.^2 - 1)./(n.^2 + 2)).^2 * ((6 + 3*rho)/(6 - 7*rho));
    
    sigma_m = s(lambda * 1e2, n) .* 1e-4; % lambda in cm and Nd in cm^-3
end