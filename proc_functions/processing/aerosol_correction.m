function [D_fr,D_nr,beta_a] = aerosol_correction(P_off,N_o3,z,C_fr,C_nr,...
        alpha_m,sigma_off,dsigma,lam_on,lam_off,aerosol)
    lidar_ratio = 60; % sr
    angstrom = 0.5;
    beta_m = alpha_m .* (3/(8*pi));
    aerosol_scale = (lam_off/lam_on)^angstrom;
    backscatter_scale = aerosol_scale / (lam_off/lam_on)^4;
    dz = z(2) - z(1);
    [~,first] = min(abs(z-aerosol.start_km*1e3));
    [~,ref] = min(abs(z-aerosol.reference_km*1e3));

    D_fr = NaN(size(N_o3));
    D_nr = NaN(size(N_o3));
    beta_a = NaN(size(N_o3));

    for i = 1:size(N_o3,2)
        ids = (first:ref)';
        valid = isfinite(P_off(ids,i)) & P_off(ids,i) > 0 & isfinite(N_o3(ids,i));
        last_gap = max([0; find(~valid)]);
        ids(1:last_gap) = [];
        if numel(ids) < 2
            continue
        end

        signal = log(P_off(ids,i) .* z(ids).^2);
        molecular = beta_m(ids,i);
        baseline = N_o3(ids,i);
        ozone = baseline;
        beta = zeros(size(ozone));
        beta(end) = aerosol.reference_beta;
        differential_cross_section = dsigma(ids,i);

        for iteration = 1:20
            extinction = alpha_m(ids,i) + sigma_off(ids,i).*ozone;
            step = -diff(signal) - dz*(extinction(1:end-1) + extinction(2:end));
            [beta,converged] = invert_backscatter(step,molecular,beta,lidar_ratio,dz);
            if ~converged
                break
            end

            % Kuang et al. (2011), Eqs. 4–5, with Rayleigh already removed.
            B = beta ./ molecular;
            log_ratio = log1p(backscatter_scale*B) - log1p(B);
            if ~isreal(log_ratio) || any(~isfinite(log_ratio))
                break
            end
            D_ext = -(aerosol_scale-1)*lidar_ratio*beta ./ differential_cross_section;
            delta_fr = smooth_derivative(log_ratio,C_fr(ids),dz) ./ (2*differential_cross_section) + D_ext;
            usable = isfinite(delta_fr);
            updated = baseline;
            updated(usable) = baseline(usable) + delta_fr(usable);

            if sum(abs(updated-ozone)) <= 0.001*sum(abs(updated))
                beta_a(ids,i) = beta;
                D_fr(ids,i) = delta_fr;
                D_nr(ids,i) = smooth_derivative(log_ratio,C_nr(ids),dz) ./ (2*differential_cross_section) + D_ext;
                break
            end
            ozone = updated;
        end
    end
end

function [beta,converged] = invert_backscatter(step,beta_m,beta,lidar_ratio,dz)
    converged = false;
    for iteration = 1:20
        previous = beta;
        for k = numel(step):-1:1
            tau = lidar_ratio*dz*(beta(k+1) + previous(k));
            beta(k) = exp(step(k)-tau)*(beta_m(k+1)+beta(k+1)) - beta_m(k);
        end
        if any(~isfinite(beta))
            return
        end
        if sum(abs(beta-previous)) <= 0.01*sum(abs(beta)) + eps*sum(beta_m)
            converged = true;
            return
        end
    end
end
