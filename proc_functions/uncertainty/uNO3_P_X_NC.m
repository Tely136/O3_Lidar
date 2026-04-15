% Uncertinty in Ozone number density due to X propapgated throguh photon counts, with no covariance terms
% From Thierry's Part 2 paper
%   Equations 32,33,36,37,43,44 in paper
function u = uNO3_P_X_NC(Pon,Poff,u_Pon_X,u_Poff_X,dsimga,dz,cp)
    arguments
        Pon         % signal counts in on channel
        Poff        % signal counts in off channel
        u_Pon_X     % uncertainty in photon counts due to source X in on channel
        u_Poff_X    % uncertainty in photon counts due to source X in off channel
        dsimga      % differential absorption cross section
        dz          % bin width
        cp          % height dependent filter coefficients
    end

    coeff = 1 ./ (abs(dsimga) * dz);

    % Altitude dimension should be 1st
    nk = size(Pon,1);
    n_prof = size(Pon,2);

    u = NaN(size(Pon));

    for prof = 1:n_prof
        for i = 1:nk
            c_temp = cp{i};
            M = size(c_temp,1);
            N = (M-1)/2;

            if i-N >= 1 && i+N <= nk
                idx = i-N:i+N;
                A = ((u_Pon_X(idx)./Pon(idx).^2 + (u_Poff_X(idx)./Poff(idx).^2)));
                u(i,prof) = sqrt(c_temp(:,2)' * A);
            end
        end
    end

    u = coeff .* u;
end