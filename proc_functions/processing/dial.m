function N = dial(Mon,Moff,C,dz,dsigma)

% data that comes in is already background corrected, gate corrected, an-pc
% merged if applicable, and time-averaged

% focus on returning signal terms and correction terms

    r = log(Moff./Mon);

    [n_row,n_col] = size(r); % height by time
    drdz = NaN(size(r));

    for i = 1:n_col
        for k = 1:n_row
            filt_diff = C{k};
            N = (length(filt_diff) - 1)/2;

            if k-N >= 1 && k+N <= n_row
                drdz(k,i) = filt_diff(:,2)' * r(k-N:k+N,i);

            % elseif k-N < 1
            % elseif k+N > size(s,1)

            end
        end
    end

    S = drdz ./ dz;

    N = S ./ (2*dsigma); % dsigma should be same shape as S, and show temperature dependence
end