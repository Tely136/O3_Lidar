function dsdz = smooth_derivative(s,C,dz)
    arguments
        s
        C
        dz
    end


    dsdz = NaN(size(s));
    
    for profile_idx = 1:size(s,2)
        for k = 1:size(s,1)
            filt_diff = C{k};
            N = (length(filt_diff) - 1)/2;
    
            if k-N >= 1 && k+N <= size(s,1)
                dsdz(k,profile_idx) = filt_diff(:,2)' * s(k-N:k+N,profile_idx);

            elseif k-N < 1
            elseif k+N > size(s,1)

            end
        end
    end
    
    dsdz = dsdz./dz;
end