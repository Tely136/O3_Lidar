function u = uNO3_MRG_NC(uNO3_up,uNO3_lo,w)
    % Eq. 95: w is LOW-range weight, unlike merge_profiles.

    u = hypot(w .* uNO3_lo,(1 - w) .* uNO3_up);
end
