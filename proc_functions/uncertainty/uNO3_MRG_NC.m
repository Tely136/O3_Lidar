function u = uNO3_MRG_NC(uNO3_up,uNO3_lo,w)
    % Eq. 95: w is the low-range weight.

    u = hypot(w .* uNO3_lo,(1 - w) .* uNO3_up);
end
