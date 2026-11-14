function u = uNO3_MRG_FC(uNO3_up,uNO3_lo,w)
    % Eq. 97: w is the low-range weight.

    u = abs(w .* uNO3_lo + (1 - w) .* uNO3_up);
end
