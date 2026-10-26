function u = uNO3_MRG_NC(uNO3_up,uNO3_lo,w)
% Eq.95: independent LOW/HIGH RANGE DET components. w is LOW-range weight.
% merge_profiles returns HIGH-range weight: pass 1-w_high here.
% Correlation is BETWEEN RANGES, not altitude bins. Merge components first.
    u=hypot(w.*uNO3_lo,(1-w).*uNO3_up);
end
