function u = uNO3_MRG_FC(uNO3_up,uNO3_lo,w)
% Eq.97: fully correlated LOW/HIGH RANGE SAT/BKG. w is LOW-range weight;
% merge_profiles returns HIGH-range weight, so pass 1-w_high.
% Atmospheric terms stay unchanged (Eqs.93-94). Merge components then totals.
    u=abs(w.*uNO3_lo+(1-w).*uNO3_up);
end
