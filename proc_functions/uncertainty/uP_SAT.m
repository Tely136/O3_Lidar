function u = uP_SAT(P,dz,L,uTau)
    % Eqs. 34-35: P in accumulated counts; L is TOTAL shots; uTau in seconds.

    u = physconst("LightSpeed") ./ (2 .* dz .* L) .* P.^2 .* uTau;
    u(~isfinite(P)) = NaN;
end
