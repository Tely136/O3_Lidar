function D = rayleigh_correction(Na,dsigma_m,dsigma_o3)
    D = Na .* dsigma_m ./ dsigma_o3; 
end