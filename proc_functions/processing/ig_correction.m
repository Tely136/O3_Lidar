function D = ig_correction(Nig,dsigma_ig,dsigma_o3)

    D = Nig .* dsigma_ig ./ dsigma_o3; 

end