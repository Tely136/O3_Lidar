function [Pon,Poff,uon,uoff,dsigma] = uncertainty_inputs(Pon,Poff,uon,uoff,dsigma)
    % Match scalar or altitude-vector inputs to the profiles.

    uon = expand_input(uon,size(Pon));
    uoff = expand_input(uoff,size(Pon));
    dsigma = expand_input(dsigma,size(Pon));
end

function x = expand_input(x,sz)
    if isvector(x) && numel(x) == sz(1)
        x = x(:);
    end
    x = x + zeros(sz);
end
