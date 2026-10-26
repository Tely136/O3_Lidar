function [total,merged] = merge_pc_uncertainty(low,high,wLow)
% Eqs.95/97 and Sec.4.11: merge components BEFORE quadrature total.
% wLow is LOW-range weight, unlike merge_profiles' returned HIGH weight.
% Caller must use the same effective weights as the retrieval (including
% its omitmissing fallback). Inputs have fields DET,BKG,SAT only.
    sz=size(low.DET);
    wLow=wLow+zeros(sz);
    if ~isequal(size(wLow),sz) || any(wLow(:)<0|wLow(:)>1)
        error('O3:MergeWeights','Weights must match profiles and lie in [0,1].');
    end
    fields={'DET','BKG','SAT'};
    for j=1:numel(fields)
        name=fields{j}; lo=low.(name); hi=high.(name);
        if ~isequal(size(lo),sz) || ~isequal(size(hi),sz)
            error('O3:MergeSize','Component profiles must have identical sizes.');
        end
        if strcmp(name,'DET'), v=uNO3_MRG_NC(hi,lo,wLow);
        else, v=uNO3_MRG_FC(hi,lo,wLow); end
        % Avoid 0*NaN poisoning when a range has zero retrieval weight.
        v(wLow==1)=lo(wLow==1); v(wLow==0)=hi(wLow==0);
        merged.(name)=v;
    end
    total=sqrt(merged.DET.^2+merged.BKG.^2+merged.SAT.^2);
end
