function u = uNO3_P_X_FC(Pon,Poff,u_Pon_X,u_Poff_X,dsigma,dz,cp)
% LITERAL Eqs.38/45: fully correlated ON/OFF AND altitude bins.
% dsigma SINGLE-PATH: denominator=2*abs(dsigma). Preserve printed inner
% absolute value and signed sum. May return NEGATIVE values; see README.
    [Pon,Poff,uon,uoff,dsigma]=uncertainty_inputs(Pon,Poff,u_Pon_X,u_Poff_X,dsigma,dz,cp);
    u=NaN(size(Pon));
    for prof=1:size(Pon,2)
        for k=1:size(Pon,1)
            f=cp{k}; c=f(:,2); n=(numel(c)-1)/2;
            if n~=floor(n), error('O3:FilterLength','Filters must have odd length.'); end
            if k-n<1 || k+n>size(Pon,1), continue; end
            idx=k-n:k+n;
            a=Pon(idx,prof); b=Poff(idx,prof); ua=uon(idx,prof); ub=uoff(idx,prof);
            if any(~isfinite(a)|~isfinite(b)|a<=0|b<=0|~isfinite(ua)|~isfinite(ub)), continue; end
            if isfinite(dsigma(k,prof)) && dsigma(k,prof)~=0
                u(k,prof)=sum(c.*abs(ua./a-ub./b))./(2.*abs(dsigma(k,prof)).*dz);
            end
        end
    end
end
