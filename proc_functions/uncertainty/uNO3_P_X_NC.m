function u = uNO3_P_X_NC(Pon,Poff,u_Pon_X,u_Poff_X,dsigma,dz,cp,verticalCorrelation)
% Independent ON/OFF; altitude x profile arrays, input u values are SDs.
% dsigma SINGLE-PATH (sigma_on-sigma_off), matching dial.m.
% 'none' default: Eq.32, independent altitude bins, squared relative SDs.
% 'full': LITERAL Eqs.36/43, signed sum of pointwise quadrature terms.
% Printed correlated formulas can return NEGATIVE values; see README.
    if nargin<8, verticalCorrelation='none'; end
    mode=validatestring(verticalCorrelation,{'none','full'});
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
            rel=hypot(ua./a,ub./b);
            if strcmp(mode,'none'), v=sqrt(sum(c.^2.*rel.^2));
            else, v=sum(c.*rel); end
            if isfinite(dsigma(k,prof)) && dsigma(k,prof)~=0
                u(k,prof)=v./(2.*abs(dsigma(k,prof)).*dz);
            end
        end
    end
end
