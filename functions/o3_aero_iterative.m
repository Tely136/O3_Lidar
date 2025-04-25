%% Iterative method to calculate the ozone number density
% 
% This program is adopted from the paper Kuang, S., Burris, J. F., Newchurch, 
% M. J., Johnson, S., & Long, S. (2010). Differential absorption lidar to measure 
% subhourly variation of tropospheric ozone profiles. IEEE Transactions on 
% Geoscience and Remote Sensing, 49(1), 557-571.
% Input:
% P         lidar signal used for aerosol correction
% h         altitude (m)
% nO3         ozone number density without aerosol correction (/m^3)
% sigma_off ozone absorption cross section of the off line (m^2/molecule)
% d_sigma   delta cross section (m^2/molecule)
% am        molecular extinction (/m)
% bm        molecular backscatter (/m/sr)
% zref      the reference height, notice the the program only calculate
%           from the zsurf to zref (m)
% zsurf     the surface height of the profile (m)
% baref     the aerosol backscatter coefficient at reference height (m)
% S         the aerosol lidar ratio
% Output:
% ba        aerosol backscatter coefficient profile (/m/sr)
% n_b       aerosol backscatter correction term for ozone number density (/m^3)
% n_e       extinction correction term for ozone number density (/m)
%% 

function [ba_out, iter, eps]=o3_aero_iterative(P, h, nO3, am, bm, zsurf, zref, baref, Sa)
%% constant 
sigmaOn=203.4*10^(-20)*(10^-2)^2;% O3 cross section at 287.2 (m^2/molecule)
sigmaOff=45.51*10^(-20)*(10^-2)^2;% O3 cross section at 299.1 (m^2/molecule)
d_sigma=(sigmaOn-sigmaOff); % delta cross section (m^2/molecule)

P_=nan(size(P));
P_(1:end-1)=P(2:end);% P_=P(r+dr)
Ln_P_ratio=log(P./P_); % ln(P(r)/P(r+dr))
% calculate aerosol correction term from the zref downwards
[~,refInd] = min(abs(h-zref));
[~,sufInd] = min(abs(h-zsurf));
r = h(sufInd:refInd);
dr = r(2)-r(1);
n_iter = 4;
ba_out = nan;
ba= nan(length(h),n_iter);
ba(refInd,1)=baref;
ba_temp = ba(:,1);
eps = nan(1,n_iter-1);
%% aerosol iteration
for k=1:n_iter % outer loop of iteration

    for i=refInd:-1:sufInd+1 % calculated from the reference height downwards
        % Calculate the αm(r+Δr/2)= (αm(r+Δr) + αm(r))/2
        am_mid=(am(i-1)+am(i))/2;
        if k==1
            % For the first iteration, αa(r+Δr/2) ≈ αa(r+Δr) = Sβa(r+Δr)
            aa_mid = Sa*(ba_temp(i));
        else
            % otherwise using αa(r+Δr/2) = S(βa(r+Δr) + βa(r))/2
            aa_mid = Sa*(ba_temp(i)+ba_temp(i-1))/2;
        end
        % equ(10) in Kuang et al. 2011 IEEE
        ba_temp(i-1)=exp(Ln_P_ratio(i-1)-2*(nO3(i-1))*sigmaOff*dr...
            -2*(am_mid+aa_mid)*dr)...
            *((h(i-1)/h(i))^2)*(bm(i)+ba_temp(i))...
            -bm(i-1);        
    end
    ba(:,k)=ba_temp;
    if k>1 % if have more than one calculated profile, compare their relative diff
        eps(k-1)= sum(abs(ba(sufInd:refInd,k-1)-ba(sufInd:refInd,k)))/sum(ba(sufInd:refInd,k-1));
        if eps(k-1) < 0.01
            % if the relative diff of βa from two iterations are less than 1% 
            formatStr= "The profile converges within %d iteration.\n The relative error: %.3f\n";
            fprintf(formatStr,k,eps(k));
            ba_out=ba(:,k);
            aa_out= Sa*ba(:,k);
            break
        else
            formatStr= "The profile doesn't converge within %d iteration.\n The relative error: %.3f\n";
            fprintf(formatStr,k,eps(k));
        end       
    end
    %
end
iter = k;
%% ozone iteration 

ae = 1.5;
absc_287=(299/287)^ae*ba_out;
bm287=(299/287)^4*bm;
totbsc_299off=ba_out+bm;
r299 = totbsc_299off./bm;
totbsc_287on=absc_287+bm287;
ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
ln_bsc=log(ratio_totbsc_onoff);

%% Using SG filter to calculate the derivative of log backscatter on-off ratio

diff1= conv(ln_bsc, factorial(1)/(dr)^1 * g1(:,2), 'same');
diff2= conv(ln_bsc, factorial(1)/(dr)^1 * g2(:,2), 'same');
diff3= conv(ln_bsc, factorial(1)/(dr)^1 * g3(:,2), 'same');
diff(1:bscdh1)=diff1(1:bscdh1);
diff(bscdh1+1:bscdh2)=diff2(bscdh1+1:bscdh2);
diff(bscdh2+1:end)=diff3(bscdh2+1:end);
N_O3_bsc=1/(2*d_sigma).*diff;

delta_aext=1e-3*((299.1/287.2)^ae-1)*aa_out;% aerosol extinction in (/m)
D_aext=delta_aext./d_sigma;

ND_O3_corr = (ND_O3(:,i) - N_O3_bsc(:,i) - D_aext(:,i));
rel_err1(k,i)= mean(abs((ND_O3(:,i) - ND_O3_corr)./ND_O3(:,i)),'omitnan');
ND_O3(:,i)=ND_O3_corr;
