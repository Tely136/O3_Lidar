%% Iterative method to calculate the ozone number density
% 
% This program is adopted from the paper Kuang, S., Burris, J. F., Newchurch, 
% M. J., Johnson, S., & Long, S. (2010). Differential absorption lidar to measure 
% subhourly variation of tropospheric ozone profiles. IEEE Transactions on 
% Geoscience and Remote Sensing, 49(1), 557-571.
% Input:
% P         lidar signal used for aerosol correction
% h         altitude (km)
% nO3         ozone number density without aerosol correction (/m^3)
% sigma_off ozone absorption cross section of the off line (m^2/molecule)
% d_sigma   delta cross section (m^2/molecule)
% am        molecular extinction (/m)
% bm        molecular backscatter (/m/sr)
% zref      the reference height, notice the the program only calculate
%           from the zsurf to zref (km)
% zsurf     the surface height of the profile (km)
% baref     the aerosol backscatter coefficient at reference height (km)
% S         the aerosol lidar ratio
% Output:
% ba        aerosol backscatter coefficient profile (/m/sr)
% n_b       aerosol backscatter correction term for ozone number density (/m^3)
% n_e       extinction correction term for ozone number density (/m)
%% 

function [ba_out, aa_out,eps_aer ]=OL_299nm_aer_iter(P, h, nO3, am, bm, zsurf, zref, Rc, Sa)
sigmaOff=45.51*10^(-20)*(10^-2)^2;% O3 cross section at 299.1 (m^2/molecule)

P_=nan(size(P));
P_(1:end-1)=P(2:end);% P_=P(r+dr)
Ln_P_ratio=log(P./P_); % ln(P(r)/P(r+dr))
% calculate aerosol correction term from the zref downwards
[~,refInd] = min(abs(h-zref));
[~,sufInd] = min(abs(h-zsurf));
r = 1e3*h;
baref=(Rc-1)*bm(refInd);
dr = r(2)-r(1);
n_iter = 3;
ba_out = nan;
ba= nan(length(h),n_iter);
ba(refInd,:)=baref;
eps_aer = nan(1,n_iter-1);
%% aerosol iteration
for k=1:n_iter % iterate n_iter times to get converged profiles

    for i=refInd-1:-1:sufInd % calculated from the reference height downwards
        % Calculate the αm(r+Δr/2)= (αm(r+Δr) + αm(r))/2
        am_mid=(am(i)+am(i+1))/2;
        if k==1
            % For the first iteration, αa(r+Δr/2) ≈ αa(r+Δr) = Sβa(r+Δr)
            aa_mid = Sa*(ba(i+1,k));
        else
            % otherwise using αa(r+Δr/2) = S(βa(r+Δr) + βa'(r))/2, βa' is from previous iteration  
            aa_mid = Sa*(ba(i+1,k)+ba(i,k-1))/2;
        end
        % equ(10) in Kuang et al. 2011 IEEE
        o3tau = (nO3(i)+nO3(i+1))*sigmaOff*dr;
        amtau= 2*(am_mid+aa_mid)*dr;
        term3= ((r(i)/r(i+1))^2)*((bm(i+1)+bm(i))/2+ba(i+1,k));      
        ba(i,k)=exp(Ln_P_ratio(i)-o3tau-amtau)*term3-bm(i);        
    end
    if k>1 % if have more than one calculated profile, compare their relative diff
        eps_aer(k)= sum(abs(ba(sufInd:refInd,k)-ba(sufInd:refInd,k-1)),'omitnan')/sum(ba(sufInd:refInd,k-1),'omitnan');
%         if eps_aer(k) < 0.02
            % if the relative diff of βa from two iterations are less than 1% 
%             formatStr= "The profile converges within %d iteration.\n The relative error: %.3f\n";
%             fprintf(formatStr,k,eps_aer(k));
            ba_out=ba(:,k);
            aa_out= Sa*ba(:,k);
%             break
%         else
%             formatStr= "The profile doesn't converge within %d iteration.\n The relative error: %.3f\n";
%             fprintf(formatStr,k,eps_aer(k));
%         end       
    end
    %
end

