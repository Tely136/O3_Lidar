% function for Gluing the LICEL A/D and Photon-Counting(PC) data
% By Yonghua Wu at CCNY 
% Last updated on Apr.19, Apr.25, 2022
% Updated by Dingdong Li Jul. 27th, 2022

function [glue, regp, regr] = adpc_glue_func_cldscr(ad,pc,z1_reg,z2_reg,zkm,cldn,cbz,ctz)
    %% Input parameters:
    % ad: A/D signals (matrix) at the channel-1
    % pc: Photon-counting (PC) data at the channel-1 
    % z1_reg: start regression height
    % z2_reg: end regression height
    % z1_mrg: start merge height
    % z2_mrg: end merge height
    % zkm: Altitude in km
    % cldn: Number of clouds
    % cbz: Cloud base height in km
    % ctz: Cloud top height in km
    % zkm0: Initial altitude for the valid signals (skip the gating-off altitude)
    %% Output parameters:
    % glue: the glued A/D and PC signal profiles
    %         (low-range:Converted A/D to PC sig.; middle or glue-range:average PC and AD; high-range: PC data)
    % regp: regression slope (a) and intercept(b) (y=ax+b)  (for diagose only)
    % regr: regression correlation coefficient (for diagose only)
    
    [len_h,len_t] = size(ad);
    regp = nan(len_t,2);
    regr = nan(len_t,1);
    glue = nan(size(ad));
    overlap_flag = false(len_t);
    inz_reg = (zkm > z1_reg & zkm < z2_reg);
    zkm_reg = zkm(inz_reg);
    indz1_mrg = find(inz_reg,1,'first');
    indz2_mrg = find(inz_reg,1,'last');
    for i = 1:len_t
        x = ad(inz_reg,i);
        y = pc(inz_reg,i);
        regp(i,:) = polyfit(x,y,1);
        r = corr(x,y);
        regr(i) = r;
        ad_reg = polyval(regp(i,:),ad(:,i));
    
        glue(1:indz1_mrg,i) = ad_reg(1:indz1_mrg);
        glue(indz2_mrg:end,i) = pc(indz2_mrg:end,i);
        glue(indz1_mrg:indz2_mrg,i) = (ad_reg(indz1_mrg:indz2_mrg)+pc(indz1_mrg:indz2_mrg,i))./2;
     
        if (cbz(1,i) >z1_reg && cbz(1,i)<z2_reg) || (ctz(1,i)>z1_reg && ctz(1,i)<z2_reg) || (cbz(2,i)>z1_reg && cbz(2,i)<z2_reg) || (ctz(2,i)>z1_reg && ctz(2,i)<z2_reg )   % Wu revised
            overlap_flag(i) = true;
        end
    end
    
    for i = 1:len_t
        ind_p = 0;
        if (overlap_flag(i) && (regr(i)<0.95)) || (regr(i)<0.85) % if there is cloud appear in the range or R<0.85 
            % use the coefficient closest to with corr higher than 0.85
            for j = i:len_t % look for the value after i
                if regr(j) >= 0.85 && ~overlap_flag(j) % if R>0.85 and no cloud
                    ind_p = j;
                    break;
                end
            end

            for j = i:-1:1 % look for the value before i
                if regr(j)>0.85 && ~overlap_flag(j) % if R>0.85 and no cloud
                    ind_p = j;
                    break;
                end
            end
            
            if ind_p == 0 % Cant find any profile with corr > 0.85 and no cloud
                % chose the regression coeff with the highest r 
                [~,ind_p] = max(regr);
                glue(:,i) = polyval(regp(ind_p,:),ad(:,i));
                warning('id = %.0d, R is too low, R =%.2f, glue signal only use the converted AD signal, convert coefficient use id=%.0d, R(%.0d) = %.2f',i,regr(i),ind_p,ind_p,regr(ind_p));
            else
                glue(:,i) = polyval(regp(ind_p,:),ad(:,i));
                if overlap_flag(i)
                    warning('id = %.0d, cloud appear in the merge range, R =%.2f, glue signal only use the converted AD signal, convert coefficient use id=%.0d, R(%.0d) = %.2f',i,regr(i),ind_p,ind_p,regr(ind_p));
                else
                    warning('id = %.0d, R is too low, R =%.2f, glue signal only use the converted AD signal, convert coefficient use id=%.0d, R(%.0d) = %.2f',i,regr(i),ind_p,ind_p,regr(ind_p));
                end            
            end
        end
    end

return;