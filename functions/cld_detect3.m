function [cldNum,cldBaseZ, cldTopZ,cldBaseZ_ind,cldTopZ_ind,cldBase_qc_flag,cldTop_qc_flag,cld_mask,pz2,d_Pz2] = cld_detect3(prof, start_hkm, end_hkm, hkm, Th_cld1, Th_cld2, framelen, grad_bg)
    [len_h,len_t] = size(prof);
    cld_mask = false(size(prof));
    pz2=nan(size(prof));
    cldNum=zeros(1,len_t);
    cldBase_qc_flag = zeros(3,len_t);
    cldTop_qc_flag = zeros(3,len_t);
    cldBaseZ_ind=nan(3,len_t);
    cldTopZ_ind=nan(3,len_t);
    cldBaseZ=nan(3,len_t);
    cldTopZ=nan(3,len_t);
    % Th_cld1=0.0045;
    % Th_cld2=0.007;
    % grad_bg =-0.002;
    Th_cldsig_lnpz2 = 2;
    Th_noncldsig_lnpz2 = 1;
    
    % determine the cloud search region 
    start_bin = find(hkm >=start_hkm,1,'first');
    end_bin = find(hkm<=end_hkm,1,'last');
    mid_bin = ceil((start_bin+end_bin)/2);
    hkm_cldsrc = hkm(start_bin:end_bin);
    len_h_cldsrc = length(hkm_cldsrc);
    
    % 1. Calculate the slope and quality of the signal, as well as the standard 
    % deviation of the background noise level using a range-uncorrected averaged 
    % vertical profile. 
    
    % std_bg = std(prof(end-dlen:end,:),0,1); % Calculate the background noise std
    
    
    %2. Calculate the ln(Pz2)
    for i=1:len_t
        pz2(:,i)=prof(:,i).*(hkm).^2;
    end
    pz2(pz2<0)=nan; 
    ln_pz2=log(pz2);
    ln_pz2_sm=ln_pz2;
    ln_pz2_cldsrc = ln_pz2_sm(start_bin:end_bin,:);
    
    %derivative by conv with differential filter
    [b,g] = sgolay(2,framelen);
    d_Pz2=nan(size(pz2));
    d_ln_pz2=nan(size(pz2));
    d2_ln_pz2=nan(size(pz2));
    for i=1:len_t
        %d_Pz2(:,i) = conv( pz2(:,i), -1*g(:,2), 'same');
        d_ln_pz2(:,i) = conv(ln_pz2_sm(:,i), -1*g(:,2), 'same');
        d2_ln_pz2(:,i) = conv(ln_pz2_sm(:,i), factorial(2)/(-1)^2 * g(:,3), 'same');
    end
    d_lnpz2_cldsrc = d_ln_pz2(start_bin:end_bin,:);
    d2_lnpz2_cldsrc = d2_ln_pz2(start_bin:end_bin,:);
    
    
    
    %
    for t=1:len_t % for each profile, search from the bottom up
        % find all the peaks in the log pz2 profile and record their peak
        % positions, width and prominence
        % base_lnpz2 = whs(ln_pz2_cldsrc(:,t),1e3);
        % temp = ln_pz2_cldsrc(:,t) - base_lnpz2;
        % 
        % [pks,locs] = findpeaks(temp,'MinPeakHeight',0.02);
        % [pks2,locs2] = findpeaks(-temp,'MinPeakHeight',0.01);
        
        cld_id = 0;
        h = 2;
        % figure
        while h <len_h_cldsrc-1
    
            if cldNum(t) == 3
                break % break out from while loop, search next profile
            end
            if ln_pz2_cldsrc(h,t)<Th_noncldsig_lnpz2 % signal too low           
                % low_sig_base = h+start_bin-1;
                % cld_mask(low_sig_base:end,t) = true ;
                % warning("id = %.0d, signal level too low to search for cloud from height: %.2f upward,move to next profile",t,hkm_cldsrc(h));              
                break;
            end 
           
            if d_lnpz2_cldsrc(h,t) > Th_cld1 && (d_lnpz2_cldsrc(h-1,t) >0) && (d_lnpz2_cldsrc(h+1,t)>0) % find a positive peak
    
                cld_id = cld_id +1;
    
                % look for the cloud base: the last zero crossing before h
                zb = find(d_lnpz2_cldsrc(1:h,t)<=0,1,'last');
                cldTop_qc_flag(cld_id,t) = 1;
                
                % subplot(1,3,1);
                % plot(ln_pz2_cldsrc(:,t),hkm_cldsrc);hold on;
                % plot(ln_pz2_cldsrc(zb,t),hkm_cldsrc(zb),'o');
                % xlabel('ln(Pz^2)'); ylabel('altitude (km)'); grid on;
                % title(sprintf('ln(Pz^2), id = %.0d',t));
                % subplot(1,3,2);
                % plot(d_lnpz2_cldsrc(:,t),hkm_cldsrc);hold on;
                % plot(d_lnpz2_cldsrc(zb,t),hkm_cldsrc(zb),'o');
                % xlabel('dln(Pz^2)/dz'); ylabel('altitude (km)'); grid on;
                % title(sprintf('dln(Pz^2)/dz, id = %.0d',t));
    
                if isempty(zb) % no cloud base find
                    warning("id = %.0d, cld_id =%.0d, find cloud rising edge, but can not find cloud base ---> set base to be start search bin",t,cld_id);
                    zb = 1;
                    cldTop_qc_flag(cld_id,t) = -1;
                end
                % look for the cloud peak positition: the first zero crossing above h
                zp = h-1+ find(d_lnpz2_cldsrc(h:end,t)<=0,1,'first');
                % subplot(1,3,1);
                % hold on;
                % plot(ln_pz2_cldsrc(zp,t),hkm_cldsrc(zp),'d');
                % subplot(1,3,2);hold on;
                % plot(d_lnpz2_cldsrc(zp,t),hkm_cldsrc(zp),'d');
                % 
                if isempty(zp) % no cloud peak find
                    warning("id = %.0d, cld_id =%.0d, find cloud rising edge, but can not find cloud peak ---> set the top to be the last search bin",t,cld_id);
                    zp = len_h_cldsrc;
                    zt = len_h_cldsrc;
                    cldTop_qc_flag(cld_id,t) = -1;
                    cldNum(t) = cldNum(t) +1 ;
                    cldBaseZ_ind(cld_id,t)=zb+start_bin-1;
                    cldTopZ_ind(cld_id,t)=zt+start_bin-1;
                    cldBaseZ(cld_id,t)=hkm_cldsrc(zb);
                    cldTopZ(cld_id,t)=hkm_cldsrc(zt);
                    cld_mask(cldBaseZ_ind(cld_id,t):cldTopZ_ind(cld_id,t),t) = true ;
                    break;% break out from while loop, search next profile
                end
                % look for cloud top: search for the last negative peak<-Th
                % within 600 m
                cld_tp_range = zp:min([zp+150,len_h_cldsrc]);
                zt_temp=zp-1+find(d_lnpz2_cldsrc(cld_tp_range,t)<-Th_cld1,1,"last");
                % subplot(1,3,1);
                % hold on;
                % plot(ln_pz2_cldsrc(zt_temp,t),hkm_cldsrc(zt_temp),'*');
                % subplot(1,3,2); hold on
                % plot(d_lnpz2_cldsrc(zt_temp,t),hkm_cldsrc(zt_temp),'*');
    
                if ~isempty(zt_temp) 
                    % cloud top: the first local minmum of second derivative
                    % above the negative peak
                    for ht = zt_temp:len_h_cldsrc-1
                        zt = ht;
                        cldTop_qc_flag(cld_id,t) =-1;
                        if(d2_lnpz2_cldsrc(ht,t) < d2_lnpz2_cldsrc(ht+1,t)) && (d2_lnpz2_cldsrc(ht,t) < d2_lnpz2_cldsrc(ht-1,t))
                            cldTop_qc_flag(cld_id,t) = 1;
                            break;
                        end
                    end
                    if ht == len_h_cldsrc-1
                        zt = len_h_cldsrc;
                        warning("id = %.0d, cld_id =%.0d, no local minmum of second derivative was found ---> set the top to be the last search bin",t,cld_id)
                    end
                else % can't find the negative peak over threshold, 
                    zt = len_h_cldsrc;
                    cldTop_qc_flag(cld_id,t) = -1;
    
                end
                % subplot(1,3,1);
                % hold on;
                % plot(ln_pz2_cldsrc(zt,t),hkm_cldsrc(zt),'^');
                % subplot(1,3,2); hold on
                % plot(d_lnpz2_cldsrc(zt,t),hkm_cldsrc(zt),'^');
                % subplot(1,3,3); hold on
                % plot(d2_lnpz2_cldsrc(:,t),hkm_cldsrc);hold on;
                % plot(d2_lnpz2_cldsrc(zt,t),hkm_cldsrc(zt),'^');
                % xlabel('d^2ln(Pz^2)/dz^2'); ylabel('altitude (km)'); grid on;
                % title(sprintf('d^2ln(Pz^2)/dz^2, id = %.0d',t));
    
                h = zt;
                cldNum(t) = cldNum(t) +1 ;
                cldBaseZ_ind(cld_id,t)=zb+start_bin-1;
                cldTopZ_ind(cld_id,t)=zt+start_bin-1;
                cldBaseZ(cld_id,t)=hkm_cldsrc(zb);
                cldTopZ(cld_id,t)=hkm_cldsrc(zt);
                cld_mask(cldBaseZ_ind(cld_id,t):cldTopZ_ind(cld_id,t),t) = true ;
            end
            h = h+1;         
        end
    
    end
    % check the signal level above the cloud
    for t=1:len_t
        if cldNum(t) > 0  % if there is cloud    
            idtop=cldTopZ_ind(1,t)-start_bin+1;
            for h=idtop:len_h_cldsrc
                if ln_pz2_cldsrc(h,t)<Th_cldsig_lnpz2 % signal too low
                    low_sig_base = h+start_bin-1;
                    cld_mask(low_sig_base:end,t) = true ;
                    break;
                end
            end
        else % if there is no cloud 
            for h=1:len_h_cldsrc
                if ln_pz2_cldsrc(h,t)<Th_noncldsig_lnpz2 % signal too low
                    cld_mask(h:end,t) = true ;
                    break;
                end
            end
    
        end
    end
    
    p = input("Do you want to show the individual cloud results? Y/N [N]:","s");
    if isempty(p)
        p = 'N';
    end
    if strcmp(p,'Y')||strcmp(p,'y')
        for ii=1:len_t
            set(0,'DefaultTextFontSize',[14]);   
            set(0,'DefaultAxesFontSize',[14]);
            if cldNum(ii) > 0
                idcld=~isnan(cldBaseZ_ind(:,ii));
                figure
                subplot(1,3,1);
                plot(ln_pz2(:,ii),hkm,'LineWidth',1.2);hold on;
                plot(ln_pz2(cldBaseZ_ind(idcld,ii),ii),cldBaseZ(idcld,ii),'o','LineWidth',1.2);
                plot(ln_pz2(cldTopZ_ind(idcld,ii),ii),cldTopZ(idcld,ii),'^','LineWidth',1.2);
                xlabel('ln(Pz^2)','FontSize',14); ylabel('altitude (km)','FontSize',14); grid on;
                title(sprintf('ln(Pz^2), id = %.0d',ii)); ylim([hkm_cldsrc(1),hkm_cldsrc(end)]);
                subplot(1,3,2);
                plot(d_ln_pz2(:,ii),hkm,'LineWidth',1.2);hold on;
                plot(d_ln_pz2(cldBaseZ_ind(idcld,ii),ii),cldBaseZ(idcld,ii),'o','LineWidth',1.2);
                plot(d_ln_pz2(cldTopZ_ind(idcld,ii),ii),cldTopZ(idcld,ii),'^','LineWidth',1.2);
                xlabel('dln(Pz^2)/dz','FontSize',14); ylabel('altitude (km)','FontSize',14); grid on;ylim([hkm_cldsrc(1),hkm_cldsrc(end)]);
                title(sprintf('dln(Pz^2)/dz, id = %.0d',ii));
                subplot(1,3,3);
                plot(d2_ln_pz2(:,ii),hkm,'LineWidth',1.2);hold on;
                plot(d2_ln_pz2(cldBaseZ_ind(idcld,ii),ii),cldBaseZ(idcld,ii),'o','LineWidth',1.2);
                plot(d2_ln_pz2(cldTopZ_ind(idcld,ii),ii),cldTopZ(idcld,ii),'^','LineWidth',1.2);
                xlabel('d^2ln(Pz^2)/dz^2','FontSize',14); ylabel('altitude (km)','FontSize',14); grid on;ylim([hkm_cldsrc(1),hkm_cldsrc(end)]);
                title(sprintf('d^2ln(Pz^2)/dz^2, id = %.0d',ii));
            end
        end
    else
        return
    end 
    
    end
