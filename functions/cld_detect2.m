function [cldNum,cldBaseZ, cldTopZ,cldBaseZ_ind,cldTopZ_ind,cldBase_qc_flag,cldTop_qc_flag,cld_mask,pz2,d_Pz2]=cld_detect2(prof,start_hkm,end_hkm,hkm, Th_cld1,Th_cld2, framelen,grad_bg)
    [len_h, len_t] = size(prof);
    cld_mask = false(size(prof));
    pz2=nan(size(prof));
    cldNum=zeros(1,len_t);
    cldBase_qc_flag = zeros(2,len_t);
    cldTop_qc_flag = zeros(2,len_t);
    cldBaseZ_ind=nan(2,len_t);
    cldTopZ_ind=nan(2,len_t);
    cldBaseZ=nan(2,len_t);
    cldTopZ=nan(2,len_t);
    
    % Calculate the pz^2
    for i=1:len_t
        pz2(:,i)=prof(:,i).*(hkm).^2;
    end
    log_pz2 = log(pz2);
    
    % derivative by conv with differential filter
    [b,g] = sgolay(2,framelen);
    d_Pz2=nan(size(pz2));
    d_ln_pz2=nan(size(pz2));
    for i=1:len_t
        d_Pz2(:,i) = conv(pz2(:,i), -1*g(:,2), 'same');
        d_ln_pz2(:,i) = conv(log_pz2(:,i), -1*g(:,2), 'same');
    end
    
    % search the cloud gradient between start_bin and end_bin
    start_bin = find(hkm > start_hkm,1,'first');
    end_bin = find(hkm > end_hkm,1,'first');
    mid_bin = ceil((start_bin+end_bin)/2);

    % smooth the profiles for different height
    d_ln_pz2_sm(start_bin:mid_bin,:)= movmean(d_ln_pz2(start_bin:mid_bin,:),40,1);
    d_ln_pz2_sm(mid_bin:end_bin,:)= movmean(d_ln_pz2(mid_bin:end_bin,:),65,1);
    cld_search_hkm = hkm(start_bin:end_bin);
    len_h = length(cld_search_hkm);
    cld_search_region=d_ln_pz2_sm(start_bin:end_bin,:);
    
    for t=1:len_t % for each profile, search from the bottom up
        hr=1;
        cld_rise_id = [];
        cld_fall_id = [];
        while hr < len_h
            if length(cld_rise_id) >1 % already found 2 layers of cloud
                break % break out from while loop, search next profile
            end
            if cld_search_region (hr,t)>Th_cld1 % if the gradient of rising edge > Th
                for hf = hr : len_h
                    if cld_search_region(hf,t)< -Th_cld1 % if find falling edge
                        cld_rise_id = [cld_rise_id;hr]; % save the index of rising and falling edge
                        cld_fall_id = [cld_fall_id;hf];
                        num_cld = length(cld_rise_id);
                        cldBase_qc_flag(num_cld,t)=1;
                        cldTop_qc_flag(num_cld,t)=1;
                        if (sum(cld_search_region(hr:hf,t)>Th_cld2))
                            cldBase_qc_flag(num_cld,t)=2; 
                            cldTop_qc_flag(num_cld,t)=2;
                        end 
                        break; % --> search next rising edge
                    end
                end
                hf_temp = hf;
                hr = hf_temp;
            end
            hr= hr+1;
        end
        % calculate the cloud base height and top height of the two layer clouds
        num_cld = length(cld_rise_id);
        cldNum(t)=num_cld;
        switch num_cld
            case 0 % no cloud, next profile
                continue;
            case 1 % one cloud layer
                % The cloud base height is where the gradient is 0 before the
                % rising edge
                temp_range = cld_search_region(1:cld_rise_id(1),t);
                if ~isnan(find(temp_range<0,1,'last'))
                    cldBaseZ_ind(1,t) = find(temp_range<0,1,'last')+start_bin-1;                
                else
                    cldBaseZ_ind(1,t) = start_bin;
                    cldBase_qc_flag(1,t)=-1;
                    warning("Cloud Base: below the minimum searching height")
                end
                cldBaseZ(1,t) = hkm(cldBaseZ_ind(1,t));
                % cloud top is gradient larger than the gradient baseline, after the falling edge
                temp_range= cld_search_region(cld_fall_id(1):end,t); % not consider the value before the falling edge
                if ~isnan(find(temp_range>grad_bg,1,'first'))
                    cldTopZ_ind(1,t) = find(temp_range>grad_bg,1,'first')+cld_fall_id(1)-1+start_bin-1;
                else
                    I = prctile(temp_range(1:min(200,length(temp_range))),75);
                    if ~isempty (find(temp_range>I,1,'first'))
                        cldTopZ_ind(1,t) = find(temp_range>I,1,'first')+cld_fall_id(1)-1+start_bin-1;
                        cldTop_qc_flag(1,t)=-1;
                        warning("Cloud Top: using the percentile")
                    else
                        cldTopZ_ind(1,t) = end_bin;
                        cldTop_qc_flag(1,t)=-1;
                    end
                    
                end
                cldTopZ(1,t) = hkm(cldTopZ_ind(1,t));
                % mask the signal between the cloud base and top
                cld_mask(cldBaseZ_ind(1,t):cldTopZ_ind(1,t),t) = true ;
            case 2% two cloud layer
                % The first cloud base height is where the gradient is 0 before the
                % first rising edge
                temp_range = cld_search_region(1:cld_rise_id(1),t);
                if ~isnan(find(temp_range<0,1,'last'))
                    cldBaseZ_ind(1,t) = find(temp_range<0,1,'last')+start_bin-1;
                else
                    cldBaseZ_ind(1,t) = start_bin;
                    cldBase_qc_flag(1,t)=-1;
                    warning("Cloud Base: below the minium searching height")
                end
                cldBaseZ(1,t) = hkm(cldBaseZ_ind(1,t));
                % The first cld Top is at first gradient which larger than
                % gradient baseline between the first falling edge and
                % second rising edge
                temp_range = cld_search_region(cld_fall_id(1):cld_rise_id(2),t);
                % grad_bg= mean(cld_search_region(max(1,cld_rise_id(1)-10):cld_rise_id(1),t));% gradient baseline
                cldTopZ_ind(1,t) = find(temp_range>grad_bg,1,'first')+cld_fall_id(1)-1+start_bin-1;
                cldTopZ(1,t) = hkm(cldTopZ_ind(1,t));
                % mask the signal between the cloud base and top
                cld_mask(cldBaseZ_ind(1,t):cldTopZ_ind(1,t),t) = true ;
                
                % The second cldBaseZ_ind is last gradient which > 0 between the first falling edge and
                % second rising edge
                cldBaseZ_ind(2,t) = find(temp_range>0,1,'last')+cld_fall_id(1)-1+start_bin-1;
                cldBaseZ(2,t)=hkm(cldBaseZ_ind(2,t));
                
                % The Second cloud top is gradient larger than the gradient baseline after the second falling edge
                % grad_bg = mean(cld_search_region(max(1,cld_rise_id(1)-10):cld_rise_id(1),t));% gradient baseline
                temp_range= cld_search_region(cld_fall_id(2):end,t); % not consider the value before the 2nd falling edge
                if ~isnan(find(temp_range>grad_bg,1,'first'))
                    cldTopZ_ind(2,t) = find(temp_range>grad_bg,1,'first')+cld_fall_id(2)-1+start_bin-1;
                else
                    I = prctile(temp_range(1:min(200,length(temp_range))),75);
                    l = find(temp_range>I,1,'first')+cld_fall_id(2)-1+start_bin-1;
                    cldTopZ_ind(2,t)=l;
                    warning("Cloud Top: using the percentile")
                    cldTop_qc_flag(2,t)=-1;
                end
                cldTopZ(2,t)=hkm(cldTopZ_ind(2,t));

                % mask the signal between the cloud base and top
                cld_mask(cldBaseZ_ind(2,t):cldTopZ_ind(2,t),t) = true ;
                
                if (cldBaseZ_ind(2,t) - cldTopZ_ind(1,t) <40) % if two layer clouds are very close to each other
                    cldTopZ_ind(1,t) = cldTopZ_ind(2,t);
                    cldTopZ(1,t) = cldTopZ(2,t);
                    cldBaseZ_ind(2,t) = nan;
                    cldBaseZ(2,t) = nan;
                    cldTopZ_ind(2,t) = nan;
                    cldTopZ(2,t) = nan;
                    cldBase_qc_flag(2,t)=0;
                    cldTop_qc_flag(2,t)=0;
                    cld_mask(cldBaseZ_ind(1,t):cldTopZ_ind(1,t),t) = true ;
                    cldNum(t)=1;
                end   
        end
    end
end
