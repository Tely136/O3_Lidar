% Check the signal above the cloud top
function [prof_qc_flag,prof_sig_level]=prof_qc_cldtop(prof,hkm,cldNum,cldBaseZ, cldTopZ,hkm_th1,hkm_th2,hkm_th3)
% hkm_th1 = 6;
% hkm_th2 =8;
% hkm_th3 =9;
[~,len_t] = size(prof);
no_cld_id = (cldNum ==0);
h_id =(hkm > hkm_th1 & hkm < hkm_th2);
Th_prof1 = mean(mean(prof(h_id,no_cld_id),'omitnan'),'omitnan');
h_id =(hkm > hkm_th2 & hkm < hkm_th3);
Th_prof2 = mean(mean(prof(h_id,no_cld_id),'omitnan'),'omitnan');
% check cld profile
prof_qc_flag = nan(2,len_t);
prof_sig_level= nan(2,len_t);
for t=1:len_t
    if cldNum(t) == 0
        continue;
    end 
    if cldNum(t) == 1 % one cloud layer
        if (cldTopZ(1,t) > 8)
            % if high cloud (>8km), no need to
            % check the profile
            prof_qc_flag(1,t) = 2;
            continue;
        end 
        % check the signal after the cloud top
        sig_level = mean(prof((hkm > cldTopZ(1,t)&hkm<8) ,t),'omitnan');
        prof_sig_level(1,t)=sig_level;
        if (sig_level>Th_prof1)
            % good signal after cldTopZ, no need to revisit
            prof_qc_flag(1,t) = 0;
            continue;
        end 
        if (sig_level>Th_prof2 & sig_level<Th_prof1)
            % fair signal after cldTopZ, need to revisit
            prof_qc_flag(1,t) = 1;
            continue;
        end
        if (sig_level<Th_prof2)% bad signal, no need to revisit
            prof_qc_flag(1,t) = 2;
            continue;
        end
    end 
    if cldNum(t) == 2 % two layer cloud
        
       sig_level1 = mean(prof((hkm > cldTopZ(1,t)&hkm<cldBaseZ(2,t)) ,t),'omitnan');
       prof_sig_level(1,t)=sig_level1;
       if (sig_level1<Th_prof2) 
           % if the signal between first and second cloud is bad
           prof_qc_flag(:,t) = 2; 
           continue;
       end
       if (sig_level1 <Th_prof1 & sig_level1 >Th_prof2)
           % if the signal between first and second cloud is fair
           prof_qc_flag(1,t) = 1; % revisit the signal between first and second cloud
           prof_qc_flag(2,t) = 2; % no need to visit the signal after second cloud
           continue;
       end
       if (sig_level1 >Th_prof1)
           
           % if the signal between first and second cloud is good 
           prof_qc_flag(1,t) = 0;
           if (cldTopZ(2,t) > 8)
               prof_qc_flag(2,t) = 2;
               continue;
           end 
           % check the signal after the 2nd cloud top
           sig_level2 = mean(prof((hkm > cldTopZ(2,t)&hkm<8) ,t),'omitnan');
           prof_sig_level(2,t)=sig_level2;
           if (sig_level2>Th_prof1)
               % good signal after cldTopZ, no need to revisit
               prof_qc_flag(2,t) = 0;
               continue;
           end
           if (sig_level2>Th_prof2 & sig_level<Th_prof1)
               % fair signal after cldTopZ, need to revisit
               prof_qc_flag(2,t) = 1;
               continue;
           end
           if (sig_level2<Th_prof2)% bad signal, no need to revisit
               prof_qc_flag(2,t) = 2;
               continue;
           end
           
       end
       
    end
    
end


