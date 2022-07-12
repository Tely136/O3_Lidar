function [cldNum, cldBaseZ, cldCenterZ,cldBaseZ_ind,cldCenterZ_ind,cld_mask,pz2,d_Pz2]=cld_detect(prof,start_hkm,end_hkm,hkm)
[len_h,len_t]=size(prof);
cld_mask = false(size(prof));
pz2=nan(size(prof));
cldNum=zeros(1,len_t);
cldBaseZ_ind=nan(3,len_t);
cldCenterZ_ind=nan(3,len_t);
cldBaseZ=nan(3,len_t);
cldCenterZ=nan(3,len_t);
Th_cld=10e-3;

%Calculate the pz^2
for i=1:len_t
  pz2(:,i)=prof(:,i).*(hkm).^2;
end
%take derivative on successive elements in the columns
d_Pz2=diff(pz2,1,1);
% search the cloud gradient between start_bin and end_bin
start_bin = find(hkm >start_hkm,1,'first');
end_bin = find(hkm >end_hkm,1,'first');
cld_search_region = d_Pz2(start_bin:end_bin,:);
N=27;% cloud gradient thickness threshold: the number of consecutive elements of the over_threshold gradient of a cloud

% find all the index of gradient profile larger than Threashold
id_cld=cld_search_region>Th_cld;
for i=1:len_t
    if nnz(id_cld(:,i))> 0
        % the index of the begining of the consecutive large gradient
        id_cld_base=strfind([0 id_cld(:,i)' 0],[0 1]);
        % the index of the end of the consecutive large gradient
        id_cld_center=strfind([0 id_cld(:,i)' 0],[1 0]) -1 ;
        %check which number of the consecutive large gradient is larger than N
        ii=(id_cld_center-id_cld_base + 1)>=N;
        if nnz(ii)>0 % if there is any gradient sequence that is longer than N
            
            cldNum(i)=min(3,nnz(ii)); % the number of cloud, at most 3;
            cldBaseInd_temp = id_cld_base(ii)+ start_bin -1;
            cldCenterInd_temp = id_cld_center(ii)+ start_bin -1;
            
            num_cld = min(3,nnz(ii)); % the number of cloud, at most 3;
            cldBaseZ_ind(1:num_cld,i)= cldBaseInd_temp(1:num_cld); % lowest three cloud base index
            cldCenterZ_ind(1:num_cld,i)= cldCenterInd_temp(1:num_cld); % lowest three cloud center index
            cldBaseZ(1:num_cld,i)= hkm(cldBaseInd_temp(1:num_cld)); % lowest three cloud base index
            cldCenterZ(1:num_cld,i)=hkm(cldCenterInd_temp(1:num_cld)); % lowest three cloud center index

        end
    end
    
end
%% 
for i=1:len_t
    % check if there is cloud, then mask the signal above the cloud base
    if cldNum(i)>0 % if there is cloud
        cld_mask(cldBaseZ_ind(1,i):end,i) = true ;
    end
end