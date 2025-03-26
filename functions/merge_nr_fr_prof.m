%
%% Description: 
% Merge the far range and near range signal using linear fit
%% Input
% prof_fr: double, len_hfr x len_t, the far range signal profiles
% prof_nr: double, len_hnr x len_t, the far range signal profiles
% hkm_fr: double, len_hfr x 1, the vertical altitude for the far range 
% hkm_nr: double, len_hnr x 1, the vertical altitude for the near range 
% z1: double, 1x1, the altitude that fitting begins (km)
% z2: double, 1x1, the altitude that fitting ends (km)
%% Output
% prof_merge:  double, len_h x len_t, the far range signal profiles
function prof_merge = merge_nr_fr_prof(prof_nr,prof_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr,z1,z2)
z1_fr_ind = find(hkm_fr>z1,1,'first');
z1_nr_ind = find(hkm_nr>z1,1,'first');
z2_fr_ind = find(hkm_fr>z2,1,'first');
z2_nr_ind = find(hkm_nr>z2,1,'first');
len_t= size(prof_fr,2);

%
prof_merge=prof_nr;

% the range below start_merge_hkm 
ind_cat_low_nr= hkm_nr<=start_merge_hkm;% the range that <= start_merge_hkm

% the range above end_merge_hkm using prof_fr
ind_cat_high_nr= hkm_nr>=end_merge_hkm;% the range that >= end_merge_hkm
ind_cat_high_fr= hkm_fr>=end_merge_hkm;
prof_merge(ind_cat_high_nr,:)= prof_fr(ind_cat_high_fr,:);
% Concanation of two channel
ind_hh_fr=hkm_fr<end_merge_hkm & hkm_fr>start_merge_hkm;
ind_hh_nr=hkm_nr<end_merge_hkm & hkm_nr>start_merge_hkm;
% merge range
h_concat=hkm_fr(ind_hh_fr);
w=(h_concat-h_concat(1))./(h_concat(end)-h_concat(1));
for i=1:len_t
   p= polyfit(prof_nr(z1_nr_ind:z2_nr_ind,i),prof_fr(z1_fr_ind:z2_fr_ind,i),1);
   % below the start_merge, using linear fit of prof_nr
   prof_merge(ind_cat_low_nr,i) = polyval(p,prof_nr(ind_cat_low_nr,i));
 % in between start_merge_hkm and end_merge_hkm using the weighted average
 % of prof_fr and linear fit using prof_nr
   prof_merge(ind_hh_nr,i) = w.*(prof_fr(ind_hh_fr,i))+(1-w).*polyval(p,prof_nr(ind_hh_nr,i));
end 