function N_O3_concat = merge_nr_fr_o3ND(N_O3_nr,N_O3_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr)
N_O3_concat=N_O3_nr;
% the range above end_merge_hkm using N_O3_fr
ind_cat_high_nr= ~(hkm_nr<end_merge_hkm);% the range that >= end_merge_hkm
ind_cat_high_fr= ~(hkm_fr<end_merge_hkm);
N_O3_concat(ind_cat_high_nr,:)= N_O3_fr(ind_cat_high_fr,:);

% Concanation of two channel
ind_hh_fr=hkm_fr<end_merge_hkm & hkm_fr>start_merge_hkm;
ind_hh_nr=hkm_nr<end_merge_hkm & hkm_nr>start_merge_hkm;
% merge range
h_concat=hkm_fr(ind_hh_fr);
% weight of high channel
w=(h_concat-h_concat(1))./(h_concat(end)-h_concat(1));
N_O3_concat(ind_hh_nr,:)=w.*(N_O3_fr(ind_hh_fr,:))+(1-w).*(N_O3_nr(ind_hh_nr,:));
