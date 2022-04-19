function N_O3_concat = merge_nr_fr_o3ND(N_O3_nr,N_O3_fr,start_merge_hkm,end_merge_hkm,hkm)
% Concanation of two channel
ind_hh=hkm<end_merge_hkm & hkm>start_merge_hkm;
N_O3_concat=N_O3_fr;
N_O3_concat(~(hkm>start_merge_hkm),:)=N_O3_nr(~(hkm>start_merge_hkm),:);
% merge range
h_concat=hkm(ind_hh);
% weight of high channel
w=(h_concat-h_concat(1))./(h_concat(end)-h_concat(1));
N_O3_concat(ind_hh,:)=w.*(N_O3_fr(ind_hh,:))+(1-w).*(N_O3_nr(ind_hh,:));
