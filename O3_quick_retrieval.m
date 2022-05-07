function [o3_dial_retr]=O3_quick_retrieval(folderpath,savepath,nbin,bgbins,td,...
                                       nAvg,start_bin_fr,start_bin_nr,....
                                       cld_start_bin,cld_end_bin)
%% Description:
% This function reads in the path of the folder which stores CCNY-O3_DIAL
% data acquisition txt files and display the O3 number density and mixing
% ratio using the air density calculated from standard atmosphere model and
% intermediate signal profiles
% 
%% Syntax:
% [o3_dial_retrieval]=O3_quick_retrieval(folderpath,savepath,nbin,bgbins,td,...
%                                        nAvg,start_bin_fr,start_bin_nr,....
%                                        cld_start_bin,cld_end_bin)
%
%% Input:
% folder_path: string, the folder path which saves the data .txt files 
% save_path: string, the path to save the results
% nbin: int 1x1, number of bins of each profile
% dz: float 1x1 the height in meter of each bin, for example, dz= 3.75m a 40MHz system
% bgbins: int 1x1, number of bins to be averaged at the end of the profile as background
% td: float 1x1 dead time constant, to calculate the true count rate Ct=Cm/(1-td*Cm)
% nAvg: int 1x1 number of bins to perform time average 
% start_bin_fr: int 1x1 ordinal number of bins that turns on the far range gate ;
% start_bin_nr: int 1x1 ordinal number of bins that turns on the near range gate ;
%% Output:
% o3_dial_retr: structure, 9 fields
% o3_dial_retr.new_sigprof : structure, 12 fields, time-averaged, deadtime-corrected, background-subtracted signal profiles and pc,ad merged signal profiles
% o3_dial_retr.hkm_fr: float (nbin-start_bin_fr+1)x1 The altitude of the far-channel sigal
% o3_dial_retr.hkm_nr: float (nbin-start_bin_nr+1)x1 The altitude of the near-channel sigal
% o3_dial_retr.DateTime_avg: datetime, Dim-(1, floor(numFiles/nAvg)), datetime of averaged profiles
% o3_dial_retr.TimeInHour_avg: float, Dim-(1, floor(numFiles/nAvg)), time in hour of averaged profiles
% o3_dial_retr.cld_mask:  float (nbin-start_bin_nr+1)xfloor(numFiles/nAvg) cloud mask
% o3_dial_retr.ND_O3:  float (nbin-start_bin_nr+1)xfloor(numFiles/nAvg) Ozone number density w/t molecular correction
% o3_dial_retr.ND_Air: float (nbin-start_bin_nr+1)xfloor(numFiles/nAvg) Air number using
% standard model
% o3_dial_retr.O3_ppb: float (nbin-start_bin_nr+1)xfloor(numFiles/nAvg) Ozone mixing ratio w/ molecular correction
dzm=3.75;%m
[OLfileName,~]=read_OL_profiles2(folderpath,savepath,nbin,dzm,bgbins,td,nAvg);
load(OLfileName);
%% Remove the signals which collect befor the gate opens
% start_bin_fr=200;
% start_bin_nr=14;
[new_sigprof,hkm_fr,hkm_nr] = gate_prof(sigprof,start_bin_fr,start_bin_nr,hkm);
%% Calculate the merge ad pc signal
[new_sigprof.merge287, ~]=adpc_glue_func(new_sigprof.an287,new_sigprof.pc287,5,10,hkm_fr,1);
[new_sigprof.merge299, ~]=adpc_glue_func(new_sigprof.an299,new_sigprof.pc299,5,10,hkm_fr,1);
[new_sigprof.merge287nr, ~]=adpc_glue_func(new_sigprof.an287nr,new_sigprof.pc287nr,3,7,hkm_nr,1);
[new_sigprof.merge299nr, ~]=adpc_glue_func(new_sigprof.an299nr,new_sigprof.pc299nr,3,7,hkm_nr,1);
%% Cloud screening 
% take derivative of the signal profile from the start bin
% the clouds are the place:1. the signal at cloud base is larger than the
% signal below the cloud 2. negative derivative < Threshold
% smooth the profiles

sgwin_len_cld=31;% must be odd
cld_screen_prof=sgolayfilt(new_sigprof.an299,1,sgwin_len_cld);

% cld_start_bin = 1;% 1km
% cld_end_bin = 10;% 10km
[~, cldBaseZ, ~,~, ~,cld_mask,~,~]=cld_detect(cld_screen_prof,cld_start_bin,cld_end_bin,hkm_fr);
%% cloud screened signal
prof_merge_287=new_sigprof.merge287; prof_merge_287(cld_mask)=nan;
prof_merge_299=new_sigprof.merge299; prof_merge_299(cld_mask)=nan;
prof_an_287_nr=new_sigprof.an287nr; 
prof_an_299_nr=new_sigprof.an299nr; 

%% calculating the ozone number density from the pon and poff
% SG filter derivative filter window length
frame_len1=31;% ~100 m
frame_len2=53; % ~200m
frame_len3=81;% 303.75m
 
% range bin sizes of the different derivative window length
h1_hkm=2; % 0-2km 1:533
h2_hkm=5;% 2-5km 534:1333

[N_O3_fr,~]=retrieve_o3ND(prof_merge_287,prof_merge_299,...
                                   frame_len1,frame_len2,frame_len3,h1_hkm,h2_hkm,hkm_fr);
[N_O3_nr,~]=retrieve_o3ND(prof_an_287_nr,prof_an_299_nr,...
                                   frame_len1,frame_len2,frame_len3,1.5,h2_hkm,hkm_nr);
                               
%% merge the near range and far range ozone number density
start_merge_hkm =0.9;
end_merge_hkm =1;
N_O3_merge = merge_nr_fr_o3ND(N_O3_nr,N_O3_fr,start_merge_hkm,end_merge_hkm,hkm_fr,hkm_nr);
%% calculate air density
% option 1: using International Standard Atmosphere (ISA)

sigmaOn=203.4*10^(-20);% O3 cross section at 287.2 (cm^2/molecule)
sigmaOff=45.51*10^(-20);% O3 cross section at 299.1 (cm^2/molecule)
d_sigma=(sigmaOn-sigmaOff)*(10^-2)^2; % delta cross section (m^2/molecule)

p0 = 1013.25*1e2; % surface pressure 1013.25 hpa = 101325pa
T0 = 288.15;% surface temp (K)
[NDAir_m3_prof,D_molex_prof] = isa_NDair_m3(hkm_nr,p0,T0,d_sigma);
len_t=length(TimeInHour_avg);
NDAir_m3= repmat(NDAir_m3_prof, 1,len_t);
D_molex = repmat(D_molex_prof, 1,len_t);

ND_O3_corr = (N_O3_merge - D_molex);
O3ppb = ND_O3_corr./NDAir_m3*1e9;
%% Plot the ND 
y_limit = [0.27,8.5]; %km
z_limit = [0,120];
fontsize= 13;
title_str = [datestr(DateTime_avg(1),'yyyy/mm/dd'),' CCNY-DIAL ozone mixing ratio (ppb)'];
NDplot(O3ppb,TimeInHour_avg,hkm_nr,title_str,y_limit,z_limit,fontsize);

y_limit = [0.27,8.5]; %km
z_limit = [0,2.5e18];
fontsize= 13;
title_str = [datestr(DateTime_avg(1),'yyyy/mm/dd'),' CCNY-DIAL ozone number density (/m^3)'];
NDplot(ND_O3_corr,TimeInHour_avg,hkm_nr,title_str,y_limit,z_limit,fontsize);
%% Save the files 
o3_dial_retr.new_sigprof =new_sigprof;
o3_dial_retr.hkm_fr=hkm_fr;
o3_dial_retr.hkm_nr=hkm_nr;
o3_dial_retr.DateTime_avg=DateTime_avg;
o3_dial_retr.TimeInHour_avg=TimeInHour_avg;
o3_dial_retr.cld_mask=cld_mask;
o3_dial_retr.cbh=cldBaseZ;
o3_dial_retr.ND_O3=ND_O3_corr;
o3_dial_retr.ND_Air=NDAir_m3;
o3_dial_retr.O3ppb = O3ppb;