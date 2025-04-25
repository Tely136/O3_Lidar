%% OL299_aero_iter_retr
% retrieve aerosol backscatter coefficient at 299nm and the correction
% terms
%
%% Input:
% no3_init_m3: Initial O3 number density with only molecular correction (molecule/m3)
% am:  molecular extinction (/m)
% prof_fr_nr: Ozone-lidar signal intensity at 299-nm (mV or photon-count)
% hkm: Ozone-lidar signal altitude in km
% Sm: Molecular extinction-to-backscatter-ratio(sr): 8*pi/3
% Sa: Aerosol extinction-to-backscatter ratio (sr)
% zsurf:     the surface altitude of the profile (km)
% zref: reference altitude in km where the Rc is assigend
% Rc: Reference value for the total backscatter to molecuar backsatter ratio (Rc=1 for pure molecular air or aerosol-free air)
%% Return:
% absc_299: aerosol backscatter coefficient profile at 299nm
% ND_O3: ozone number density corrected by aerosol backscatter, molecular
% and aerosol extinction.
% N_O3_bsc: ozone number density correction term constributed by the aerosol gradiant change
% D_aext:ozone number density correction term contributed by aerosol
% extinction

%%
function [absc_299,N_O3_bsc,D_aext,no3,rel_erro3]=OL299_aero_iter_retr(no3_init_m3,am,prof_fr_nr,hkm,zsurf,zref,Rc,Sa)



Sm= 8*pi/3;
bm = am/Sm;
% zsurf=0.3;
% zref=7;
% Rc=1.005;  % reference value for the total-backscatter to molecular backscatter ratio
% Sa= 60;

no3_init_m3(hkm<zsurf)=nan;% the O3 number density below 300m set to be nan
N_iter=2;
absc_299 = nan(length(prof_fr_nr),1);
aext_299 = nan(size(absc_299));
N_O3_bsc = nan(size(absc_299));
D_aext= nan(size(absc_299));
rel_erro3=[];
no3_cm3= no3_int_m3./1e6;
bm_km=bm299*1e3;
[aero_para]=aer_iterative_func(hkm_nr,prof_299_fr_nr,bm_km,no3_cm3,Sa,Sm,Rc,zref);
for k = 1 : N_iter
    if k==1
        no3 = no3_init_m3;
    end
    [absc_299,aext_299,eps]=OL_299nm_aer_iter(prof_fr_nr, hkm, no3, am, bm, zsurf, zref, Rc, Sa);
    
    [N_O3_bsc,D_aext] = OL_299nm_o3_iter(hkm,absc_299,bm,aext_299);
    
    no3_new = no3_init_m3-fillmissing(N_O3_bsc,'constant',0)-fillmissing(D_aext,'constant',0);
    
    eps2= sum(abs(no3-no3_new),'omitnan')/sum(no3,'omitnan');
    
    rel_erro3=[rel_erro3,eps2];   
    
    no3 = no3_new;
    
    if eps2<0.01   
        break
    end 
end