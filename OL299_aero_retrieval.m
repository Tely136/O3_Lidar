function [absc_299_retri,ND_O3,N_O3_bsc,D_aext]=OL299_aero_retrieval(N_O3,NDAir_m3_prof,prof_299,hkm)

%% ozone absorption cross section constant
sigmaOn=203.4*10^(-20);% O3 cross section at 287.2 (cm^2/molecule)
sigmaOff=45.51*10^(-20);% O3 cross section at 299.1 (cm^2/molecule)
d_sigma=(sigmaOn-sigmaOff)*(10^-2)^2; % delta cross section (m^2/molecule)
[len_h,len_t]=size(prof_299);
%% Calculate the molecular extinction at 299nm
S2= 8*pi/3;
lamda_1=299; %% unit: nm
WAVE_1=lamda_1/(1e+3);  %% wavelength unit: micro, um
PI=3.1415926; 
ns= (7.247249e+18)* 1013.5/288.15; %% /* unit: molecule/cm3 */
    % 7.247249e+18=Ns*Ys/Ps (S: surface, Ns=2.546899e+19 /cm3, P3=1013.25mbar, Ts=288.15K % 
    
ms1=1.0+(1e-8)*(6432.8 + 2949810/(146- 1.0/WAVE_1/WAVE_1) + 25540/(41-1/WAVE_1/WAVE_1));
sigma1=(6+3*0.035)/(6-7*0.035)*8*PI*PI*PI*(ms1*ms1-1)*(ms1*ms1-1)/3.0/WAVE_1/WAVE_1/WAVE_1/WAVE_1/ns/ns*(1e+16);
               %%% Rayleigh scattering cross section, unit: cm2 
               %%% Corrected by depolarization ratio delta=0.035
am299=sigma1*(NDAir_m3_prof)*(0.1);   %% unit: km^-1   %% Z/km, ext/km, ncn3, Pmb
am287=(287/299)^(-4)* am299;
delta_molex=am287-am299;
D_molex=delta_molex./d_sigma;

% molecular backscatter at 299nm
bm299=am299/S2;


absc_299 = nan(size(prof_299));
aext_299 = nan(size(prof_299));
r299 = nan(size(prof_299));
%% Set up the lidar ratio, aerosol backscatter retrieval range, cloud height 

 
% range bin sizes of the different derivative window length
h1_hkm=2; % 0-2km 1:533
h2_hkm=5;% 2-5km 534:1333
bscframeLen1=31;% ~100 m
bscframeLen2=53; % ~200m
bscframeLen3=81;% 303.75m
bscdh1=find(hkm>h1_hkm,1,'first');
bscdh2=find(hkm>h2_hkm,1,'first');
[b,g1] = sgolay(2,bscframeLen1);
[b,g2] = sgolay(2,bscframeLen2);
[b,g3] = sgolay(2,bscframeLen3);
dz=15;


N_O3_bsc=nan(len_h,len_t);
diff=nan(len_h,1);

% Initial ozone number density

ND_O3 = (N_O3 - D_molex);
D_aext = nan(size(ND_O3));

Sa = 30;% assume S ratio to be 30sr for 299
Z1 = 7; Z2=8;
Z1_c=15; Z2_c=20;
Rc =1;
ae=1.5;
num_iter=5;
rel_err1=zeros(num_iter,len_t);
    
for i =1:len_t % for each profile at certain time
    for k=1:num_iter  % for each iteration epoch
        % ozone absoption extinction
        a_o3=ND_O3(:,i)*(sigmaOff*(1e-4))*1e3;% /km
        % using fernald method to calculate the
        [aero_para]=fernald_func_aercld_o3(hkm,prof_299(:,i),bm299,Sa,Rc,Z1,Z2,Z1_c,Z2_c,a_o3);
        absc_299(:,i)=aero_para(:,3);
        aext_299(:,i)=aero_para(:,2);
        
        %% Calculate the aerosol correction term
        
        absc_287=(299/287)^ae*absc_299(:,i);
        bm287=(299/287)^4*bm299;
        totbsc_299off=absc_299(:,i)+bm299;
        r299(:,i) = totbsc_299off./bm299;
        totbsc_287on=absc_287+bm287;
        ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
        ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
        ln_bsc=log(ratio_totbsc_onoff);
        
        %% Using SG filter to calculate the derivative of log backscatter on-off ratio
        
        diff1= conv(ln_bsc, factorial(1)/(dz)^1 * g1(:,2), 'same');
        diff2= conv(ln_bsc, factorial(1)/(dz)^1 * g2(:,2), 'same');
        diff3= conv(ln_bsc, factorial(1)/(dz)^1 * g3(:,2), 'same');
        diff(1:bscdh1)=diff1(1:bscdh1);
        diff(bscdh1+1:bscdh2)=diff2(bscdh1+1:bscdh2);
        diff(bscdh2+1:end)=diff3(bscdh2+1:end);
        N_O3_bsc(:,i)=1/(2*d_sigma).*diff;
        
        delta_aext=1e-3*((299.1/287.2)^ae-1)*aext_299(:,i);% aerosol extinction in (/m)
        D_aext(:,i)=delta_aext./d_sigma;
        ND_O3_corr = (ND_O3(:,i) - N_O3_bsc(:,i) - D_aext(:,i));
        rel_err1(k,i)= mean(abs((ND_O3(:,i) - ND_O3_corr)./ND_O3(:,i)),'omitnan');
        ND_O3(:,i)=ND_O3_corr;
        if rel_err1(k,i) < 0.01
            break;          
        end           
    end
end 
absc_299_retri.hkm=hkm;
absc_299_retri.bsc=absc_299;
absc_299_retri.ext=aext_299;
absc_299_retri.r=r299;