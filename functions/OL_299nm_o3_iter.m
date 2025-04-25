%% ozone iteration 
function [N_O3_bsc,D_aext] = OL_299nm_o3_iter(hkm,absc_299,bm299,aext_299)
%% constant 
sigmaOn=203.4*10^(-20)*(10^-2)^2;% O3 cross section at 287.2 (m^2/molecule)
sigmaOff=45.51*10^(-20)*(10^-2)^2;% O3 cross section at 299.1 (m^2/molecule)
d_sigma=(sigmaOn-sigmaOff); % delta cross section (m^2/molecule)

hm=hkm*1e3;
% range bin sizes of the different derivative window length
h1_hkm=2; % 0-2km 1:533
h2_hkm=5;% 2-5km 534:1333
bscframeLen1=21;% ~78.5 m
bscframeLen2=53; % ~200m
bscframeLen3=81;% 303.75m
bscdh1=find(hkm>h1_hkm,1,'first');
bscdh2=find(hkm>h2_hkm,1,'first');
[b,g1] = sgolay(2,bscframeLen1);
[b,g2] = sgolay(2,bscframeLen2);
[b,g3] = sgolay(2,bscframeLen3);
dz=hm(2)-hm(1);

ae = 1.5;
absc_287=(299/287)^ae*absc_299;
bm287=(299/287)^4*bm299;
totbsc_299off=absc_299+bm299;
r299= totbsc_299off./bm299;
totbsc_287on=absc_287+bm287;
ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
ln_bsc=log(ratio_totbsc_onoff);

% Using SG filter to calculate the derivative of log backscatter on-off ratio
diff=nan(size(ln_bsc));
diff1= conv(ln_bsc, factorial(1)/(dz)^1 * g1(:,2), 'same');
diff2= conv(ln_bsc, factorial(1)/(dz)^1 * g2(:,2), 'same');
diff3= conv(ln_bsc, factorial(1)/(dz)^1 * g3(:,2), 'same');
diff(1:bscdh1)=diff1(1:bscdh1);
diff(bscdh1+1:bscdh2)=diff2(bscdh1+1:bscdh2);
diff(bscdh2+1:end)=diff3(bscdh2+1:end);
N_O3_bsc=1/(2*d_sigma).*diff;
delta_aext=((299/287)^ae-1)*aext_299;% aerosol extinction in (/m)
D_aext=delta_aext./d_sigma;

end 