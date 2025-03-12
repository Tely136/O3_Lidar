
function [N_O3_SG,ratio_P]=retrieve_o3ND(pon,poff,frame_len1,frame_len2,frame_len3,h1_hkm,h2_hkm,hkm)
sigmaOn=203.4*10^(-20);% O3 cross section at 287.2 (cm^2/molecule)
sigmaOff=45.51*10^(-20);% O3 cross section at 299.1 (cm^2/molecule)
h1 =find(hkm>h1_hkm,1,'first');
h2 =find(hkm>h2_hkm,1,'first');
d_sigma=(sigmaOn-sigmaOff)*(10^-2)^2; % delta cross section (m^2/molecule)
[len_h,len_t]=size(pon);
ratio_P=pon./poff;
ratio_P(ratio_P<0)=nan;
% ratio_P = sgolayfilt(sgolayfilt(ratio_P,1,sm_len),1,sm_len);
Ln_ratio_P=log(ratio_P);

% Using covolution with the SG filter to calculate the first order
% derivative 

[b1,g1] = sgolay(2,frame_len1);
[b2,g2] = sgolay(2,frame_len2);
[b3,g3] = sgolay(2,frame_len3);
dR=3.75;
N_O3_SG=nan(size(ratio_P));
diff0=nan(size(ratio_P,1),1);

for i=1:len_t
  diff1= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g1(:,2), 'same');
  diff2= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g2(:,2), 'same');
  diff3= conv(Ln_ratio_P(:,i), factorial(1)/(dR)^1 * g3(:,2), 'same');
  diff0(1:h1-1)=diff1(1:h1-1); diff0(h1)=(diff1(h1)+diff2(h1))/2;
  diff0(h1+1:h2-1)=diff2(h1+1:h2-1); diff0(h2)=(diff1(h2)+diff2(h2))/2;
  diff0(h2+1:end)=diff3(h2+1:end);
  N_O3_SG(:,i)=1/(2*d_sigma).*diff0;
end

