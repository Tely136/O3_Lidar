function [N_absc,D_aext]=aero_corr(wl_nm,ae,aero_bsa,aero_ext,am,dsigma_hkm,frame_len,dz)

%% Calculate the aerosol correction term
% prompt = 'Do you want to use the Aeronet Angstrom Exponent? (if no aeronet AE available or choosing no, ae=1.5 will be used)? Y/N [N]:';
% str = input(prompt,'s');
% if isempty(str)
%     str = 'N';
% end
% if str=='N'|str=='n'|isnan(ae_1020_340)|ae_1020_340<1
%     ae=1.5;
% else
%     [file,path] = uigetfile('/Users/Tinker/Documents/MATLAB/ozonelidar/AOD/*.txt',...
%                         'Select an AOD File');
%     aodfilename = [path,file];                
%     aeronet_aod=loadAODdata(aodfilename);
%     ae=ae_1020_340;
% end
[len_h,len_t]=size(aero_bsa);
S1= 8*pi/3;
absc_299=(wl_nm/299)^ae*aero_bsa;
absc_287=(wl_nm/287)^ae*aero_bsa;
mbsc_299=(wl_nm/299)^4/S1*am;
mbsc_287=(wl_nm/287)^4/S1*am;
totbsc_299off=absc_299+mbsc_299;
totbsc_287on=absc_287+mbsc_287;
ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
ln_bsc=log(ratio_totbsc_onoff);

% Using second order polyfit to calculate the derivative of log backscatter on-off ratio

N_absc=nan(size(ln_bsc));


hf0 = round((frame_len(1)-1)/2);
hf1 = round((frame_len(end)-1)/2);
for i=1:len_t
    for j=1+hf0:len_h-hf1
        hf = round((frame_len(j)-1)/2);
        [~,g0] = sgolay(2,frame_len(j));
        diff0= ln_bsc(j-hf:j+hf,i)'*g0(:,2);
        N_absc(j,i)= -(1/2)*diff0/dsigma_hkm(j)/dz;
    end
end

delta_aext=1e-3*((wl_nm/287.2)^ae-(wl_nm/299.1)^ae)*aero_ext;% aerosol extinction in (/m)
D_aext=delta_aext./double(dsigma_hkm);