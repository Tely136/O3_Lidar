%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fernald Function: Inversion for retrieving aerosol extinction, backscatter and scattering ratio 
%% Refernence: Fernald, F. G., 1984. Analysis of atmospheric lidar observations: some comments. Appl. Opt. 23, 652-653.
%% Written by Yonghua Wu at CCNY 
%% Modified by Dingdong Li at CCNY 2022 May, including the ozone absorption 
%% and retrieve the signal of 299nm for Ozone DIAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [aero_para]=fernald_func_aercld_o3(zkm,lidar_p,bm,S1,Rc,Z1,Z2,Z1_c,Z2_c,ao3) 
% Input variables (array): 
%    zkm: Range or altitude in km
%    lidar_p: lidar signal profile subtracting from background noise
%    bm: moleclar backscatter coefficient in /km/sr (wavelength conistency with the lidar signals)
%    S1: aerosol extinction-backscatter-ratio or lidar-ratio (e.g. 50-sr at 532-nm, 40-sr at 1064-nm)
%    Rc: reference value of aerosol scattering-ratio (e.g. Rc =1 at 299-nm, Rc=1.01 at 532-nm, 1.06-1.08 at 1064-nm)
%    Z1, Z2: range window in searching for free-aerosol layer in km
%    ao3: ozone extinction in /km
%    Z1_c, Z2_c: thin-cloud base and top in km. Under clear sky: Z1_c and Z2_c can be set as higher altitude, e.g. 15-km, 20-km etc.

% Output or return [aero_para]:
%    aero_para]=zkm, ext(/km), Ba(/km/sr), R(scattering-ratio), Pz2 (Range-corrected signal)

S2=8*3.1415916; % Molecular extinction-to-backscatter ratio
S1_cloud=18.5;  % Cloud extinction-to-backscatter ratio

% Calculate the p(z)*z2/mol_bs/Tm2
x_dim=size(lidar_p,1);
% x_dim=size(xx,1);
% xx: the range corrected signal (lidar_p.*zkm^2)* attenuation due to ozone
T_o3=zeros(x_dim,1)+1;
for i=2:x_dim
    a_temp_o3=trapz(zkm(1:i),ao3(1:i));
    T_o3(i)=exp(-2*a_temp_o3);
end  
% xx=lidar_p.*zkm.*zkm./T_o3;
xx=lidar_p.*zkm.*zkm;
% lidar_p=xx./(zkm.^2);

% Random noise Pb for the reference: last 200-points (267*3.75m=1000 m )
pts=200;
sd=std(lidar_p(x_dim-pts:x_dim),'omitnan');
ave=mean(lidar_p(x_dim-pts:x_dim),'omitnan');
SNR=lidar_p/sd;

% Tm2=zeros(x_dim,1)+1;
% for i=2:x_dim
%     a=trapz(zkm(1:i),bm(1:i)*S2);
%     Tm2(i)=exp(-2*a);
% end   
% b=[zkm,Tm2];
% ww=xx./bm./Tm2;   %% from Sassno's paper in Appl Opt., 1985
ww=xx./bm;

% Search for the reference-altitude Zc where there is the minimum "ww" (i.e. clear air layer or free-aerosol layer)
kk=(zkm>=Z1 & zkm<=Z2);
xx2=ww(kk);
yy2=zkm(kk);
pp=lidar_p(kk);
aa=min(xx2);

kk=(xx2==aa);
Zc=yy2(kk);
Pc=pp(kk);
SNR_Zc=SNR(kk);

kk=(zkm>=(Zc-0.5) & zkm<=(Zc+0.5)) ;
hh=ww(kk);
%fprintf('** S1=%4.2f  Zc=%4.4f  Pc-mv=%5.1e SNR=%4.2f  minX=%4.1e maxX=%4.1e Xc=%4.1e Rc=%4.3f *\n',S1,Zc,Pc,SNR_Zc,min(hh),max(hh),min(xx2),Rc);

r=[];   b=[]; ext=[];
kk=find(zkm==Zc);
r(kk)=Rc;
b(kk)=(r(kk)-1.0)*bm(kk);  %% aerosol backscatter coffecient at the reference altitude
ext(kk)=b(kk)*S1;  %% aerosol extinction coffecient at the reference altitude (free-aerosol or less aerosol)

% Integration for retrieving aerosol
% Backward integration order-by-order for aerosol backscatter coefficient
for i=kk:-1:2
      if( zkm(i-1)>=Z1_c & zkm(i-1)<=Z2_c)
          S1_aer=S1_cloud;
      else 
          S1_aer=S1;
      end
      aa=(S1_aer-S2)*(bm(i)+bm(i-1))*(zkm(i)-zkm(i-1));
      mm=xx(i-1)*exp(aa);
      b(i-1)=mm/(xx(i)/(b(i)+bm(i))+S1_aer*(xx(i)+xx(i-1)*exp(aa) )*(zkm(i)-zkm(i-1)))-bm(i-1);
      r(i-1)=1.0 +b(i-1)/bm(i-1);
      ext(i-1)=b(i-1)*S1_aer;
end
%% Forward integration
for i=kk:x_dim-1
      if(zkm(i+1)>Z1_c & zkm(i+1)<Z2_c) 
          S1_aer=S1_cloud;
      else
          S1_aer=S1;
      end    
      aa=(S1_aer-S2)*(bm(i)+bm(i+1))*(zkm(i)-zkm(i-1));
      mm=xx(i+1)*exp(-aa);
      b(i+1)=mm/(xx(i)/(b(i)+bm(i))-S1_aer*(xx(i)+xx(i+1)*exp(-aa) )*(zkm(i)-zkm(i-1)))-bm(i+1);
      r(i+1)=1.0+b(i+1)/bm(i+1);
      ext(i+1)=b(i+1)*S1_aer;      
 end  
 aero_para=[zkm, ext', b',r',xx]; %% z-km, ext/km, backscatter/km/sr, scattering ratio, pz2
  
return;
