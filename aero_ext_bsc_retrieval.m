%% function of CHM15k aerosol backscatter retrieval using the radiometer and aeronet Aod
% load the ceilometer attenuated backscatter coefficient
%clear all; close all;
function [chm15kRetrieval,N_O3_bsc_interp,D_aext_interp]=aero_ext_bsc_retrieval(NDAir_m3_prof,d_sigma,hkm,TimeInHour_avg,timediff,aodfile)
[chm15data,pathname]=newgetbschm15k_longrange('choose at least two files','chm15k data');
ind=1;
    TimeArr=chm15data.AveTime{ind};
    DateTimeArr=chm15data.AveDateTime{ind};
    DataAve=chm15data.Avedata{ind};
    datechosed=chm15data.Date{ind};
    pbltopAveTime=chm15data.pbltop{ind};
    pblfp2=chm15data.pblfp2{ind};
    cbhAveTime=chm15data.cbh{ind};
    figure
    Iave=imagesc(TimeArr,chm15data.Height,DataAve,[0,0.1]);
    set(gca,'YDir','normal');%,'ColorScale','log');
    ylim([0,15])
    colormap('jet');
    colorbar;
    xlabel('Time (UTC Hour)','FontSize',13);
    ylabel('Altitude (km)','FontSize',13);
    hold on
    plot(TimeArr(~pblfp2),pbltopAveTime(~pblfp2),'m+','MarkerSize',8);
    plot(TimeArr(pblfp2),pbltopAveTime(pblfp2),'rx','MarkerSize',8);
    plot(TimeArr,cbhAveTime,'go','MarkerSize',8);
    title([datestr(chm15data.Date{ind},'yyyymmdd'),' CCNY Ceilometer RCS (a.u.) at 1064nm'],'FontSize',13);
    legend('PBLH','CldBase','FontSize',11);
    

prompt = 'What is the maximum retrieval range(km)? ';
x = input(prompt)
if isempty(x)
    x = 7.5;
end
m = floor(x*1000/15)

prompt = 'What is the maximum range of the extinction integrated AOD(km)? ';
y = input(prompt)
if isempty(y)
    y = 4;
end
h1 = floor(y*1000/15)

% load the ceilometer attenuated backscatter coefficient
rv0=chm15data.Height(1:m);
file_datestr=datestr(datetime(datechosed,'InputFormat','dd-MMM-yyyy'),'yyyymmdd');
study_date=file_datestr;
%  Load the radiosonde-measured temperature and pressure profiles
 
%% Function: Calculate molecular extinction at 1064-nm
lamda_1=1064; %% unit: nm
WAVE_1=lamda_1/(1e+3);  %% wavelength unit: micro, um
PI=3.1415926; 
ns= (7.247249e+18)* 1013.5/288.15; %% /* unit: molecule/cm3 */
    % 7.247249e+18=Ns*Ys/Ps (S: surface, Ns=2.546899e+19 /cm3, P3=1013.25mbar, Ts=288.15K % 
    
ms1=1.0+(1e-8)*(6432.8 + 2949810/(146- 1.0/WAVE_1/WAVE_1) + 25540/(41-1/WAVE_1/WAVE_1));
sigma1=(6+3*0.035)/(6-7*0.035)*8*PI*PI*PI*(ms1*ms1-1)*(ms1*ms1-1)/3.0/WAVE_1/WAVE_1/WAVE_1/WAVE_1/ns/ns*(1e+16);
               %%% Rayleigh scattering cross section, unit: cm2 
               %%% Corrected by depolarization ratio delta=0.035
molex_1=sigma1*(NDAir_m3_prof)*(0.1);   %% unit: km^-1   %% Z/km, ext/km, ncn3, Pmb

hm=hkm;  %% Z-km
am=molex_1;  %% mol_ext
S2=(8*pi)/3; 

am_interp=interp1(hm,am,rv0,'linear','extrap');% interp the molecular extinction to the height of ceilometer
bm_interp=am_interp/S2;  % mol.bs at 1064-nm, km-1
% molecular two-way transmittance
tm2=zeros(size(am_interp))+1;

for mk=2:length(rv0)
    taum=trapz(rv0(1:mk),am_interp(1:mk));
    tm2(mk)=exp(-2*taum);
end

 %% load AOD data

% aodfile='/Users/Tinker/Documents/MATLAB/Ceilometer/AOD/20210923_20211011_CCNY/20210923_20211011_CCNY.txt';
aod_aeronet=loadAODdata(aodfile);

ind_aod=aod_aeronet.date_utc==datetime(datechosed,'InputFormat','dd-MMM-yyyy');
aod=aod_aeronet.aod_1020(ind_aod);
aod_time=aod_aeronet.TimeInHour_utc(ind_aod);
ae_1020_340=mean(aod_aeronet.ae_1020_340(ind_aod),'omitnan');

% movmean of 30min 
pbltopAveTime(pblfp2)=nan;
pbl_aero=movmean(pbltopAveTime,3,'omitnan');
cbh_aero=cbhAveTime;
%attbs=movmean(chm15data.Avedata{ind}(1:m,:),3,2,'omitnan');
attbs=chm15data.Avedata{ind}(1:m,:);
attbs(1:13,:)=ones(13,1)*attbs(13,:);
time_aeroprof=chm15data.AveTime{ind};


%range movmean <1.5km 100, 1.5-3km 200m, >3km 300m
movnum1=7;
movnum2=13;
movnum3=20;
attbs_sm=nan(m,length(time_aeroprof));

temp=movmean(attbs,movnum1,1);
attbs_sm(1:100,:)=temp(1:100,:);
temp=movmean(attbs,movnum2,1);
attbs_sm(100:200,:)=temp(100:200,:);
temp=movmean(attbs,movnum3,1);
attbs_sm(200:end,:)=temp(200:end,:);
attbs_sm = movmean(attbs_sm,movnum1,1);


% attbs_sm=nan(m,length(time_aeroprof));
aero_zkm=nan(m,length(time_aeroprof));
aero_ext=nan(m,length(time_aeroprof));
aero_bsa=nan(m,length(time_aeroprof));
aero_aod_cal=nan(1,length(time_aeroprof));
aero_r=nan(m,length(time_aeroprof));
aero_lidarratio=nan(1,length(time_aeroprof));
% 

cl_chm15=162;
%% fix the lidar constant with aod
prompt = 'Do you want to use the Aeronet AOD to constrain the lidar ratio? (if no aeronet AOD available or choosing no, S=40sr will be used)? Y/N [N]:';
str = input(prompt,'s');
if isempty(str)
    str = 'N';
end
if str=='N'|isempty(aod_time)
    lidar_r=40;
    lidar_r_array=40*ones(1,length(time_aeroprof));
else
    lidar_r_array=nan(1,length(time_aeroprof));
   [~,ind_start]=min(abs(time_aeroprof-aod_time(1)));
   [~,ind_end]=min(abs(time_aeroprof-aod_time(end)));
   for i=ind_start:ind_end
       if cbh_aero(i)<rv0(h1)|cbh_aero(i-1)<rv0(h1)|cbh_aero(i+1)<rv0(h1)
           continue;
       else
           zc1=cbh_aero(i);
           zc2=20;
       end
        lidar_r=[20:1:70];
        aero_aod=nan(1,length(lidar_r));
        for j=1:length(lidar_r)
            [aero_para,~]=iterative_func2(rv0,attbs_sm(:,i)./cl_chm15,bm_interp,tm2,lidar_r(j),zc1,zc2,18.5);
            aero_aod(j)=trapz(aero_para(:,1),aero_para(:,3));
        end
        [~,ind_i]=min(abs(aod_time-time_aeroprof(i)));
        [~,ind_lidar_r]=min(abs(aero_aod-aod(ind_i)));
        lidar_r_array(i)=lidar_r(ind_lidar_r);
   end
   lidar_r_array(lidar_r_array<22|lidar_r_array>68)=nan;
   lidar_r=mean(lidar_r_array,'omitnan');
   lidar_r_array(isnan(lidar_r_array))=lidar_r;
end
            

% lidar_r_array_temp=lidar_r_array;
for i=1:length(time_aeroprof)% for day time
%     check if there is low cloud, if there is low cloud do not calculate,
%        if cbh_aero(i)<rv0(h1)|cbh_aero(i-1)<rv0(h1)|cbh_aero(i+1)<rv0(h1)
%            continue;
%        else
        zc1=cbh_aero(i);
        zc2=cbh_aero(i)+1;
%         zc1=cbh_aero(i);
%         zc1=15;
%         zc2=20;
        
%        end
    
    
    [aero_para ntimes]=iterative_func2(rv0,attbs_sm(:,i)./cl_chm15,bm_interp,tm2,lidar_r_array(i),zc1,zc2,18.5); 
    aero_para(rv0>(cbh_aero(i)-0.12),:)=nan;
    aero_para(rv0>(cbh_aero(max(i-1,1))-0.12)& rv0<(cbh_aero(max(i-1,1))+0.2),:)=nan;
    aero_para(rv0>(cbh_aero(min(i+1,length(time_aeroprof)))-0.12)& rv0<(cbh_aero(min(i+1,length(time_aeroprof)))+0.2),:)=nan;
    aero_zkm(:,i)=aero_para(:,1);
    aero_ext(:,i)=aero_para(:,3);
    aero_bsa(:,i)=aero_para(:,2);
    aero_r(:,i)=1+aero_para(:,2)./bm_interp;
    aero_aod_cal(i)= trapz(aero_para(:,1),aero_para(:,3));

end

prompt = 'Do you want to display the CHM15k retrieval results? Y/N [Y]:';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end
if str=='Y'
    

    figure
    subplot(2,2,1)
    Iaero=imagesc(time_aeroprof,rv0,aero_bsa,[0,1e-3]);
    hold on
    plot(time_aeroprof,cbh_aero,'go')
    plot(time_aeroprof,pbl_aero,'m+')
    set(Iaero,'AlphaData',~isnan(aero_bsa))
    set(gca,'YDir','normal');
    colormap('jet');
    colorbar;
    legend('CldBase','PBLH')
    xlabel('Time (UTC Hour)');
    ylabel('Altitude (km)');
    title([file_datestr,' CHM15k aerosol backscatter (/km/sr) at 1064nm']);
    
    
    subplot(2,2,2)
    Iaero=imagesc(time_aeroprof,rv0,aero_ext,[0,5e-2]);
    set(Iaero,'AlphaData',~isnan(aero_bsa))
    set(gca,'YDir','normal');
    colormap('jet');
    colorbar;
    xlabel('Time (UTC Hour)');
    ylabel('Altitude (km)');
    title([file_datestr,' CHM15k aerosol extinction (/km) at 1064nm']);
    
    subplot(2,2,3)
    Iaero=imagesc(time_aeroprof,rv0,aero_r,[0,10]);
    set(Iaero,'AlphaData',~isnan(aero_bsa))
    set(gca,'YDir','normal');
    colormap('jet');
    colorbar;
    xlabel('Time (UTC Hour)');
    ylabel('Altitude (km)');
    title([file_datestr,' CHM15k scatter ratio (/km) at 1064nm']);
    
    subplot(2,2,4)
    plot(time_aeroprof,aero_aod_cal,'.');hold on
    plot(aod_time,aod,'r.');
    legend('AOD at 1064 from CHM15k aerosol extinction retrieval','AOD at 1020nm from AERONET');
    xlabel('Time (UTC Hour)');
    ylabel('AOD');
    grid on
    title([file_datestr,'AOD of CHM15k retrieval and AERONENT, S=',num2str(lidar_r,'%.1f'),'\pm',num2str(std(lidar_r_array),'%.1f')])
    
end
%% Calculate the aerosol correction term
prompt = 'Do you want to use the Aeronet Angstrom Exponent? (if no aeronet AE available or choosing no, ae=1.5 will be used)? Y/N [N]:';
str = input(prompt,'s');
if isempty(str)
    str = 'N';
end
if str=='N'|str=='n'|isnan(ae_1020_340)|ae_1020_340<1
    ae=1.5;
else
    ae=ae_1020_340;
end
S1= 8*pi/3;
absc_299=(1064/299)^ae*aero_bsa;
absc_287=(1064/287)^ae*aero_bsa;
mbsc_299=(1064/299)^4/S1*am_interp;
mbsc_287=(1064/287)^4/S1*am_interp;
totbsc_299off=absc_299+mbsc_299;
totbsc_287on=absc_287+mbsc_287;
ratio_totbsc_onoff=totbsc_287on./totbsc_299off;
ratio_totbsc_onoff(ratio_totbsc_onoff<0)=nan;
ln_bsc=log(ratio_totbsc_onoff);

% Using second order polyfit to calculate the derivative of log backscatter on-off ratio
bscframeLen1=5;% 75m 0-2km 1:133
bscdh1=133;
bscframeLen2=13;% 195m  2-5km 134:333
bscdh2=333;
bscframeLen3=21;% 315m  >5km  334:end

% Using SG filter to calculate the derivative of log backscatter on-off ratio

[b,g1] = sgolay(2,bscframeLen1);
[b,g2] = sgolay(2,bscframeLen2);
[b,g3] = sgolay(2,bscframeLen3);
dz=15;
N_O3_bsc=nan(size(ln_bsc));
diff=nan(size(ln_bsc,1),1);
for i=1:length(time_aeroprof)
  diff1= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g1(:,2), 'same');
  diff2= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g2(:,2), 'same');
  diff3= conv(ln_bsc(:,i), factorial(1)/(dz)^1 * g3(:,2), 'same');
  diff(1:bscdh1)=diff1(1:bscdh1);
  diff(bscdh1+1:bscdh2)=diff2(bscdh1+1:bscdh2);
  diff(bscdh2+1:end)=diff3(bscdh2+1:end);
  N_O3_bsc(:,i)=1/(2*d_sigma).*diff;
end

delta_aext=1e-3*((1064/287.2)^ae-(1064/299.1)^ae)*aero_ext;% aerosol extinction in (/m)
D_aext=delta_aext./d_sigma;

%% interpolate the aersosl correction terms to the ozone number density
t_chm15k=time_aeroprof - timediff;% time difference between the local time and utc time (EDT = 4; EST =5)
h_chm15k=rv0; % in km

[X_time,Y_height] = meshgrid(t_chm15k,h_chm15k'); % grid of the original matrix
[Xq_time,Yq_height] = meshgrid(TimeInHour_avg,hkm);% grid of the query matrix: machting the ozone dial
N_O3_bsc_interp=interp2(X_time,Y_height,N_O3_bsc,Xq_time,Yq_height);
D_aext_interp=interp2(X_time,Y_height,D_aext,Xq_time,Yq_height);

%% save the data
chm15kRetrieval.h1=h1;
chm15kRetrieval.m=m;
chm15kRetrieval.study_date=study_date;
chm15kRetrieval.height=rv0;
chm15kRetrieval.time_aeroprof=time_aeroprof;
chm15kRetrieval.time_array_utc=DateTimeArr;

chm15kRetrieval.attbs=attbs;
chm15kRetrieval.attbs_sm=attbs_sm;
chm15kRetrieval.cbh_aero=cbh_aero;
chm15kRetrieval.pbl_aero=pbl_aero;

chm15kRetrieval.cl_chm15=cl_chm15;
chm15kRetrieval.am_interp=am_interp;
chm15kRetrieval.bm_interp=bm_interp;
chm15kRetrieval.lidar_r_array=lidar_r_array;
chm15kRetrieval.lidar_ratio=lidar_r;

chm15kRetrieval.aero_bsa1064=aero_bsa;
chm15kRetrieval.aero_ext1064=aero_ext;
chm15kRetrieval.aero_r1064=aero_r;
chm15kRetrieval.aod1064=aero_aod_cal;
chm15kRetrieval.ae_1020_340=ae_1020_340;
chm15kRetrieval.aod_aeronet=aod_aeronet;


end
% save('20201014_absc_forward.mat','rv0','aero_bsa','aero_ext','aero_r','aero_aod_cal','h1','h2','aod','aod_time','attbs','attbs_sm','bm_00','tm2','time_aeroprof','lidar_r_array','cbh_aero','pbl_aero','cl_chm15','study_date');
