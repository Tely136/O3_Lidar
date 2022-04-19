%% function of CHM15k aerosol backscatter retrieval using the radiometer and aeronet Aod
% load the ceilometer attenuated backscatter coefficient
%clear all; close all;
function chm15kRetrieval=CHM15kForwardRetrieval_new(mwr,aodfile)
[chm15data,pathname]=newgetbschm15k_longrange('choose files','chm15k data');
for i=1
    TimeArr=chm15data.AveTime{i};
    DateTimeArr=chm15data.AveDateTime{i};
    DataAve=chm15data.Avedata{i};
    datechosed=chm15data.Date{i};
    pbltopAveTime=chm15data.pbltop{i};
    pbltfp=chm15data.pbltfp{i};
    pblfp2=chm15data.pblfp2{i};
    cbhAveTime=chm15data.cbh{i};
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
    plot(TimeArr,cbhAveTime,'go','MarkerSize',8);
    title([datestr(chm15data.Date{i},'yyyymmdd'),' CCNY Ceilometer RCS (a.u.) at 1064nm'],'FontSize',13);
    legend('PBLH','CldBase','FontSize',11);
end

ind2=1;


%h1=266;%4km
% h2=333;%5km
% h3=400;%6km
% m=500;%7.5km
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
% folderPath='/Users/Tinker/Documents/MATLAB/ozone lidar';
% mwrfilename='2021-09-27_00-04-06_lv2.csv';
% mwr=loadMWRdata([folderPath,'/',mwrfilename]);
% mwr=loadMWRdata(mwrfile);
 
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
molex_1=sigma1*(mwr.NDAir_m3)*(0.1);   %% unit: km^-1   %% Z/km, ext/km, ncn3, Pmb

hm=mwr.height;  %% Z-km
am=molex_1;  %% mol_ext
S2=(8*pi)/3; 
[X_time,Y_height] = meshgrid(mwr.TimeInHour_utc,hm');
temp_time_array=[mwr.TimeInHour_utc(1);TimeArr(2:end)];
[Xq_time,Yq_height] = meshgrid(temp_time_array,rv0);
am_interp=interp2(X_time,Y_height,am,Xq_time,Yq_height);
bm_interp=am_interp/S2;  % mol.bs at 1064-nm, km-1
% molecular two-way transmittance
 tm2=zeros(size(am_interp))+1;    
 for j=1:length(TimeArr)
     for mk=2:length(rv0)
         taum=trapz(rv0(1:mk),am_interp(1:mk,j));
         tm2(mk,j)=exp(-2*taum);
     end
 end
 %% load AOD data

% aodfile='/Users/Tinker/Documents/MATLAB/Ceilometer/AOD/20210923_20211011_CCNY/20210923_20211011_CCNY.txt';
aod_aeronet=loadAODdata(aodfile);

ind_aod=aod_aeronet.date_utc==datetime(datechosed,'InputFormat','dd-MMM-yyyy');
aod=aod_aeronet.aod_1020(ind_aod);
aod_time=aod_aeronet.TimeInHour_utc(ind_aod);
ae_1020_340=mean(aod_aeronet.ae_1020_340(ind_aod),'omitnan');
% movmean of 30min 
pbl_aero=movmean(chm15data.pbltop{ind2},3,'omitnan');
cbh_aero=cbhAveTime;
%% !!!!!!!!!! dec need to be divided by 10 but oct data no,
attbs=movmean(chm15data.Avedata{ind2}(1:m,:),3,2,'omitnan');
attbs(1:13,:)=ones(13,1)*attbs(13,:);
time_aeroprof=chm15data.AveTime{ind2};


%range movmean <1.5km 100, 1.5-3km 200m, >3km 300m
movnum1=7;
movnum2=13;
movnum3=20;
attbs_sm=nan(m,length(time_aeroprof));
% attbs_sm(1:100,:)=movmean(attbs(1:100,:),movnum1,1);
% attbs_sm(100:200,:)=movmean(attbs(100:200,:),movnum2,1);
% attbs_sm(200:end,:)=movmean(attbs(200:end,:),movnum3,1);

temp=movmean(attbs,movnum1,1);
attbs_sm(1:100,:)=temp(1:100,:);
temp=movmean(attbs,movnum2,1);
attbs_sm(100:200,:)=temp(100:200,:);
temp=movmean(attbs,movnum3,1);
attbs_sm(200:end,:)=temp(200:end,:);



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
            [aero_para,~]=iterative_func2(rv0,attbs_sm(:,i)./cl_chm15,bm_interp(:,i),tm2(:,i),lidar_r(j),zc1,zc2,18.5);
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
    
    
    [aero_para ntimes]=iterative_func2(rv0,attbs_sm(:,i)./cl_chm15,bm_interp(:,i),tm2(:,i),lidar_r_array(i),zc1,zc2,18.5); 
    aero_para(rv0>(cbh_aero(i)-0.12),:)=nan;
    aero_para(rv0>(cbh_aero(max(i-1,1))-0.12)& rv0<(cbh_aero(max(i-1,1))+0.2),:)=nan;
    aero_para(rv0>(cbh_aero(min(i+1,length(time_aeroprof)))-0.12)& rv0<(cbh_aero(min(i+1,length(time_aeroprof)))+0.2),:)=nan;
    aero_zkm(:,i)=aero_para(:,1);
    aero_ext(:,i)=aero_para(:,3);
    aero_bsa(:,i)=aero_para(:,2);
    aero_r(:,i)=1+aero_para(:,2)./bm_interp(:,i);
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

% chm15kRetrieval.chm15k=chm15k;
chm15kRetrieval.aod_aeronet=aod_aeronet;
% chm15kRetrieval.mwr=mwr;

end
% save('20201014_absc_forward.mat','rv0','aero_bsa','aero_ext','aero_r','aero_aod_cal','h1','h2','aod','aod_time','attbs','attbs_sm','bm_00','tm2','time_aeroprof','lidar_r_array','cbh_aero','pbl_aero','cl_chm15','study_date');
