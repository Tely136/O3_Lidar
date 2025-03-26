% load the radiometer data

function mwr=loadMWRdata4lidar(mwrfilename,N_O3_2,DateTime_avg,height)
%% read the microwave radiometer whole file
opts = delimitedTextImportOptions;
opts.Delimiter=',';
% mwrfilename='2021-09-27_00-04-06_lv2.csv';
rm=readtable(mwrfilename,opts);

%parse the height
h=1e3*cell2mat(cellfun(@str2num,table2array(rm(7,5:end-5)),'UniformOutput',false));% [m] altitude of temperature profile

%parse the code for temperature and pressure
Tcode=cell2mat(cellfun(@str2num,table2array(rm(:,3)),'UniformOutput',false));% code to find temp
ind_temp=(Tcode==401);% rows index for temperature 
ind_p=(Tcode==201);% row index for surface pressure

% parse the time
time_array_utc=datetime(table2array(rm(ind_temp,2)),'inputFormat','MM/dd/yy HH:mm:ss');% time array for radiometer

% change the utc hour to local hour, EDT: -hours(4), EST: -hours(5)
rm_time=time_array_utc-hours(4);

% get the date(utc) of the first data 
mwrdate=datestr(time_array_utc(1),'yyyymmdd');

% parse the temperature and pressure using their code
temp_array=cell2mat(cellfun(@str2num,table2array(rm(ind_temp,5:end-5)),'UniformOutput',false));% temperature profiles (each row)
p_array=100*cell2mat(cellfun(@str2num,table2array(rm(ind_p,6)),'UniformOutput',false));% surface pressure (Pa)

% match the temperature profile with the ozone profile
[~,N]=size(N_O3_2);
temp_avg=nan(length(h),N);
pres0_avg=nan(1,N);
temp_avg_interp=nan(length(height),N);
for i=1:N
    indt=isbetween(rm_time,DateTime_avg(i)-minutes(1),DateTime_avg(i)+minutes(1));
    temp_avg(:,i)=mean(temp_array(indt,:),1,'omitnan')';
    pres0_avg(i)=mean(p_array(indt),'omitnan')';
   temp_avg_interp(:,i)=interp1(h,temp_avg(:,i),height./1000);
    
end 
% P=P0(1-L(h-h0)/T0)^(gM/RL)
% Constant
R=8.31446;%ideal gas constant(J/mol.K)
g=9.80665;% gravitational acceleration (m/s^2)
M=0.0289652;% molar mass of dry air 0.0289652 (kg/mol)
Na=6.02214*10^23;
p = polyfit(h,mean(temp_array,1),1);
L=-p(1);% K/m
T0=p(2);
h0=h(1);
Pres=nan(length(height),N);
ND_air=nan(length(height),N);% m^-3
for i=1:N

    Pres(:,i)=pres0_avg(i).*(1-L*(height-h0)./T0).^(g*M/(R*L));%Pa
    ND_air(:,i)= Na/R.*Pres(:,i)./temp_avg_interp(:,i);
end 
 
mwr.time_array_utc=time_array_utc;
mwr.rm_time=rm_time;
mwr.mwrdate=mwrdate;
mwr.temp_array=temp_array;
mwr.p_array=p_array;
mwr.temp_avg_interp=temp_avg_interp;
mwr.pres0_avg=pres0_avg;
mwr.pres=Pres;
mwr.h=h;
mwr.NDAir_m3=ND_air;
savefileName=['MWR',mwrdate,'.mat'];
save(savefileName,'mwr');
end
