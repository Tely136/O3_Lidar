function aeronet_aod=loadAODdata(aodfilename)
opts = delimitedTextImportOptions;
opts.Delimiter=',';
%aodfilename='/Users/Tinker/Documents/MATLAB/Ceilometer/AOD/20210923_20211011_CCNY/20210923_20211011_CCNY.txt';
tb=readtable(aodfilename,opts);

% parse the time
tb_year=cell2mat(cellfun(@(x) str2double(x(7:10)),table2array(tb(8:end,1)),'UniformOutput',false));
tb_dayofyear=cell2mat(cellfun(@str2double,table2array(tb(8:end,4)),'UniformOutput',false));% days of year

time_array_utc=datetime(tb_year,1,0)+days(tb_dayofyear);
[h,m,s] = hms(time_array_utc);
TimeInHour_utc=h+m./60+s./3600;

% parse the aod of 1020 and angstrom exponent 
% AOD at 1020nm 
aod_1020=cell2mat(cellfun(@str2double,table2array(tb(8:end,6)),'UniformOutput',false));
aod_340=cell2mat(cellfun(@str2double,table2array(tb(8:end,26)),'UniformOutput',false));
aod_380=cell2mat(cellfun(@str2double,table2array(tb(8:end,25)),'UniformOutput',false));

aod_1020(aod_1020<0)=nan;
aod_340(aod_340<0)=nan;
aod_380(aod_380<0)=nan;
% Angstrom exponent 
aod_ae_1020_340=-log(aod_1020./aod_340)./log(1020/340);% aod1/aod2=(w1/w2)^(-a)
aod_ae_380_340=-log(aod_380./aod_340)./log(380/340);% aod1/aod2=(w1/w2)^(-a)

aeronet_aod.time_array_utc=time_array_utc;
aeronet_aod.date_utc=datetime(datestr(time_array_utc,'dd-mmm-yyyy'));
aeronet_aod.TimeInHour_utc=TimeInHour_utc;
aeronet_aod.aod_1020=aod_1020;
aeronet_aod.aod_340=aod_340;
aeronet_aod.ae_1020_340=aod_ae_1020_340;
aeronet_aod.aod_380=aod_380;
aeronet_aod.ae_380_340=aod_ae_380_340;

