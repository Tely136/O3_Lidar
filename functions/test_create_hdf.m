%% Function: Save CCNY O3DIAL results in the HDF4 format 
% INPUT:
% o3result: O3 result structure
% o3result.nAvg : number of profiles to average [profile]
% o3result.DateTime : date time array of the original time resolution in UTC [yyyy-MMM-dd HH:mm:ss]
% o3result.DateTime_avg : averaged date time array in UTC [yyyy-MMM-dd HH:mm:ss]
% o3result.TimeInHour: time in hour of the original time resolution in UTC [hour]
% o3result.TimeInHour_avg: averaged time in hour  in UTC [hour]
% o3result.alt: altitude [m] 
% o3result.o3_nd: o3 number density combined uncertainty [molec/m^3] 
% o3result.o3_nd_comb_uncer : o3 number density combined uncertainty [molec/m^3] 
% o3result.o3_nd_rand_uncer : o3 number density random (statistical) uncertainty
% o3result.o3_nd_sys_uncer : o3 number density systematic uncertainty
% o3result.o3_nd_vrt_res : o3 number density vertical resolution [m]
% o3result.o3_mr : o3 mixing ratio [ppmv] 
% o3result.o3_mr_comb_uncer: o3 mixing ratio combined uncertainty [ppmv]
% o3result.o3_mr_rand_uncer: o3 mixing ratio random (statistical) uncertainty [ppmv]
% o3result.o3_mr_sys_uncer : o3 mixing ratio systematic uncertainty [ppmv]
% o3result.p: pressure profile [hpa]
% o3result.t: temperature profile [K]
% o3result.p_source: Source of the pressure profile
% o3result.t_source: Source of the temperature profile
% OUTPUT:
% hdf4 file
function test_create_hdf(save_path,o3result)
hdfml('closeall');
DateTime_avg = o3result.DateTime_avg;
DateTime = o3result.DateTime;
nAvg = o3result.nAvg;
start_DateTime = DateTime_avg - minutes(nAvg/2);
stop_DateTime = DateTime_avg + minutes(nAvg/2);

% Modified 2024 Jan 14: fill the sysmatic uncer with 0 and recalculate the
% combined uncer
o3result.o3_nd_comb_uncer = sqrt((o3result.o3_nd_rand_uncer).^2 + fillmissing(o3result.o3_nd_sys_uncer,'constant',0).^2); 
o3result.o3_mr_comb_uncer = sqrt((o3result.o3_mr_rand_uncer).^2 + fillmissing(o3result.o3_mr_sys_uncer,'constant',0).^2); 

% covert the time array to mjd2k 
dt_conv = juliandate(DateTime_avg)-juliandate(datetime('20000101 00:00:00','InputFormat','yyyyMMdd HH:mm:ss'));
% get the start time array for each integrated o3 profiles
% dt_start_conv = juliandate(DateTime(1:nAvg:length(DateTime)))-juliandate(datetime('20000101 00:00:00','InputFormat','yyyyMMdd HH:mm:ss'));
% % get the stop time array for each integrated o3 profiles
% if mod(length(DateTime),nAvg)==0 
%     dt_stop_conv = juliandate(DateTime(nAvg:nAvg:length(DateTime))+minutes(1))-juliandate(datetime('20000101 00:00:00','InputFormat','yyyyMMdd HH:mm:ss'));
% else
%     dt_stop_conv = juliandate(DateTime([nAvg:nAvg:length(DateTime),length(DateTime)])+minutes(1))-juliandate(datetime('20000101 00:00:00','InputFormat','yyyyMMdd HH:mm:ss'));
% end

dt_start_conv = juliandate(start_DateTime)-juliandate(datetime('20000101 00:00:00','InputFormat','yyyyMMdd HH:mm:ss'));
dt_stop_conv = juliandate(stop_DateTime)-juliandate(datetime('20000101 00:00:00','InputFormat','yyyyMMdd HH:mm:ss'));

dt_int_conv= 24.*(dt_stop_conv - dt_start_conv);
% groundbased_lidar.o3_ccny001_hires_nyc.ny_20220218t212713z_20220218t231424_001.hdf
hdf_file_name = sprintf('groundbased_lidar.o3_ccny001_hires_new.york.ny_%st%sz_%st%sz_001.hdf'...
                ,datestr(start_DateTime(1),'yyyymmdd'),datestr(start_DateTime(1),'HHMMss'),...
                datestr(stop_DateTime(end),'yyyymmdd'),datestr(stop_DateTime(end),'HHMMss'));

import matlab.io.hdf4.*

%% Create HDF4 File
flnmo3 = strjoin([save_path,hdf_file_name],'');
sdID = sd.start(flnmo3,'create');
DATA_START_DATE = sprintf('%sT%sZ',datestr(start_DateTime(1),'yyyymmdd'),datestr(start_DateTime(1),'HHMMss'));% 20220218T212713Z
DATA_STOP_DATE = sprintf('%sT%sZ',datestr(stop_DateTime(end),'yyyymmdd'),datestr(stop_DateTime(end),'HHMMss'));%'20220218T231424Z'
FILE_GENERATION_DATE = sprintf('%sT%sZ',datestr(datetime('now'),'yyyymmdd'),datestr(datetime('now'),'HHMMss'));%'20220725T000000Z'

%% Set global attributes
sd.setAttr(sdID,'PI_NAME','Moshary;Fred');
sd.setAttr(sdID,'PI_AFFILIATION','City College of New York;CCNY');
sd.setAttr(sdID,'PI_ADDRESS','160 Convent Ave; New York, NY, 10031; UNITED STATES');
sd.setAttr(sdID,'PI_EMAIL','moshary@ccny.cuny.edu');
sd.setAttr(sdID,'DO_NAME','Wu;Yonghua');
sd.setAttr(sdID,'DO_AFFILIATION','City College of New York;CCNY');
sd.setAttr(sdID,'DO_ADDRESS','160 Convent Ave; New York, NY, 10031; UNITED STATES');
sd.setAttr(sdID,'DO_EMAIL','yhwu@ccny.cuny.edu');
sd.setAttr(sdID,'DS_NAME','Moshary;Fred');
sd.setAttr(sdID,'DS_AFFILIATION','City College of New York;CCNY');
sd.setAttr(sdID,'DS_ADDRESS','160 Convent Ave; New York, NY, 10031; UNITED STATES');
sd.setAttr(sdID,'DS_EMAIL','moshary@ccny.cuny.edu');
sd.setAttr(sdID,'DATA_DESCRIPTION','Quasi-routine tropospheric ozone profiles from CCNY New York Tropospheric O3 Lidar (NYTOL) at CCNY');
sd.setAttr(sdID,'DATA_DISCIPLINE','ATMOSPHERIC.CHEMISTRY;REMOTE.SENSING;GROUNDBASED');
sd.setAttr(sdID,'DATA_GROUP','EXPERIMENTAL;PROFILE.STATIONARY');
sd.setAttr(sdID,'DATA_LOCATION','NEW.YORK.NY');
sd.setAttr(sdID,'DATA_SOURCE','LIDAR.O3_CCNY001_HIRES');
sd.setAttr(sdID,'DATA_VARIABLES','LATITUDE.INSTRUMENT;LONGITUDE.INSTRUMENT;ALTITUDE.INSTRUMENT;DATETIME;DATETIME.START;DATETIME.STOP;INTEGRATION.TIME;ALTITUDE;O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL;O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD;O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD;O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD;O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM;O3.MIXING.RATIO.VOLUME_DERIVED;O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD;O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD;O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD;PRESSURE_INDEPENDENT;TEMPERATURE_INDEPENDENT;PRESSURE_INDEPENDENT_SOURCE;TEMPERATURE_INDEPENDENT_SOURCE');
sd.setAttr(sdID,'DATA_START_DATE',DATA_START_DATE);
sd.setAttr(sdID,'DATA_STOP_DATE',DATA_STOP_DATE);
sd.setAttr(sdID,'DATA_FILE_VERSION','001');
sd.setAttr(sdID,'DATA_MODIFICATIONS','')
sd.setAttr(sdID,'DATA_CAVEATS','');
sd.setAttr(sdID,'DATA_RULES_OF_USE','')
sd.setAttr(sdID,'DATA_ACKNOWLEDGEMENT','Notify PI that data is being used and ask for proper form of acknowledgement');
sd.setAttr(sdID,'DATA_QUALITY','Nominal'); 
sd.setAttr(sdID,'DATA_TEMPLATE','GEOMS-TE-LIDAR-O3-005');
sd.setAttr(sdID,'DATA_PROCESSOR',''); 
sd.setAttr(sdID,'FILE_NAME',hdf_file_name);
sd.setAttr(sdID,'FILE_GENERATION_DATE',FILE_GENERATION_DATE);
sd.setAttr(sdID,'FILE_ACCESS','TOLNET');
sd.setAttr(sdID,'FILE_PROJECT_ID','');
sd.setAttr(sdID,'FILE_ASSOCIATION','')
sd.setAttr(sdID,'FILE_DOI','');
sd.setAttr(sdID,'FILE_META_VERSION','04R062;CUSTOM');


%% Create HDF4 Data Set

% LATITUDE.INSTRUMENT
lat=single(40.821);
psID0 = sd.create(sdID,'LATITUDE.INSTRUMENT','single',1);
sd.writeData(psID0,lat);
VAR_NAME = 'LATITUDE.INSTRUMENT';
VAR_DESCRIPTION = 'Instrument geolocation (+ for north; - for south)';
VAR_NOTES = ' ';
VAR_SIZE = '1';
VAR_DEPEND = 'CONSTANT';
VAR_DATA_TYPE = "REAL";
VAR_UNITS = 'deg';
VAR_SI_CONVERSION= '0.0;1.74533E-2;rad';
VAR_VALID_MIN = -90.0;
VAR_VALID_MAX = 90.0;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID0,'VAR_NAME',VAR_NAME);
sd.setAttr(psID0,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID0,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID0,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID0,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID0,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID0,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID0,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID0,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID0,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID0,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID0);

% LONGITUDE.INSTRUMENT
lon=single(-73.948);
psID1 = sd.create(sdID,'LONGITUDE.INSTRUMENT','single',1);
sd.writeData(psID1,lon);
VAR_NAME = 'LONGITUDE.INSTRUMENT';
VAR_DESCRIPTION = 'Instrument geolocation (+ for east; - for west)';
VAR_NOTES = ' ';
VAR_SIZE = '1';
VAR_DEPEND = 'CONSTANT';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'deg';
VAR_SI_CONVERSION= '0.0;1.74533E-2;rad';
VAR_VALID_MIN = -180.0;
VAR_VALID_MAX = 180.0;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID1,'VAR_NAME',VAR_NAME);
sd.setAttr(psID1,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID1,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID1,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID1,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID1,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID1,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID1,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID1,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID1,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID1,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID1);

% ALTITUDE.INSTRUMENT
alt_intrm=single(90);
psID2 = sd.create(sdID,'ALTITUDE.INSTRUMENT','single',1);
sd.writeData(psID2,alt_intrm);
VAR_NAME = 'ALTITUDE.INSTRUMENT';
VAR_DESCRIPTION = 'Instrument altitude above sea level (m)';
VAR_NOTES = ' ';
VAR_SIZE = '1';
VAR_DEPEND = 'CONSTANT';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'm';
VAR_SI_CONVERSION= '0.0;1.0;m';
VAR_VALID_MIN = -300.0;
VAR_VALID_MAX = 20000.0;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID2,'VAR_NAME',VAR_NAME);
sd.setAttr(psID2,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID2,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID2,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID2,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID2,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID2,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID2,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID2,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID2,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID2,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID2);

% DATETIME
dt=fillmissing(dt_conv,'constant',-9999.0);
psID3 = sd.create(sdID,'DATETIME','double',length(dt));
sd.writeData(psID3,dt);
VAR_NAME = 'DATETIME';
VAR_DESCRIPTION = 'Weighted datetime of the LIDAR measurement.';
VAR_NOTES = ' ';
VAR_SIZE = num2str(length(dt),'%d');
VAR_DEPEND = 'DATETIME';
VAR_DATA_TYPE = 'DOUBLE';
VAR_UNITS = 'MJD2K';
VAR_SI_CONVERSION= '0.0;86400.0;s';
VAR_VALID_MIN = -2375.0;
VAR_VALID_MAX = 365243.0;
VAR_FILL_VALUE = -9999.0;
sd.setAttr(psID3,'VAR_NAME',VAR_NAME);
sd.setAttr(psID3,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID3,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID3,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID3,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID3,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID3,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID3,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID3,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID3,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID3,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID3);

%DATETIMET.START
dt_start=fillmissing(dt_start_conv,'constant',-9999.0);
psID4 = sd.create(sdID,'DATETIME.START','double',length(dt_start));
sd.writeData(psID4,dt_start);
VAR_NAME = 'DATETIME.START';
VAR_DESCRIPTION = 'Start date and time of measurement (MJD2000).';
VAR_NOTES = ' ';
VAR_SIZE = num2str(length(dt_start),'%d');
VAR_DEPEND = 'DATETIME';
VAR_DATA_TYPE = 'DOUBLE';
VAR_UNITS = 'MJD2K';
VAR_SI_CONVERSION= '0.0;86400.0;s';
VAR_VALID_MIN = -2375.0;
VAR_VALID_MAX = 365243.0;
VAR_FILL_VALUE = -9999.0;
sd.setAttr(psID4,'VAR_NAME',VAR_NAME);
sd.setAttr(psID4,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID4,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID4,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID4,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID4,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID4,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID4,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID4,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID4,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID4,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID4);

%DATETIME.STOP
dt_stop=fillmissing(dt_stop_conv,'constant',-9999.0);
psID5 = sd.create(sdID,'DATETIME.STOP','double',length(dt_stop));
sd.writeData(psID5,dt_stop);
VAR_NAME = 'DATETIME.STOP';
VAR_DESCRIPTION = 'End date and time of measurement (MJD2000).';
VAR_NOTES = ' ';
VAR_SIZE = num2str(length(dt_stop),'%d');
VAR_DEPEND = 'DATETIME';
VAR_DATA_TYPE = 'DOUBLE';
VAR_UNITS = 'MJD2K';
VAR_SI_CONVERSION= '0.0;86400.0;s';
VAR_VALID_MIN = -2375.0;
VAR_VALID_MAX = 365243.0;
VAR_FILL_VALUE = -9999.0;
sd.setAttr(psID5,'VAR_NAME',VAR_NAME);
sd.setAttr(psID5,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID5,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID5,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID5,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID5,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID5,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID5,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID5,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID5,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID5,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID5);

% INTEGRATION.TIME
dt_int=single(fillmissing(dt_int_conv,'constant',-9999.0));
psID6 = sd.create(sdID,'INTEGRATION.TIME','single',length(dt_int));
sd.writeData(psID6,dt_int);
VAR_NAME = 'INTEGRATION.TIME';
VAR_DESCRIPTION = 'Effective integration time of measurement (hours).';
VAR_NOTES = ' ';
VAR_SIZE = num2str(length(dt_int),'%d');
VAR_DEPEND = 'DATETIME';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'h';
VAR_SI_CONVERSION= '0.0;3600.0;s';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX = 24.0;
VAR_FILL_VALUE = -9999.0;
sd.setAttr(psID6,'VAR_NAME',VAR_NAME);
sd.setAttr(psID6,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID6,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID6,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID6,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID6,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID6,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID6,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID6,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID6,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID6,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID6);

% ALTITUDE

alt=single(fillmissing(o3result.alt,'constant',-999.0));
psID7 = sd.create(sdID,'ALTITUDE','single',length(alt));
sd.writeData(psID7,alt);
VAR_NAME = 'ALTITUDE';
VAR_DESCRIPTION = 'Altitude above sea level (m)';
VAR_NOTES = 'Actual altitude measurement grid (amsl)';
VAR_SIZE = num2str(length(alt),'%d');
VAR_DEPEND = 'ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'm';
VAR_SI_CONVERSION= '0.0;1.0;m';
VAR_VALID_MIN = -300.0;
VAR_VALID_MAX = 120000.0;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID7,'VAR_NAME',VAR_NAME);
sd.setAttr(psID7,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID7,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID7,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID7,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID7,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID7,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID7,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID7,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID7,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID7,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID7);

% O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL
o3_nd_abs_diff=single(fillmissing(o3result.o3_nd,'constant',-9.99e+19));
psID8 = sd.create(sdID,'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL','single',size(o3_nd_abs_diff));
sd.writeData(psID8,o3_nd_abs_diff);
VAR_NAME = 'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL';
VAR_DESCRIPTION = 'Ozone number density (molec m-3)';
VAR_NOTES = 'This is the native lidar measurement, to be used preferably over mixing ratio';
VAR_SIZE = sprintf('%d;%d',size(o3_nd_abs_diff,2),size(o3_nd_abs_diff,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'molec m-3';
VAR_SI_CONVERSION= '0.0;1.66054E-24;mol m-3';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  2.5e+19;
VAR_FILL_VALUE = -9.99e+19;
sd.setAttr(psID8,'VAR_NAME',VAR_NAME);
sd.setAttr(psID8,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID8,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID8,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID8,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID8,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID8,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID8,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID8,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID8,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID8,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID8);


% O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD
o3_nd_abs_diff_uncr_rand=single(fillmissing(o3result.o3_nd_rand_uncer,'constant',-999.0));
psID10 = sd.create(sdID,'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD','single',size(o3_nd_abs_diff_uncr_rand));
sd.writeData(psID10,o3_nd_abs_diff_uncr_rand);
VAR_NAME = 'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.RANDOM.STANDARD';
VAR_DESCRIPTION = 'Random uncertainty: NDACC-lidar-standardized definition';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_nd_abs_diff_uncr_rand,2),size(o3_nd_abs_diff_uncr_rand,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'molec m-3';
VAR_SI_CONVERSION= '0.0;1.66054E-24;mol m-3';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  2.0e+19;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID10,'VAR_NAME',VAR_NAME);
sd.setAttr(psID10,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID10,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID10,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID10,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID10,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID10,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID10,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID10,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID10,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID10,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID10);

% O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD
o3_nd_abs_diff_uncr_sys=single(fillmissing(o3result.o3_nd_sys_uncer,'constant',-999.0));
psID11 = sd.create(sdID,'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD','single',size(o3_nd_abs_diff_uncr_sys));
sd.writeData(psID11,o3_nd_abs_diff_uncr_sys);
VAR_NAME = 'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.SYSTEMATIC.STANDARD';
VAR_DESCRIPTION = 'Systematic uncertainty: NDACC-lidar-standardized definition';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_nd_abs_diff_uncr_sys,2),size(o3_nd_abs_diff_uncr_sys,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'molec m-3';
VAR_SI_CONVERSION= '0.0;1.66054E-24;mol m-3';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  2.0e+19;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID11,'VAR_NAME',VAR_NAME);
sd.setAttr(psID11,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID11,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID11,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID11,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID11,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID11,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID11,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID11,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID11,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID11,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID11);

% O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD
o3_nd_abs_diff_uncr_comb=single(fillmissing(o3result.o3_nd_comb_uncer,'constant',-999.0));
psID9 = sd.create(sdID,'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD','single',size(o3_nd_abs_diff_uncr_comb));
sd.writeData(psID9,o3_nd_abs_diff_uncr_comb);
VAR_NAME = 'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD';
VAR_DESCRIPTION = 'Combined uncertainty: NDACC-lidar-standardized definition';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_nd_abs_diff_uncr_comb,2),size(o3_nd_abs_diff_uncr_comb,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'molec m-3';
VAR_SI_CONVERSION= '0.0;1.66054E-24;mol m-3';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  2.0e+19;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID9,'VAR_NAME',VAR_NAME);
sd.setAttr(psID9,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID9,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID9,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID9,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID9,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID9,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID9,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID9,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID9,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID9,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID9);

% O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM
o3_nd_abs_diff_res_alt=single(fillmissing(o3result.o3_nd_vrt_res,'constant',-999.0));
psID12 = sd.create(sdID,'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM','single',size(o3_nd_abs_diff_res_alt));
sd.writeData(psID12,o3_nd_abs_diff_res_alt);
VAR_NAME = 'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE.FWHM';
VAR_DESCRIPTION = 'Vertical resolution: Full-width at half-maximum (FWHM) of impulse response';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_nd_abs_diff_res_alt,2),size(o3_nd_abs_diff_res_alt,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'm';
VAR_SI_CONVERSION= '0.0;1.0;m';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  25000.0;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID12,'VAR_NAME',VAR_NAME);
sd.setAttr(psID12,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID12,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID12,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID12,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID12,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID12,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID12,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID12,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID12,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID12,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID12);

% O3.MIXING.RATIO.VOLUME_DERIVED
o3_mr=single(fillmissing(o3result.o3_mr,'constant',-999.0));
psID13 = sd.create(sdID,'O3.MIXING.RATIO.VOLUME_DERIVED','single',size(o3_mr));
sd.writeData(psID13,o3_mr);
VAR_NAME = 'O3.MIXING.RATIO.VOLUME_DERIVED';
VAR_DESCRIPTION = 'Derived Ozone mixing ratio';
VAR_NOTES = 'Mixing ratio calculated from ozone number density profile.';
VAR_SIZE = sprintf('%d;%d',size(o3_mr,2),size(o3_mr,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'ppmv';
VAR_SI_CONVERSION= '0.0;1.0E-6;1';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  20;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID13,'VAR_NAME',VAR_NAME);
sd.setAttr(psID13,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID13,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID13,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID13,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID13,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID13,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID13,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID13,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID13,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID13,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID13);

% O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD
o3_mr_uncr_comb=single(fillmissing(o3result.o3_mr_comb_uncer,'constant',-999.0));
psID14 = sd.create(sdID,'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD','single',size(o3_mr_uncr_comb));
sd.writeData(psID14,o3_mr_uncr_comb);
VAR_NAME = 'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD';
VAR_DESCRIPTION = 'Uncertainty: NDACC-lidar-standardized definition';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_mr_uncr_comb,2),size(o3_mr_uncr_comb,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'ppmv';
VAR_SI_CONVERSION= '0.0;1.0E-6;1';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  20;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID14,'VAR_NAME',VAR_NAME);
sd.setAttr(psID14,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID14,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID14,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID14,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID14,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID14,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID14,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID14,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID14,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID14,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID14);

% O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD
o3_mr_uncr_rand=single(fillmissing(o3result.o3_mr_rand_uncer,'constant',-999.0));
psID15 = sd.create(sdID,'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD','single',size(o3_mr_uncr_rand));
sd.writeData(psID15,o3_mr_uncr_rand);
VAR_NAME = 'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.RANDOM.STANDARD';
VAR_DESCRIPTION = 'Random uncertainty: NDACC-lidar-standardized definition';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_mr_uncr_rand,2),size(o3_mr_uncr_rand,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'ppmv';
VAR_SI_CONVERSION= '0.0;1.0E-6;1';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  20;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID15,'VAR_NAME',VAR_NAME);
sd.setAttr(psID15,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID15,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID15,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID15,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID15,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID15,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID15,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID15,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID15,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID15,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID15);


% O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD
o3_mr_uncr_sys=single(fillmissing(o3result.o3_mr_sys_uncer,'constant',-999.0));
psID16 = sd.create(sdID,'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD','single',size(o3_mr_uncr_sys));
sd.writeData(psID16,o3_mr_uncr_sys);
VAR_NAME = 'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.SYSTEMATIC.STANDARD';
VAR_DESCRIPTION = 'Systematic uncertainty: NDACC-lidar-standardized definition';
VAR_NOTES = 'See DOI:10.5194/amt-9-4051-2016 and http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm';
VAR_SIZE = sprintf('%d;%d',size(o3_mr_uncr_sys,2),size(o3_mr_uncr_sys,1));
VAR_DEPEND = 'DATETIME;ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'ppmv';
VAR_SI_CONVERSION= '0.0;1.0E-6;1';
VAR_VALID_MIN = 0.0;
VAR_VALID_MAX =  20;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID16,'VAR_NAME',VAR_NAME);
sd.setAttr(psID16,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID16,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID16,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID16,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID16,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID16,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID16,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID16,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID16,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID16,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID16);

%  PRESSURE_INDEPENDENT
p=single(fillmissing(o3result.p','constant',-999.0));
psID17 = sd.create(sdID,'PRESSURE_INDEPENDENT','single',size(p));
sd.writeData(psID17,p);
VAR_NAME = 'PRESSURE_INDEPENDENT';
VAR_DESCRIPTION = 'Pressure profile used to derive mixing ratio';
VAR_NOTES = 'May be a single source, or a combination of multiple sources';
VAR_SIZE = sprintf('%d',length(p));
VAR_DEPEND = 'ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'hPa';
VAR_SI_CONVERSION= '0.0;1.0E2;kg m-1 s-2';
VAR_VALID_MIN = 0.0001;
VAR_VALID_MAX =  1100;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID17,'VAR_NAME',VAR_NAME);
sd.setAttr(psID17,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID17,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID17,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID17,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID17,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID17,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID17,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID17,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID17,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID17,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID17);


% TEMPERATURE_INDEPENDENT
t=single(fillmissing(o3result.t','constant',-999.0));
psID18 = sd.create(sdID,'TEMPERATURE_INDEPENDENT','single',size(t));
sd.writeData(psID18,t);
VAR_NAME = 'TEMPERATURE_INDEPENDENT';
VAR_DESCRIPTION = 'Temperature profile used to derive mixing ratio';
VAR_NOTES = 'May be a single source, or a combination of multiple sources';
VAR_SIZE = sprintf('%d',length(t));
VAR_DEPEND = 'ALTITUDE';
VAR_DATA_TYPE = 'REAL';
VAR_UNITS = 'K';
VAR_SI_CONVERSION= '0.0;1.0;K';
VAR_VALID_MIN = 50.0;
VAR_VALID_MAX =  450.0;
VAR_FILL_VALUE = -999.0;
sd.setAttr(psID18,'VAR_NAME',VAR_NAME);
sd.setAttr(psID18,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID18,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID18,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID18,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID18,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID18,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID18,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID18,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID18,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID18,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID18);

% PRESSURE_INDEPENDENT_SOURCE
p_source=repelem(o3result.p_source,length(alt),1)';
psID19 = sd.create(sdID,'PRESSURE_INDEPENDENT_SOURCE','char',size(p_source));
sd.writeData(psID19,p_source);
VAR_NAME = 'PRESSURE_INDEPENDENT_SOURCE';
VAR_DESCRIPTION = 'Pressure profile source (e.g. Lidar; NCEP; Sonde; ECMWF etc.)';
VAR_NOTES = 'May be a single source, or a combination of multiple sources';
VAR_SIZE = sprintf('%d',length(alt));
VAR_DEPEND = 'ALTITUDE';
VAR_DATA_TYPE = 'STRING';
VAR_UNITS = ' ';
VAR_SI_CONVERSION= ' ';
VAR_VALID_MIN = '';
VAR_VALID_MAX =  '';
VAR_FILL_VALUE = '';
sd.setAttr(psID19,'VAR_NAME',VAR_NAME);
sd.setAttr(psID19,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID19,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID19,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID19,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID19,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID19,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID19,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID19,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID19,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID19,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID19);

% % TEMPERATURE_INDEPENDENT_SOURCE
t_source=repelem(o3result.t_source,size(alt,1),1)';
psID20 = sd.create(sdID,'TEMPERATURE_INDEPENDENT_SOURCE','char',size(t_source));
sd.writeData(psID20,t_source);
VAR_NAME = 'TEMPERATURE_INDEPENDENT_SOURCE';
VAR_DESCRIPTION = 'Temperature profile source (e.g. Lidar; NCEP; Sonde; ECMWF etc.)';
VAR_NOTES = 'May be a single source, or a combination of multiple sources';
VAR_SIZE = sprintf('%d',length(alt));
VAR_DEPEND = 'ALTITUDE';
VAR_DATA_TYPE = 'STRING';
VAR_UNITS = ' ';
VAR_SI_CONVERSION= ' ';
VAR_VALID_MIN = '';
VAR_VALID_MAX =  '';
VAR_FILL_VALUE = '';
sd.setAttr(psID20,'VAR_NAME',VAR_NAME);
sd.setAttr(psID20,'VAR_DESCRIPTION', VAR_DESCRIPTION);
sd.setAttr(psID20,'VAR_NOTES', VAR_NOTES);
sd.setAttr(psID20,'VAR_SIZE',VAR_SIZE);
sd.setAttr(psID20,'VAR_DEPEND',VAR_DEPEND);
sd.setAttr(psID20,'VAR_DATA_TYPE',VAR_DATA_TYPE);
sd.setAttr(psID20,'VAR_UNITS',VAR_UNITS);
sd.setAttr(psID20,'VAR_SI_CONVERSION',VAR_SI_CONVERSION);
sd.setAttr(psID20,'VAR_VALID_MIN',VAR_VALID_MIN);
sd.setAttr(psID20,'VAR_VALID_MAX',VAR_VALID_MAX);
sd.setAttr(psID20,'VAR_FILL_VALUE',VAR_FILL_VALUE);
sd.endAccess(psID20);

%%Test
%% Close HDF4 File
sd.close(sdID);
end 