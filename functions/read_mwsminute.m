%% Read weather station data
% # Function: read data the csv files from the MSWminute.csv
% #           time, temp, windspeed, RH, wind direction
% Input: 
% fullfile: the full file name including the file path and name
% Output:
% weather_stat: weather data structure
%           field:
%           wv : average wind speed [m/s]
%           wd : wind direction from true north [deg]
%           tc : air temperature [deg C]
%           rh : relative humidity [%]
%           dt : datetime [local time]
%% 
% # DATA DATE REPRESENTATION
% # Older data, prior to 01/01/2021:
% # 	Data (observation date and time) is collected and displayed in EST (UTC-5:00)
% # Newer data, starting on 01/01/2021:
% # 	Data (observation date and time) is colleected and displayed in UTC (UTC-0:00)
% #
% # For files in named format of YYYYMMDD"MWSminute": Data Columns for MWS (Marshak Weather Station) date format
% #0 A.	"101" ingore
% #1 B.	Year of data observation
% #2 C.	DOY - calender day of data observation
% #3 D.	Time - Hour and minute of data observation
% #4 E.	N/A
% #5 F.	Solar, w/m2, Instantaneous
% #6 G.	Solar, w/m2, Average
% #7 H.	Solar, w/m2, Maximum
% #8 I.	Solar, w/m2, Minimum
% #9 J.	Solar, w/m2, Standard Deviation
% #10 K.	Barometric Pressure in kP
% #11 L.	Wind Speed, Instantaneous, m/s
% #12 M.	Wind Speed, Average, m/s
% #13 N.	Wind Speed, Maximum, m/s
% #14 O.	Wind Speed, Minimum, m/s
% #15 P.	Wind Speed, Standard Deviation, m/s
% #16 Q.	Wind Direction, deg from true North
% #17 R.	Wind Speed, m/s, S_WVT
% #18 S.	Wind Speed, m/s, U_WVT
% #19 T.	Wind Direction, m/s, DU_WVT
% #20 U.	Wind Direction, deg, SDU_WVT
% #21 V.	Air Temperature, deg C.
% #22 W.	Relative Humidity, %
% #23 X.	Dew Point, deg C.
% #24 Y.	Heat Index, deg C.
%%
function weather_stat= read_mwsminute(fullfile)
% path = '/Users/Tinker/Documents/MATLAB/Ceilometer/WeatherStation/';
% filename = '20220606MWSminute.csv';
% fullfile = [path,filename];
T = readtable(fullfile,'Delimiter',',','ReadVariableNames',false);
yr=T{:,2};
doy = T{:,3};
mint= mod(T{:,4},100);
hrs = floor(T{:,4}./100);
mon = month(datetime(yr(1),1,0) + days(doy(1)));
if mon >=11||mon<3
    dt_diff = 5; % if month is 11, 12, 1, 2 convert from utc to est: edt= utc -4
else
    dt_diff =4; % otherwise convert to est
end
dt =NaT(size(yr));
utcdt =NaT(size(yr));
for ii=1:length(yr)
    utcdt(ii) = datetime(yr(ii),1,0) + days(doy(ii))+hours(hrs(ii))+minutes(mint(ii));
    dt(ii) = datetime(yr(ii),1,0) + days(doy(ii))+hours(hrs(ii)-dt_diff)+minutes(mint(ii));
end
weather_stat.wv = T{:,13};
weather_stat.wd = T{:,17};
weather_stat.tc = T{:,22};
weather_stat.rh = T{:,23};
weather_stat.dt = dt;
weather_stat.utc = utcdt;
