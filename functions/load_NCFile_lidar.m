%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read lidar data in the .NC format
% Input: choose lidar NetCdf file
% Output: a structure with fields:
% field:    
% Global Attributes:
%            location            = 'New York City (City College of New York)'
%            latitude            = '40.821 deg North'
%            longitude           = '73.949 deg West'
%            elevation           = '100 m'
%            zenith_angle        = '0 deg'
%            laser0_polarization = 'Linear polarization'
%            laser0_divergence   = '0.5 mrad'
%            laser0_prf          = '30 Hz'
%            telesc0_fov         = '2 mrad'
%            telesc0_diameter    = '0.5 meter'
%            range_resol         = '3.750000e+00 meter'
%            time_resol          = '5 minute'
%            year                = 2018
%            month               = 10
%            day                 = 30
% Variables:
%     time       
%            Size:       77x1
%            Dimensions: time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Lidar_Data_Time in hour, New York Local time on 20181030'
%                        standard_name = 'Lidar_time'
%                        units         = 'hour'
%     range      
%            Size:       2040x1
%            Dimensions: range
%            Datatype:   single
%            Attributes:
%                        units     = 'km'
%                        long_name = 'Range from Telescope to each range gate'
%                        axis      = 'z'
%     abs1064    
%            Size:       1773x77
%            Dimensions: extradim01,time
%            Datatype:   double
%            Attributes:
%                        units         = 'mV*km^2'
%                        long_name     = 'Range-corrected backcatter signals (a.u.) at 1064-nm'
%                        missing_value = -999
%                        wavelength    = 1064
%     ext532     
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'km^(-1)'
%                        long_name     = 'Aerosol extinction coefficient (1/km) at 532-nm'
%                        missing_value = -999
%                        wavelength    = 532
%     ext355     
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'km^(-1)'
%                        long_name     = 'Aerosol extinction coefficient (1/km) at 355-nm'
%                        missing_value = -999
%                        wavelength    = 355
%     ext1064    
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'km^(-1)'
%                        long_name     = 'Aerosol extinction coefficient (1/km) at 1064-nm'
%                        missing_value = -999
%                        wavelength    = 1064
%     bs532      
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'km^(-1)sr^(-1)'
%                        long_name     = 'Aerosol backscatter coefficient (1/km/sr) at 532-nm'
%                        missing_value = -999
%                        wavelength    = 532
%     bs355      
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'km^(-1) sr^(-1)'
%                        long_name     = 'Aerosol backscatter coefficient (1/km/sr) at 355-nm'
%                        missing_value = -999
%                        wavelength    = 355
%     bs1064     
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'km^(-1) sr^(-1)'
%                        long_name     = 'Aerosol backscatter coefficient (1/km/sr) at 1064-nm'
%                        missing_value = -999
%                        wavelength    = 1064
%     extang     
%            Size:       2040x77
%            Dimensions: range,time
%            Datatype:   double
%            Attributes:
%                        units         = 'unitness'
%                        long_name     = 'Extinction-related Angstrom exponent at 532-1064-nm'
%                        missing_value = -999
%                        wavelength    = [532  1064]
%     Sratio     
%            Size:       77x16
%            Dimensions: time,extradim02
%            Datatype:   double
%            Attributes:
%                        units         = 'unitless'
%                        long_name     = 'hour_time,sp_aod532,sp_aod355,sp_aod1064,best_aod532,best_aod355,best_aod1064,best_aod532_pbl,best_aod355_pbl,best_aod1064_pbl,best_aod532_overlap,best_aod355_overlap,best_aod1064_overlap,best_s532,best_s355,best_s1064'
%                        missing_value = -999
%                        wavelength    = [532   355  1064]
%     bm1064     
%            Size:       2040x1
%            Dimensions: range
%            Datatype:   double
%            Attributes:
%                        units         = '1/km/sr'
%                        long_name     = 'Molecular backscatter coefficient (1/km/sr) at at 1064-nm'
%                        missing_value = -999
%                        wavelength    = 1064
%     range_resol
%            Size:       1x1
%            Dimensions: 
%            Datatype:   single
%            Attributes:
%                        units     = 'm'
%                        long_name = 'Range resolution'
%     time_resol 
%            Size:       1x1
%            Dimensions: 
%            Datatype:   single
%            Attributes:
%                        units     = 'minute'
%                        long_name = 'Time resolution or average time'
%     elevation  
%            Size:       1x1
%            Dimensions: 
%            Datatype:   single
%            Attributes:
%                        units     = 'm'
%                        long_name = 'Elevation above mean sea level'
%     latitude   
%            Size:       1x1
%            Dimensions: 
%            Datatype:   single
%            Attributes:
%                        units     = 'degree_north'
%                        valid_min = -90
%                        valid_max = 90
%     longitude  
%            Size:       1x1
%            Dimensions: 
%            Datatype:   single
%            Attributes:
%                        units     = 'degree_east'
%                        valid_min = -180
%                        valid_max = 180
%
% Tested by Dingdong Li at CCNY, Matlab R2018b 
% Last updated on Nov 08, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NC,pathname]=load_NCFile_lidar(instruction,startingfolder)

if nargin==0  % nargin is the number of arguments
    instruction='Pick a NC file to load';
end

if ~(nargin==2)
    startingfolder= cd;
end
%Set up to look for .nc files
[filename, pathname, filterindex] = uigetfile( ...
{  '*.nc', 'WinSpec/32 Data File, version 2.5 header'; ...
   '*.*',  'All Files (*.*)', ... % Add to this to support other files (change filterindex number below)
   },...
   instruction, ...
   startingfolder,...
   'MultiSelect', 'on');

NC.filenames=filename; % Create cell for filename string
NC.bs1064=double(ncread([pathname filename],'bs1064')); % read aerosol backscatter at wavelength of 1064nm
NC.abs1064=ncread([pathname filename],'abs1064'); % read attenuated backscatter at wavelength of 1064nm
NC.bm1064=ncread([pathname filename],'bm1064'); % read molecular backscatter at wavelength of 1064nm
NC.xtime=ncread([pathname filename],'time');% read the local time in hour
NC.Date=datetime(ncreadatt([pathname filename],'/','year'), ncreadatt([pathname filename],'/','month'), ncreadatt([pathname filename],'/','day')); %
NC.zkm=double(ncread([pathname filename],'range'));% The range height in km
k=NC.bs1064<-900;  % filled value: nan
NC.bs1064(k)=NaN;
k=NC.abs1064<-900;  % filled value: nan
NC.abs1064(k)=NaN;
k=NC.bm1064<-900;  % filled value: nan
NC.bm1064(k)=NaN;
NC.Sratio=double(ncread([pathname filename],'Sratio'));% The Sratio
NC.Sratio=NC.Sratio(:,end);
k=NC.Sratio<0;
NC.Sratio(k)=NaN;





