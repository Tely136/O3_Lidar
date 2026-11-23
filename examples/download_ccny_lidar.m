clc;clearvars;
%%
root_dir = "C:\Users\tely1\PhD";
addpath(genpath(fullfile(root_dir, "O3_Lidar\")));

data_path = fullfile(root_dir, "lidar_data");

dates = datetime("2026-01-01", "InputFormat", "uuuu-MM-dd"):datetime("2027-01-01", "InputFormat", "uuuu-MM-dd"); %[control:datePicker:551c]{"position":[9,60]} %[control:datePicker:93cd]{"position":[61,112]}

for i = 1:length(dates) %[output:group:8883e5a9]
    download_ccny_lidar_data(dates(i),data_path) %[output:6a80fc8c]
end %[output:group:8883e5a9]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[control:datePicker:551c]
%   data: {"defaultValue":"2026-10-01","displayFormat":"dd-MMM-uuuu","label":"dates","run":"Nothing"}
%---
%[control:datePicker:93cd]
%   data: {"defaultValue":"2026-10-31","displayFormat":"dd-MMM-uuuu","label":"Date picker","run":"Nothing"}
%---
%[output:6a80fc8c]
%   data: {"dataType":"text","outputData":{"text":"pathconf: Not a directory\n\npathconf: Not a directory\n\npathconf: Not a directory\n\n","truncated":false}}
%---
