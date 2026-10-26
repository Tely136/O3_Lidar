clc; clearvars;

root_dir = "C:\Users\tely1\PhD";

addpath(genpath(fullfile(root_dir, "O3_Lidar\")));

% DIAL Processing Functions
% addpath(genpath(fullfile(root_dir, "O3_Lidar\proc_functions")));

% Signal Analysis and Modeling Functions
addpath(genpath(fullfile(root_dir, "signal_analysis\")));
%%
start_day       = datetime(datetime("2024-10-01", "InputFormat", "uuuu-MM-dd"), "InputFormat", "uuuu-MM-dd"); %[control:datePicker:8a9f]{"position":[28,79]}
end_day         = datetime(datetime("2024-10-31", "InputFormat", "uuuu-MM-dd"), "InputFormat", "uuuu-MM-dd"); %[control:datePicker:4f06]{"position":[28,79]}
data_path       = fullfile(root_dir, "lidar_data");
result_path     = fullfile(root_dir, "O3_processed_data");
log_path        = fullfile(root_dir, "O3_processed_data");

log_file        = fullfile(log_path, "o3_process.log");

if ~exist(result_path, "dir")
    mkdir(result_path)
end

if ~exist(log_path, "dir")
    mkdir(log_path)
end

if ~exist(log_file, "file")
    fid = fopen(log_file, "a");
    if fid ~= -1
        fclose(fid);
    end
end

folders = dir(data_path);
folders = folders([folders.isdir]);

for i = 1:size(folders, 1)
    tmp_folder = folders(i).name;

    try
        folderDate = datetime(tmp_folder, "InputFormat", "uuuuMMdd", "Format", "uuuuMMdd");
    catch
        continue
    end

    if folderDate >= start_day && folderDate <= end_day
        fprintf("processing %s\n", tmp_folder);
        output_filename = sprintf("o3_processed_%04d%02d%02d.mat", year(folderDate), month(folderDate), day(folderDate));

        try
            out = process(fullfile(data_path, tmp_folder), fullfile(result_path, output_filename));
        catch ME
            warning("Processing failed for %s: %s", tmp_folder, ME.message);

            fid = fopen(log_file, "a");
            if fid ~= -1
                fprintf(fid, "%s | %s | %s | %s\n", ...
                    char(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
                    tmp_folder, ...
                    "FAILED", ...
                    ME.message);
                fclose(fid);
            end
        end
    end
end
%%
plt = "mg"; %[control:dropdown:11a3]{"position":[7,11]}
switch plt
    case "an"
        q_nr_plt = out.q_o3_an_nr;
        q_fr_plt = out.q_o3_an_fr;
        q_merge_plt = out.q_o3_an_merged;
    case "pc"
        q_nr_plt = out.q_o3_pc_nr;
        q_fr_plt = out.q_o3_pc_fr;
        q_merge_plt = out.q_o3_pc_merged;
    case "mg"
        q_nr_plt = out.q_o3_mg_nr;
        q_fr_plt = out.q_o3_mg_fr;
        q_merge_plt = out.q_o3_mg_merged;
end

times_avg = out.times_avg;
z_km_vec = out.z_vec.*1e-3;

ND_plot(times_avg,z_km_vec,q_nr_plt,[0 120],[0 3],"","Altitude (km)", "NR Ozone Mixing Ratio", "",12)
ND_plot(times_avg,z_km_vec,q_fr_plt,[0 120],[0 10],"","Altitude (km)", "FR Ozone Mixing Ratio", "",12)
ND_plot(times_avg,z_km_vec,q_merge_plt,[0 120],[0 10],"","Altitude (km)", "Ozone Concentration (ppb)", "",12)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":18.6}
%---
%[control:datePicker:8a9f]
%   data: {"defaultValue":"2026-09-15","displayFormat":"uuuu-MM-dd","label":"Date picker","run":"Nothing"}
%---
%[control:datePicker:4f06]
%   data: {"defaultValue":"2026-10-15","displayFormat":"uuuu-MM-dd","label":"Date picker","run":"Nothing"}
%---
%[control:dropdown:11a3]
%   data: {"defaultValue":"\"an\"","itemLabels":["an","pc","mg"],"items":["\"an\"","\"pc\"","\"mg\""],"label":"plt","run":"Section"}
%---
