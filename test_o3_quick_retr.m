folderpath='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar/20220429/txt/';
savepath='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar_results/20220429/';
nbin=6000;
bgbins=100;
td=280;
nAvg=10;
start_bin_fr=200;
start_bin_nr=14;
cld_start_bin = 1;% 1km
cld_end_bin = 10;% 10km
[o3_dial_retr]=O3_quick_retrieval(folderpath,savepath,nbin,bgbins,td,...
                                       nAvg,start_bin_fr,start_bin_nr,....
                                       cld_start_bin,cld_end_bin);
                                   