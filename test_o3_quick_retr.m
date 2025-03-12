folder_path='C:\Users\OzoneLidar\Documents\Lidar Data\20220517\txt';
save_path='C:\Users\OzoneLidar\Documents\MATLAB\ozone lidar\ozone lidar data\ozone lidar results\';
nbin=8000;
bgbins=100;
td=280;
nAvg=10;
start_bin_fr=200;
start_bin_nr=14;
cld_start_bin = 1;% 1km
cld_end_bin = 10;% 10km
[o3_dial_retr]=O3_quick_retrieval(folder_path,save_path,nbin,bgbins,td,...
                                       nAvg,start_bin_fr,start_bin_nr,....
                                       cld_start_bin,cld_end_bin);
                                   