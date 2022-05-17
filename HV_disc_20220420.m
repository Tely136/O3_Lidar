% Analyse the data of different HV and discriminator levels 
%% Load the data files, get time average, dead-time correction
flist = ["750","775","800","825","850","875","875_2"];

folder_path='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar/20220420/';
save_path='/Users/Tinker/Documents/MATLAB/ozonelidar/ozonelidar_repo/ozone_lidar_results/20220420/';
nbin=6000;
dzm=3.75;
bgbins=100;
td=280;
nAvg=5;
for i= 1: length(flist)
    fullpath = strcat(folder_path,flist(i),'txt/');
    [ol{i},save_path]=read_OL_profiles2(fullpath,save_path,nbin,dzm,bgbins,td,nAvg);
end 

for i= 1: length(flist)
    od{i}=importdata(ol{i});
end 
le={"disc: 4","disc: 6","disc: 8","disc: 10","disc: 12","disc: 14"};
for i= 1: length(flist)
%     pc_ad_r287_fr=od{i}.sigprof.pc287/od{i}.sigprof.an287;
%     pc_ad_r287_nr=od{i}.sigprof.pc287nr/od{i}.sigprof.an287nr;
%     pc_ad_r299_fr=od{i}.sigprof.pc299/od{i}.sigprof.an299;
%     pc_ad_r299_nr=od{i}.sigprof.pc299nr/od{i}.sigprof.an299nr;
%     f1=figure('Position', [10 10 900 500]);
%     subplot(2,2,1)
%     plot(od{i}.hkm,pc_ad_r287_fr);set(gca,'LineWidth',1,'FontSize',14);
%     xlabel('Altitude (km)');
%     ylabel('PC/AD ratio');
%     xlim([0,15]);
%     ylim([0,200]);
%     grid on; grid minor
%     legend(le);
%     title(strcat('Far-287 HV(V): ',flist(i)))
%     
%     subplot(2,2,2)
%     plot(od{i}.hkm,pc_ad_r299_fr);set(gca,'LineWidth',1,'FontSize',14);
%     xlabel('Altitude (km)');
%     ylabel('PC/AD ratio');
%     xlim([0,15]);
%     ylim([0,200]);
%     grid on; grid minor
%     legend(le);
%     title(strcat('Far-299 HV(V): ',flist(i)))
%     
%     subplot(2,2,3)
%     plot(od{i}.hkm,pc_ad_r287_nr);set(gca,'LineWidth',1,'FontSize',14);
%     xlabel('Altitude (km)');
%     ylabel('PC/AD ratio');
%     xlim([0,5]);
%     ylim([0,300]);
%     grid on; grid minor
%     legend(le);
%     title(strcat('Near-287 HV(V): ',flist(i)))
%     
%     subplot(2,2,4)
%     plot(od{i}.hkm,pc_ad_r299_nr);set(gca,'LineWidth',1,'FontSize',14);
%     xlabel('Altitude (km)');
%     ylabel('PC/AD ratio');
%     xlim([0,5]);
%     ylim([0,300]);
%     grid on; grid minor
%     legend(le);
%     title(strcat('Near-299 HV(V): ',flist(i)))
%     figname=strcat(save_path,'pc_ad_ratio',flist(i),'.jpg');
%     saveas(f1,figname)

    f1=figure('Position', [10 10 900 500]);
    subplot(2,2,1)
    plot(od{i}.hkm,od{i}.sigprof.pc287);set(gca,'LineWidth',1,'FontSize',14);
    xlabel('Altitude (km)');
    ylabel('PC/AD ratio');
    xlim([0,15]);
    ylim([0,200]);
    grid on; grid minor
    legend(le);
    title(strcat('Far-287 HV(V): ',flist(i)))
    
    subplot(2,2,2)
    plot(od{i}.hkm,od{i}.sigprof.pc287);set(gca,'LineWidth',1,'FontSize',14);
    xlabel('Altitude (km)');
    ylabel('PC/AD ratio');
    xlim([0,15]);
    ylim([0,200]);
    grid on; grid minor
    legend(le);
    title(strcat('Far-299 HV(V): ',flist(i)))
    
    subplot(2,2,3)
    plot(od{i}.hkm,pc_ad_r287_nr);set(gca,'LineWidth',1,'FontSize',14);
    xlabel('Altitude (km)');
    ylabel('PC/AD ratio');
    xlim([0,5]);
    ylim([0,300]);
    grid on; grid minor
    legend(le);
    title(strcat('Near-287 HV(V): ',flist(i)))
    
    subplot(2,2,4)
    plot(od{i}.hkm,pc_ad_r299_nr);set(gca,'LineWidth',1,'FontSize',14);
    xlabel('Altitude (km)');
    ylabel('PC/AD ratio');
    xlim([0,5]);
    ylim([0,300]);
    grid on; grid minor
    legend(le);
    title(strcat('Near-299 HV(V): ',flist(i)))
    figname=strcat(save_path,'pc_ad_ratio',flist(i),'.jpg');
    saveas(f1,figname)

end 