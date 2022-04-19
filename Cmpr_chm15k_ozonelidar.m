%% compare the ceilometer data and the ozone lidar 
[chm15k,path]=newgetbschm15k('choose a file','Documents');

for i=2
    TimeArr=chm15k.AveTime{i};
    DataAve=chm15k.Avedata{i};
    datechosed=chm15k.Date{i};
    cbhAveTime=chm15k.cbh{i};
    
    figure
    I=imagesc(TimeArr,chm15k.Height,DataAve,[0,0.5]);
    set(gca,'YDir','normal');
    colormap(jet);
    colorbar;
    xlabel('UTC Time (hour)','FontSize',13);
    ylabel('Altitude (km)','FontSize',13);
    hold on;
    plot(TimeArr,cbhAveTime,'mo','MarkerSize',4);
    title([datestr(datechosed,'yy/mm/dd'),'CCNY ceilometer RCS (a.u.) at 1064nm']);
    legend('CldBase')
end

id=74;
TimeDiff=4;
chm15kTime_local=TimeArr-TimeDiff+3/60;
ind_t=abs(chm15kTime_local-TimeInHour(id))<0.5/60;
cTime=chm15kTime_local(ind_t);
cProfile=DataAve(:,ind_t);
figure 
plot(cProfile,chm15k.Height,'g');hold on
plot(profile_287_an(:,id),height/1000,'r');
plot(profile_299_an(:,id),height/1000,'b');
legend('CHM15k RCS','287nm','299nm');
xlabel('Signal or Range-corrected Signal (a.u.)');
ylabel('Altitude (km)');
title(['Ozone lidar and CHM15k profile at ',datestr(DateTime(id))])
ylim([0.05,5])