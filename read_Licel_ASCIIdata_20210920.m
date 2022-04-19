% Matlab codes for Reading out the ASCII datafile
% te2191616.444837
% CCNY     16/09/2021 16:44:48 16/09/2021 16:44:48 0100 -073.950000 0040.821000 00.0 00.0
% 0000008 0020 0000000 0020 08 0000000 0020 0000000 0000
% 1 0 1 04000 1 0000 3.75 00287.o 0 0 01 000 16 000008 0.500 BT0
% 1 1 1 04000 1 0000 3.75 00287.o 0 0 00 000 00 000008 3.1746 BC0
% 1 0 1 04000 1 0000 3.75 00299.o 0 0 01 000 16 000008 0.500 BT1
% 1 1 1 04000 1 0000 3.75 00299.o 0 0 00 000 00 000008 3.1746 BC1
% 1 0 1 14000 1 0000 3.75 01064.o 0 0 01 000 16 000008 0.500 BT2
% 1 1 1 14000 1 0000 3.75 01064.o 0 0 00 000 00 000008 3.1746 BC2
% 1 0 1 14000 1 0000 3.75 00387.o 0 0 01 000 16 000008 0.500 BT3
% 1 1 1 14000 1 0000 3.75 00387.o 0 0 00 000 00 000008 3.1746 BC3
% 287.000 .o Analog 0 	287.000 .o Photon Counting 0 	299.000 .o Analog 1 	299.000 .o Photon Counting 1 	1064.000 .o Analog 2 	1064.000 .o Photon Counting 2 	387.000 .o Analog 3 	387.000 .o Photon Counting 3 
% 5.26818	0	5.25864	0	5.10986	0	5.6487	0
% 5.23671	0	5.23289	0	5.11368	0	5.65633	0


clear all; close all;
   % set various constants, plot parameters, etc
   set(0,'DefaultTextFontSize',[15]);   % change font size on all plots
   set(0,'DefaultAxesFontSize',[15]);
   
dz=3.75/1000;  % 

%fl_path='C:\O3DIAL\20210916ozonelidarAlignment\';
fl_path='C:\O3DIAL\20210920_test\';

cd(fl_path);
fl_dirlist=dir('*.txt');
cd ..;

mn=length(fl_dirlist);
ch1_data=[]; ch2_data=[];
ch1ch2_bg=[];

bin_max=4000;
bin_bg=200;
xtime_hr=[];
ch1_pz2=[]; ch2_pz2=[];

zkm=(1:bin_max)*dz;
zkm=zkm';

for i=1:mn
    
   fl_nm=fl_dirlist(i).name
   x= str2num(fl_nm(8:9))+str2num(fl_nm(11:12))/60;
   xtime_hr=[xtime_hr,x];
   xdate=fl_nm(1:7);
   
   fl_nm2=[fl_path fl_nm];
   fileID=fopen(fl_nm2,'r');    
 %  for j=1:12
  for j=1:8
     g=fgetl(fileID); % read the whole 1-st line
            
  end
%    formatSpec = '%f %f %f %f %f %f %f %f\n';
%    [A,count] = fscanf(fileID, formatSpec,[8 bin_max]);
   formatSpec = '%f %f %f %f\n';
   [A,count] = fscanf(fileID, formatSpec,[4 bin_max]);
   data0=A';
   fclose(fileID);

   % calculate the background noise
   data0_bg=nanmean(data0(bin_max-bin_bg:bin_max,:),1);
   
   p1=data0(:,1)-data0_bg(1);
   p2=data0(:,3)-data0_bg(3);
   ch1ch2_bg=[ch1ch2_bg; data0_bg(1) data0_bg(3)];
   ch1_data=[ch1_data,p1];
   ch2_data=[ch2_data,p2];
   ch1_pz2=[ch1_pz2,p1.*zkm.*zkm];
   ch2_pz2=[ch2_pz2,p2.*zkm.*zkm];
end


id=10;

zmax=2.;

% plot the data
figure(1); hold on;
plot(zkm,ch1_data(:,id),'b-'); hold on;
plot(zkm,ch2_data(:,id),'r-'); hold on;
legend('ch1','ch2',1);
%xlim([0 zmax]);
xlim([0 2.4]);
ylim([0 600]);
box on;
grid on;
set(gca,'yscale','log');
ylabel('Lidar signal (mV)');
xlabel('Altitude (km)');
