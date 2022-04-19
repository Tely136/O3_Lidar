 % Function iterative_func2(): estimate aerosol backscatter using a forward iterative algorithm
 % Please note that the attenuated backscatter is the calibrated P(z)z2 with the calibration constant
 
  function [aero_para ntimes]=iterative_func2(zkm,attbs_532,bm_532,tm2,S1_aer,Z1_c,Z2_c,S1_cld)
  %% Input:  Zkm-height, attbs_532-attenuated backscatter coefficient, bm_532-molecular backscatter, tm2-molecular transmittance,
  %%         S1_aer: lidar-ratio, Z1_c, Z2_c, S1_cld: cloud-base, top and cloud lidar-ratio
  %% Output: aero_para: zkm, Bs(1/km/sr)   
  %%         ntimes: zkm,counts, rel_diff%
  
  mn=length(zkm);
  bs0=zeros(mn,1);
  bs=zeros(mn,1);
  ext=zeros(mn,1);
  tau=zeros(mn,1);
  
  ta2=zeros(mn,1)+1; % all equal to 1.0 with the aerosol extinction=0.0
  
  ntimes=zeros(mn,3);
  
  % Threshold for ending the iterative calculation: relative difference<0.01% or counts=30
  diff0=0.0001;  % contraint threshold: relative difference 0.01%
  Iterative_counts=30;
  
  dz=zkm(3)-zkm(2); %range-interval in km
  
  % for the 1-st point: 
  % for the forward iterative process, 1-st point is very close to the real value due to the transmittance=1.0  
  tmp1= attbs_532(1)/tm2(1)/ta2(1);
  bs(1)=tmp1-bm_532(1);  
  ext(1)=bs(1)*S1_aer;
  ext1=ext(1);
  %fprintf('1- %6.5f  %6.5f %6.5e\n',ta2(1), ext1, bs(1));  
  tau(1)=ext(1)*zkm(1); % optical depth  
  ta2(1)=exp(-2*tau(1)); % 2-way transmittance   
  tmp1= attbs_532(1)/tm2(1)/ta2(1);
  bs(1)=tmp1-bm_532(1);
  ext(1)=bs(1)*S1_aer;
  ext2=ext(1);  
  %fprintf('2- %6.5f  %6.5f  %6.5e \n',ta2(1), ext2,bs(1));
  
  for i=2:mn  % height      
      if(zkm(i)>=Z1_c & zkm(i)<=Z2_c );
           S1=S1_cld;
      else
           S1=S1_aer;
      end      
        
      tau(i)=tau(i-1)+ext(i)*dz;
      ta2(i)=exp(-2*tau(i));
      tmp1= attbs_532(i)/tm2(i)/ta2(i);
      bs0(i)=tmp1-bm_532(i);               
      ext(i)=bs0(i)*S1;
      
      times0=0; 
      rel_diff=0;
      
      while  1 % forward iterative calculation      
         tau(i)=tau(i-1)+ext(i)*dz;
         ta2(i)=exp(-2*tau(i));          
         tmp1= attbs_532(i)/tm2(i)/ta2(i);
         bs(i)=tmp1-bm_532(i);               
      
         times0=times0+1;         
         if times0>=Iterative_counts % max.=30
             break;
         end
        
         rel_diff=abs(bs(i)-bs0(i))/bs(i);         
         if(rel_diff<=diff0)            
            break;            
         else             
            ext(i)=bs(i)*S1;
            bs0(i)=bs(i);           
         end  %if(rel_diff<=diff0)  
         
      end % while
      
     ntimes(i,:)=[zkm(i),times0,rel_diff*100];
             
  end
  aero_para=[zkm,bs,ext];  

return;
  
  
     
          
      
      
      
      
      
      
      