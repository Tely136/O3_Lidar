function plot_cld_screen_result(DateTime_avg,TimeInHour_avg,hkm,pz2,d_Pz2,cldBaseZ,cldBaseZ_ind,cldTopZ,cldTopZ_ind)
%% cloud screen test  result
figure
subplot(1,2,1)
I=imagesc(TimeInHour_avg,hkm,pz2); hold on
plot(TimeInHour_avg,cldBaseZ(1,:),'o');
plot(TimeInHour_avg,cldTopZ(1,:),'^');
ylim([0,12]);
set(gca,'YDir','normal','FontSize',14,'ColorScale','log');colormap('jet');colorbar
xlabel('Local Time (hour)');ylabel('Altitude (km)')
title(['Ozone lidar smoothed Pz2(a.u.) at ',datestr(DateTime_avg(1),'yy/mm/dd')]);

subplot(1,2,2)
I=imagesc(TimeInHour_avg,hkm(1:end-1),d_Pz2);ylim([0,12]);
set(gca,'YDir','normal','FontSize',14,'ColorScale','log');colormap('jet');colorbar
xlabel('Local Time (hour)');ylabel('Altitude (km)')
title(['Ozone lidar smoothed P derivative (a.u.) at ',datestr(DateTime_avg(1),'yy/mm/dd')]);

%% Test--- plot each cloud profile
% id_all = find(~isnan(cldBaseZ(1,:)));
% idh =(hkm>1&hkm<12);
% 
% for id=id_all
% sig1= pz2(idh,id);
% dsig1= d_Pz2(idh,id);
% d_ln_pz2 = diff(log(pz2(:,id)),1,1);
% 
% figure
% subplot(1,2,1)
% plot(pz2(idh,id),hkm(idh),'b');hold on
% plot(pz2(cldBaseZ_ind(1,id),id),cldBaseZ(1,id),'o');
% plot(pz2(cldTopZ_ind(1,id),id),cldTopZ(1,id),'o');
% if ~isnan(cldBaseZ_ind(2,id))
% plot(pz2(cldBaseZ_ind(2,id),id),cldBaseZ(2,id),'^');
% plot(pz2(cldTopZ_ind(2,id),id),cldTopZ(2,id),'^');
% end
% set(gca,'XScale','log');legend('299nm Far','Cloud Base','Cloud Top');ylim([0,12]);
% xlabel('Range corrected signal (a.u.)');ylabel('Altitude (km)');
% title(['Ozone lidar P\cdotz^2 at ',datestr(DateTime_avg(id))])
% 
% subplot(1,2,2)
% plot(movmean(d_ln_pz2(idh(1:end-1)),40),hkm(idh),'b');hold on
% plot(d_ln_pz2(cldBaseZ_ind(1,id)),cldBaseZ(1,id),'o');
% plot(d_ln_pz2(cldTopZ_ind(1,id)),cldTopZ(1,id),'o');
% if ~isnan(cldBaseZ_ind(2,id))
% plot(d_ln_pz2(cldBaseZ_ind(2,id)),cldBaseZ(2,id),'^');
% plot(d_ln_pz2(cldTopZ_ind(2,id)),cldTopZ(2,id),'^');
% end
% set(gca,'XScale','linear');legend('299nm Far','Cloud Base','Cloud Top');ylim([0,12]);
% xlabel('derivative range corrected signal (a.u.)');ylabel('Altitude (km)');
% title(['Ozone lidar derivative P\cdotz^2 at ',datestr(DateTime_avg(id))])
% end

% for ii=1:len_t
%     if cldNum(ii) > 0
%         idcld=~isnan(cldBaseZ_ind(:,ii));
%         figure
%         subplot(1,3,1);
%         plot(ln_pz2(:,ii),hkm);hold on;
%         plot(ln_pz2(cldBaseZ_ind(idcld,ii),ii),cldBaseZ(idcld,ii),'o');
%         plot(ln_pz2(cldTopZ_ind(idcld,ii),ii),cldTopZ(idcld,ii),'^');
%         xlabel('ln(Pz^2)'); ylabel('altitude (km)'); grid on;
%         title(sprintf('ln(Pz^2), id = %.0d',ii)); ylim([hkm_cldsrc(1),hkm_cldsrc(end)]);
%         subplot(1,3,2);
%         plot(d_ln_pz2(:,ii),hkm);hold on;
%         plot(d_ln_pz2(cldBaseZ_ind(idcld,ii),ii),cldBaseZ(idcld,ii),'o');
%         plot(d_ln_pz2(cldTopZ_ind(idcld,ii),ii),cldTopZ(idcld,ii),'^');
%         xlabel('dln(Pz^2)/dz'); ylabel('altitude (km)'); grid on;ylim([hkm_cldsrc(1),hkm_cldsrc(end)]);
%         title(sprintf('dln(Pz^2)/dz, id = %.0d',ii));
%         subplot(1,3,3);
%         plot(d2_ln_pz2(:,ii),hkm);hold on;
%         plot(d2_ln_pz2(cldBaseZ_ind(idcld,ii),ii),cldBaseZ(idcld,ii),'o');
%         plot(d2_ln_pz2(cldTopZ_ind(idcld,ii),ii),cldTopZ(idcld,ii),'^');
%         xlabel('d^2ln(Pz^2)/dz^2'); ylabel('altitude (km)'); grid on;ylim([hkm_cldsrc(1),hkm_cldsrc(end)]);
%         title(sprintf('d^2ln(Pz^2)/dz^2, id = %.0d',ii));
%     end
% end