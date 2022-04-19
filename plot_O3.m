figure
plot((1064/287)^-4*t0prof(:,7),t0prof(:,1),'r');hold on;
legend('1064')
xlabel('Molecular ext (km^{-1})');
ylabel('Altitude (km)');
title(['Molecular ext profile on 14:26 (EST) Nov.18, 2021'])
ylim([0.05,5])
xlim([0,0.05])

figure
subplot(1,3,1)
plot(t0prof(:,2),t0prof(:,1),'r');hold on;
plot(t0prof(:,3),t0prof(:,1),'b');
legend('287nm','299nm')
xlabel('Totoal Backscatter (km^{-1}sr^{-1})');
ylabel('Altitude (km)');
title(['Totoal Backscatter profile on 14:26 (EST) Nov.18, 2021'])
ylim([0.05,5])
xlim([0,0.05])

subplot(1,3,2)
plot(log(t0prof(:,4)),t0prof(:,1),'g');
xlabel('Ln(\beta_{299}/\beta_{287}) (a.u.)');
ylabel('Altitude (km)');
title(['ln(\beta_{299}/\beta_{287}) profile at on 14:26 (EST) Nov.18, 2021'])
ylim([0.05,5])

subplot(1,3,3)
plot(log(t0prof(:,9)),t0prof(:,1),'g');
xlabel('$\frac{1}{2\delta\sigma_{O_3}}\cdot\frac{d ln(\beta_{299}/\beta_{287})}{dz} (/cm^3)$ ','Interpreter','latex');
ylabel('Altitude (km)');
title(['$\frac{1}{2\delta\sigma_{O_3}}\cdot\frac{d ln(\beta_{299}/\beta_{287})}{dz} (/cm^3)$  profile on 14:26 (EST) Nov.18, 2021'],'Interpreter','latex')
ylim([0.05,5])