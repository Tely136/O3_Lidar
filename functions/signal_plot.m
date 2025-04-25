function signal_plot(title_str,legend_str,id,pc287,an287,pc299,an299,hkm)
    figure
    semilogx(pc287(:,id),hkm,an287(:,id),hkm,pc299(:,id),hkm,an299(:,id),hkm);

    title(title_str)
    legend(legend_str);
    % xlim([0,200]);
    ylim([0,8]);
    xlabel('Signal');
    ylabel('Altitude (km)');
end

