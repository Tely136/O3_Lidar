function ND_plot(x,y,z,c_lim,y_lim,x_label,y_label,c_label,title_str,fs)
    figure;
    I = imagesc(x, y, z);
    cbar_tolnet=colorbar_tolnet_func;
    colormap(cbar_tolnet);
    ax = gca;
    ax.CLim = c_lim;
    set(gca,'YDir','normal');
    set(I,'AlphaData',~isnan(z));
    ylim(y_lim)
    c = colorbar;
    c.Label.String = c_label;
    xlabel(x_label)
    ylabel(y_label)
    title(title_str)
    fontsize(fs,"points")
end