function ND_plot_ax(ax,x,y,z,c_lim,y_lim,title_str)
    I = imagesc(ax, x, y, z);
    cbar_tolnet = colorbar_tolnet_func;
    colormap(ax,cbar_tolnet);
    % % ax = gca;
    ax.CLim = c_lim;
    set(ax,'YDir','normal');
    set(I,'AlphaData',~isnan(z));
    ylim(ax,y_lim)
    colorbar(ax)
    title(ax,title_str)
end