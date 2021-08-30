function plot_data(xx,yy,zz)
% used to plot the input phase and azimuth data
    figure
    numcol=200;   % number of colors in colormap
    cmap = jet(numcol);
    cmax = max(abs(zz));
    colormap(cmap);
    clima=[-cmax cmax];
    for i=1:length(zz)
        ind(i)=round(abs((zz(i)-clima(1))/(clima(2)-clima(1)))*(numcol-1))+1;
        line(xx(i),yy(i),'LineStyle','none','Marker','.', ...
        'MarkerSize',10,'MarkerFaceColor',cmap(ind(i),:), ...
        'MarkerEdgeColor',cmap(ind(i),:)), hold on
    end
    colorbar('Fontsize',18)
    set(gca,'box','on','CLim',clima,'Fontsize',18);