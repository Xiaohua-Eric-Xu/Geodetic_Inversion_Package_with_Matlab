function subsample(path,file,strain,uplimit,downlimit,type,write_or_not,plot_or_not)
% this function is used for subsampling data
% type: 1 for phase, other for azimuth offset

name = [path file '.grd'];

xvec = ncread(name,'x');
yvec = ncread(name,'y');
%xvec = ncread(name,'lon');
%yvec = ncread(name,'lat');
zz = ncread(name,'z');

if plot_or_not == 1
    figure
    %imagesc(flipud(zz')),colorbar
    imagesc(zz),colorbar
    colormap(jet)
    %[a,b] = size(zz);
end

zz = zz';
[xx,yy] = meshgrid(xvec,yvec);

% subsample the data
if type == 1
    disp('Quatree Sampling with mean...')
    [x,y,z,n] = quatree(xx,yy,zz,strain,uplimit,downlimit);
else
    disp('Quatree Sampling with median...')
    [x,y,z,n] = quatree_median(xx,yy,zz,strain,uplimit,downlimit);
end

n = 1./sqrt(n);

if plot_or_not == 1
    figure
    numcol=200;   % number of colors in colormap
    cmap = jet(numcol);
    colormap(cmap);
    clima=[-max([abs(z)]) max([abs(z)])];
    %clima=[min([z;zz]) max([z;zz])];

    for i=1:length(z)
        ind(i)=round(abs((z(i)-clima(1))/(clima(2)-clima(1)))*(numcol-1))+1;
            line(x(i),y(i),'LineStyle','none','Marker','.', ...
            'MarkerSize',12,'MarkerFaceColor',cmap(ind(i),:), ...
            'MarkerEdgeColor',cmap(ind(i),:)), hold on
    end
    hold off
    axis equal

    set(gca,'box','on','CLim',clima,'Fontsize',18);
    colorbar('Fontsize',18)
    xlabel('longitude')
    ylabel('latitude')
end

disp(['number of subsampled data: ' num2str(length(x))])
% write the result into a file
if write_or_not == 1
    fid = fopen([path '/' file '.llde'],'w');
    for j = 1:1:length(x)
        fprintf(fid,'%.9f\t%.9f\t%.9f\t%.9f\n',x(j),y(j),z(j),n(j));
    end
    fclose(fid);
end



