function clima = subsample_outlines(path,file,strain,uplimit,downlimit,type,write_or_not,plot_or_not,Axes,clima2)
% this function is used for subsampling data and ploting as rectangles
% type: 1 for phase, other for azimuth offset
% Axes: used for subplot controls
% clima2: used for homogenizing colorbars for multiple subplots

name = [path file '.grd'];

xvec = ncread(name,'x');yvec = ncread(name,'y');
%cxvec = ncread(name,'lon');yvec = ncread(name,'lat');
zz = ncread(name,'z');

%if plot_or_not == 1
%    figure('WindowStyle','docked');axis equal;
%    %imagesc(flipud(zz')),colorbar
%    imagesc(zz);colorbar;colormap(jet);axis equal;
%    %[a,b] = size(zz);
%end

zz = zz';[xx,yy] = meshgrid(xvec,yvec);
[xmin,xmax,ymin,ymax] = findindex(zz); %find the index for finite(NaN) data
zz = zz(xmin:xmax,ymin:ymax);xx = xx(xmin:xmax,ymin:ymax);yy = yy(xmin:xmax,ymin:ymax);
%figure('WindowStyle','docked');
axis equal;hold on;
% subsample the data
if type == 1
    disp('Quatree Sampling with mean...')
    [x,y,z,n,X1,Y1,w,h,wn,hm] = quatree_outlines(xx,yy,zz,strain,uplimit,downlimit);
else
    disp('Quatree Sampling with median...')
    [x,y,z,n] = quatree_median(xx,yy,zz,strain,uplimit,downlimit);
end
%Area = w.*h;
%Density = n./Area;
n = 1./sqrt(n);

if plot_or_not == 1
    %figure('WindowStyle','docked');axis equal;
    numcol=200;  cmap = jet(numcol); colormap(cmap);% number of colors in colormap
    %cmap2 = gray(numcol);
    
    clima=max(abs(z));
    if clima2 > clima         clima = clima2;end

    %ind = ones(length(z),1);ind2 = ind;
    %clima=[min([z;zz]) max([z;zz])];
    ind = round((z+clima)/2/clima*199)+1;
    for i=1:length(z)
        %ind(i)=round(abs((z(i)-clima(1))/(clima(2)-clima(1)))*(numcol-1))+1;
        %ind2(i)=round((Density(i)-min(Density))/(max(Density)-min(Density))*(numcol-1))+1;
        %rectangle('Position',[X1(i),Y1(i),w(i),h(i)]);  hold on;
        rectangle('Position',[X1(i)-360,Y1(i),w(i),h(i)],'FaceColor',cmap(ind(i),:));  hold on;
        %line(x(i),y(i),'LineStyle','none','Marker','.','MarkerSize',12,'MarkerFaceColor',cmap(ind(i),:),'MarkerEdgeColor',cmap(ind(i),:)); 
    end
    
    set(Axes,'box','on','CLim',[-clima clima],'Fontsize',18);
    colorbar('Fontsize',18)
    xlabel('Longitude')
    ylabel('Latitude')
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

