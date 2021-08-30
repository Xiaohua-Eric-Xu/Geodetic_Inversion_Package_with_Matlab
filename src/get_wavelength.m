function a = get_wavelength(img,dl,dw,plot_or_not)
    [x,y] =  find_local_minmax(img,10,10,plot_or_not);
    x = [x;0;0;size(img,1);size(img,1);0;floor(size(img,1)/2);size(img,1);floor(size(img,1)/2)];
    y = [y;0;size(img,2);0;size(img,2);floor(size(img,2)/2);0;floor(size(img,2)/2);size(img,2)];
    tri = delaunay(x,y);
    z = zeros(size(x));
    iz = z;
    for ii = 1:1:size(tri,1)
        i1 = tri(ii,1);
        i2 = tri(ii,2);
        i3 = tri(ii,3);
        
        z(i1) = z(i1)+sqrt((x(i1)-x(i2))^2*dw^2+(y(i1)-y(i2))^2*dl^2)+sqrt((x(i1)-x(i3))^2*dw^2+(y(i1)-y(i3))^2*dl^2);
        z(i2) = z(i2)+sqrt((x(i2)-x(i1))^2*dw^2+(y(i2)-y(i1))^2*dl^2)+sqrt((x(i2)-x(i3))^2*dw^2+(y(i2)-y(i3))^2*dl^2);
        z(i3) = z(i3)+sqrt((x(i3)-x(i2))^2*dw^2+(y(i3)-y(i2))^2*dl^2)+sqrt((x(i3)-x(i1))^2*dw^2+(y(i3)-y(i1))^2*dl^2);
        
        iz(i1) = iz(i1)+2;
        iz(i2) = iz(i2)+2;
        iz(i3) = iz(i3)+2;
    end
    
    for ii = 1:1:length(z)
        z(ii) = z(ii)/iz(ii);
    end
    z = z/1000;
    if plot_or_not ~= 0
        figure, trisurf(tri,x,y,z),colorbar
    end
    
    
    fid = fopen('tmp.xyz','w');
    for i = 1:1:length(x)
        fprintf(fid,'%f\t%f\t%f\n',x(i),y(i),z(i));
    end
    
    
    fclose(fid);
    str = ['!gmt surface tmp.xyz -R0/',num2str(size(img,1)),'/0/',num2str(size(img,2)),' -I1/1 -T0.1 -r -Gtmp.grd'];
    eval(str);
    a = ncread('tmp.grd','z')';
    a = imgaussfilt(a,21);
    %ft = fspecial('gaussian',101,21);
    %a = conv2(a,ft,'same');
    %w = ones(size(z));
    %a = xyz2surface(x,y,z,w,[1:1:size(img,1)],[1:1:size(img,2)],0.1,'smoothness',0.0001,'regularizer','laplacian');
    if plot_or_not ~= 0
        figure,imagesc(a);axis xy;colormap(jet),colorbar,caxis([0 max(max(a))])
    end
    