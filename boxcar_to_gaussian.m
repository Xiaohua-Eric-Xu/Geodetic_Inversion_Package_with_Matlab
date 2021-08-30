function gxy = boxcar_to_gaussian(xwid,ywid)
    % convert boxcar filter to semi-equivalent gaussian filter
    %clear all
    %xwid = 19;
    %ywid = 7;
    r = 0.3;
    wid_ext = 0.5;
    
    x = [zeros(round(ywid*wid_ext),1);ones(xwid,1)/xwid;zeros(round(xwid*wid_ext),1)];
    fx = abs(fftshift(fft(x)));
    figure, plot(fx,'k--.'),hold on
    gx = sum(fspecial('gaussian',length(x),xwid*r));
    fgx = abs(fftshift(fft(gx)));
    plot(fgx,'r--x'),grid on
    
    y = [zeros(round(ywid*wid_ext),1);ones(ywid,1)/ywid;zeros(round(ywid*wid_ext),1)];
    gxy = zeros(length(x),length(y));
    for i = 1:1:length(x)
        for j = 1:1:length(y)
            gxy(i,j) = exp(-((i-(floor(length(x)/2)+1))^2/2/(xwid*r)^2+(j-(floor(length(y)/2)+1))^2/2/(ywid*r)^2));
        end
    end
    %gxy = fspecial('gaussian',[length(x) length(y)], [xwid*r ywid*r]);
    %figure,imagesc(gxy),colormap(jet),colorbar
    
    
    
    
