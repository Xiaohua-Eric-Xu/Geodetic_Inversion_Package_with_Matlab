function [x,y] = find_local_minmax(img,nx,ny,plot_or_not)

    x = [];y = [];
    a = FastPeakFind(img);
    b = FastPeakFind(-img);
    for i = 1:2:length(a)
        idx = a(i)-nx:1:a(i)+nx;
        idy = a(i+1)-ny:1:a(i+1)+ny;
        idx(idx<=0 | idx>size(img,1)) = [];
        idy(idy<=0 | idy>size(img,2)) = [];
        if img(a(i+1),a(i)) >= max(max(img(idy,idx)))*0.999
            x = [x;a(i)];
            y = [y;a(i+1)];
        %else
            %str = ['img: ',num2str(img(a(i+1),a(i))),'  max: ',num2str(max(max(img(idy,idx))))];
            %disp(str)
        end
    end
  
    for i = 1:2:length(b)
        idx = b(i)-nx:1:b(i)+nx;
        idy = b(i+1)-ny:1:b(i+1)+ny;
        idx(idx<=0 | idx>size(img,1)) = [];
        idy(idy<=0 | idy>size(img,2)) = [];
        if -img(b(i+1),b(i)) >= max(max(-img(idy,idx)))*0.999
            x = [x;b(i)];
            y = [y;b(i+1)];
        %else
        	%str = ['img: ',num2str(-img(b(i+1),b(i))),'  max: ',num2str(max(max(-img(idy,idx))))];
            %disp(str)
        end
    end

    ii = 1;
    while ii<=length(x)
        kk = [];
        for jj = 1:1:length(x)
            if ii ~= jj && abs(x(ii)-x(jj)) <= nx && abs(y(ii)-y(jj)) <= ny
                kk = [kk,jj];
            end
        end
        if isempty(kk) 
            ii = ii+1;
        else
            x(kk) = [];
            y(kk) = [];
        end
    end
    if plot_or_not ~= 0
        figure, imagesc(img),axis xy, hold on, plot(a(1:2:end),a(2:2:end),'mo'),plot(x,y,'ko')
    end
            





%{
    x = [];y = [];
    for ii = 1:1:size(img,1)
        for jj = 1:1:size(img,2)
            idx = ii-nx:1:ii+nx;
            idy = jj-ny:1:jj+ny;
            idx(idx<=0 | idx>size(img,1)) = [];
            idy(idy<=0 | idy>size(img,2)) = [];
            if img(ii,jj) >= max(max(img(idx,idy)))
                x = [x;ii];
                y = [y;jj];
            end
        end
    end
%}