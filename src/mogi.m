function green = mogi(x,y,xs,ys,dpth,dv,nu,look)
    % x and y are the location on the surface where deformation to be computed
    % dtph is the depth of the source
    % dv is the volume change
    % nu is possion's ration, usually we use 0.25
    % look is an vector on the look angle
    
    green = [];
    for i = 1:1:length(xs)
        R = sqrt((x-xs(i)).^2+(y-ys(i)).^2+(y*0+dpth(i)).^2);
        ratio = dv*(1-nu)/pi;
        ux = ratio*(x-xs(i))./R.^3;
        uy = ratio*(y-ys(i))./R.^3;
        uz = ratio*(y*0+dpth(i))./R.^3;
        %size(ux),size(look')
        green = [green,ux.*look(:,1)+uy.*look(:,2)+uz.*look(:,3)];
    end
    

