function patch = trace2patch(lon,lat,xo,yo)
    % used to generate patch from trace
    %[x0,y0] = utm2ll(xo,yo,0,1);
    %[x,y] = utm2ll(lon,lat,0,1);
    [x0,y0]=utm2ll(xo,yo,0,1);
    [x,y] = utm2ll(lon,lat,0,1);
    x = x - x0; 
    y = y - y0;
    
    patch(1:length(lon)/2) = struct('x',0,'y',0,'z',0,'len',0,'wid',20e3,'dip',90,'strike',0);
    
    for j = 1:1:length(x)/2
        patch(j).x = (x((j-1)*2+1)+x(j*2))/2;
        patch(j).y = (y((j-1)*2+1)+y(j*2))/2;
        patch(j).z = 0;
        patch(j).len = sqrt((x((j-1)*2+1)-x(j*2))^2+(y((j-1)*2+1)-y(j*2))^2);
        patch(j).strike = 90 - atan2(y(j*2)-y((j-1)*2+1),x(j*2)-x((j-1)*2+1))/pi*180;
        fprintf('        {trace%d}\n',j)
        fprintf('         x = %e\n',patch(j).x);
        fprintf('         y = %e\n',patch(j).y);
        fprintf('         z = %e\n',patch(j).z);
        fprintf('       len = %e\n',patch(j).len);
        fprintf('       wid = %e\n',patch(j).wid);
        fprintf('       dip = %e\n',patch(j).dip);
        fprintf('    strike = %e\n',patch(j).strike);
        fprintf('\n')

    end
    plot_trace(patch)
    hold on, plot(0,0,'b*')
    axis equal
    grid on
    hold off
        
    
    
    
    