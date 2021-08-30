function plot_trace(patch)
    % this function is used to plot trace on the map
    N = length(patch)
    
    for j = 1:1:N
        x(j) = patch(j).x;
        y(j) = patch(j).y;
        len(j) = patch(j).len;
        strike(j) = patch(j).strike;
    end
    
    hold on
    for j = 1:1:N
        plot(x,y,'bo');
        x1 = x(j) - 0.5*len(j).*sind(strike(j));
        y1 = y(j) - 0.5*len(j).*cosd(strike(j));
        x2 = x(j) + 0.5*len(j).*sind(strike(j));
        y2 = y(j) + 0.5*len(j).*cosd(strike(j));
        plot(x1,y1,'b.');
        plot(x2,y2,'b.');
        plot([x1,x2],[y1,y2],'k')
    end