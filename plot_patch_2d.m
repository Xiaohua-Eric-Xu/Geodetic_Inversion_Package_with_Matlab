function plot_patch_2d(str,mode,ratio,sav)
    a = load(str);
    ip = a(:,2);
    il = a(:,3);
    xs = a(:,4);
    ys = a(:,5);
    zs = a(:,6);
    
    L = a(:,7)/1000;
    W = a(:,8)/1000;
    dip = a(:,9)/180*pi;
    strike = a(:,10)/180*pi;
    
    if mode == 1
        slip = a(:,11);
    elseif mode == 2
        slip = a(:,12);
    elseif mode == 3
        slip = a(:,13);
    elseif mode == 12
        slip = sqrt(a(:,11).^2+a(:,12).^2);
    elseif mode == 13
        slip = sqrt(a(:,11).^2+a(:,13).^2);
    elseif mode == 23
        slip = sqrt(a(:,12).^2+a(:,13).^2);
    else
        disp('ERROR: Please choose a correct mode for plotting')
        return
    end
    
    ncol = 200;
    clima = [min(slip),max(slip)];
    cmap = jet(ncol);
    
    for j = 1:1:length(unique(ip))
        h = figure;
        hold on
        ii = find(ip == j);
        dpth = sort(unique(W(ii)));
        dpth = cumsum(dpth);
        
        for k = 1:1:length(ii);
            if k ~= 1 && zs(ii(k-1)) ~= zs(ii(k))
                xx = -L(ii(k));
            elseif k == 1;
                xx = -L(ii(k));
            end
            
            ind=round(abs((slip(ii(k))-clima(1))/(clima(2)-clima(1)))*(ncol-1))+1;
            xx = xx + L(ii(k));
            yy = -dpth(il(ii(k)));
            rectangle('Position',[xx,yy,L(ii(k)),W(ii(k))],'FaceColor',cmap(ind,:),'EdgeColor','black');
            quiver(xx+0.5*L(ii(k)),yy+0.5*W(ii(k)),a(ii(k),11)*ratio/1000,a(ii(k),12)*ratio/1000,0,'Color','k','MaxHeadSize',0.8,'LineWidth',1.2,'SelectionHighlight','off');
        end
        axis equal,axis off
        if sav ~= 0
            saveas(h,['patch' num2str(j) '.eps']);
        end
    end
    colormap(cmap);        
            
        
        
        
        
        
        