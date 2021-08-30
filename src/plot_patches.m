function [sp,dp] = plot_patches(str,mode,ratio)
    % this function is used to plot the result in 3D
    % mode 1 for strike slip, mode 2 for dip slip, mode 3 for normal slip
    a = load(str);
    ip = a(:,2);
    xs = a(:,4);
    ys = a(:,5);
    zs = a(:,6);

    L = a(:,7);
    W = a(:,8);
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
    
    %slip = log10(slip);
    
    figure('Position', [100, 100, 800, 600])
    hold on
    grid on
    for j = 1:length(slip)
        tx = [xs(j) , xs(j) + L(j)*sin(strike(j));
            xs(j) + W(j)*cos(dip(j))*cos(strike(j)), xs(j) + L(j)*sin(strike(j)) + W(j)*cos(dip(j))*cos(strike(j))];
        ty = [ys(j) , ys(j) + L(j)*cos(strike(j)) ;
            ys(j) - W(j)*cos(dip(j))*sin(strike(j)), ys(j) + L(j)*cos(strike(j)) - W(j)*cos(dip(j))*sin(strike(j))];
        tz = [zs(j) , zs(j);
            zs(j) - W(j)*sin(dip(j)), zs(j) - W(j)*sin(dip(j))];
        ts = repmat(slip(j),2,2);
        %if ip(j) == 7
        %    continue
        %end
        surf(tx,ty,tz,ts);
        x1 = a(j,11)*sin(strike(j))-a(j,12)*cos(dip(j))*cos(strike(j));
        x2 = a(j,11)*cos(strike(j))+a(j,12)*cos(dip(j))*sin(strike(j));
        x3 = a(j,12)*sin(dip(j));
        quiver3(mean(mean(tx))-1,mean(mean(ty))-1,mean(mean(tz)),(x1-1)*ratio,(x2-1)*ratio,x3*ratio,0,'k')
        quiver3(mean(mean(tx))+1,mean(mean(ty))+1,mean(mean(tz)),(x1+1)*ratio,(x2+1)*ratio,x3*ratio,0,'k')
    end

    hold off
    colormap(jet)
    view(-45,45)
    colorbar
    
    mu = 30e9;
    %m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    %m = 2/3*log10(mu*(sqrt(sum(a(:,11)/100.*L.*zs)^2+sum(a(:,12)/100.*L.*zs)^2)))-6.07;
    %m = 2/3*log10(mu*(sum(sqrt((a(:,11)/100.*L.*W).^2+(a(:,12)/100.*L.*W).^2))))-6.07;
    m = 2/3*log10(mu*(sqrt(sum(a(:,11)/100.*L.*W)^2+sum(a(:,12)/100.*L.*W)^2)))-6.07;

    disp(m)
    
    % plot slip versus depth
    figure
    hold on
    grid on
    sp = zeros(max(a(:,3)),1);
    dp = sp;
    for j = 1:1:max(a(:,3))
        ii = find(a(:,3) == j);
        if mode == 1 || mode == 13
            sp(j) = abs(a(ii,11))'*L(ii)/100/1000;      % strike slip
        else
            sp(j) = abs(a(ii,12))'*L(ii)/100/1000;      % dip slip
        end
        dp(j) = mean(zs(ii)-W(ii).*sin(dip(ii))/2)/1000;
    end
    plot(sp,dp,'b-s','LineWidth',2,'MarkerSize',10);
    
    