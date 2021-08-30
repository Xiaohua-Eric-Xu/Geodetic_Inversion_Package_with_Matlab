function v = patch_smooth(p,str,Nl,Nw,mode,plot_or_not)
    % For converting a okada patches solution to a spectral solution
    % p is patch paramaters that has xc,yc,zc,width,length,dip,strike
    % xs, ys and zs are the divided patch source
    % u is the corresponding slip
    % ll, ww are the corresponding length and width for the patch source
    % Nl Nw are the num of divisions of length and width for p
    
    a = load(str);
    xs = a(:,4);
    ys = a(:,5);
    zs = a(:,6);
    ll = a(:,7);
    ww = a(:,8);
    idpth = a(:,3);
    maxdpth = max(idpth);
    
    if mode == 1
        u = a(:,11);
    elseif mode == 2
        u = a(:,12);
    elseif mode == 3
        u = a(:,13);
    elseif mode == 12
        u = sqrt(a(:,11).^2+a(:,12).^2);
    elseif mode == 13
        u = sqrt(a(:,11).^2+a(:,13).^2);
    elseif mode == 23
        u = sqrt(a(:,12).^2+a(:,13).^2);
    else
        disp('ERROR: Please choose a correct mode for plotting')
        return
    end
    
    x = zeros(Nl,Nw);
    y = x;
    z = x;
    
    xc = xs + 0.5*ll*sind(p.strike) + 0.5*ww*cosd(p.dip)*cosd(p.strike);
    yc = ys + 0.5*ll*cosd(p.strike) - 0.5*ww*cosd(p.dip)*sind(p.strike);
    zc = zs - 0.5*ww*sind(p.dip);
    
    x0 = p.x - 0.5*p.len*sind(p.strike);
    y0 = p.y - 0.5*p.len*cosd(p.strike);
    z0 = p.z;
    
    ext_x = [];
    ext_y = [];
    ext_z = [];
    ext_u = [];
    ext = 1;
    for j = 1:1:maxdpth
        ii = find(idpth==j);
        jj = min(ii);
        ext_x = [ext_x; xs(jj) - ext*ll(jj)*sind(p.strike)];
        ext_y = [ext_y; ys(jj) - ext*ll(jj)*cosd(p.strike)];
        ext_z = [ext_z; zs(jj)];
        ext_u = [ext_u; u(jj)];
        if j == 1
            ext_x = [ext_x; xs(jj) - ext*ll(jj)*sind(p.strike) - ext*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) - ext*ll(jj)*cosd(p.strike) + ext*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) + ext*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
            jj = ii;
            ext_x = [ext_x; xs(jj) - ext*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) + ext*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) + ext*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
        end
        if j == maxdpth
            ext_x = [ext_x; xs(jj) - ext*ll(jj)*sind(p.strike) + (1+ext)*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) - ext*ll(jj)*cosd(p.strike) - (1+ext)*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) - (1+ext)*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
            jj = ii;
            ext_x = [ext_x; xs(jj) + (1+ext)*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) - (1+ext)*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) - (1+ext)*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
        end
        
        jj = max(ii);
        ext_x = [ext_x; xs(jj) + (1+ext)*ll(jj)*sind(p.strike)];
        ext_y = [ext_y; ys(jj) + (1+ext)*ll(jj)*cosd(p.strike)];
        ext_z = [ext_z; zs(jj)];
        ext_u = [ext_u; u(jj)];
        if j == 1
            ext_x = [ext_x; xs(jj) + (1+ext)*ll(jj)*sind(p.strike) - ext*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) + (1+ext)*ll(jj)*cosd(p.strike) + ext*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) + ext*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
        end
        if j == maxdpth
            ext_x = [ext_x; xs(jj) + ll(jj)*sind(p.strike) + (1+ext)*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) + ll(jj)*cosd(p.strike) - (1+ext)*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) - (1+ext)*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
            ext_x = [ext_x; xs(jj) + (1+ext)*ll(jj)*sind(p.strike) + (1+ext)*ww(jj)*cosd(p.dip)*cosd(p.strike)];
            ext_y = [ext_y; ys(jj) + (1+ext)*ll(jj)*cosd(p.strike) - (1+ext)*ww(jj)*cosd(p.dip)*sind(p.strike)];
            ext_z = [ext_z; zs(jj) - (1+ext)*ww(jj)*sind(p.dip)];
            ext_u = [ext_u; u(jj)];
        end 
    end
    
    dl = p.len/Nl;
    dw = p.wid/Nw;
    
    stepxl = dl*sind(p.strike);
    stepxw = dw*cosd(p.dip)*cosd(p.strike);
    stepyl = dl*cosd(p.strike);
    stepyw = -dw*cosd(p.dip)*sind(p.strike);
    stepz = -dw*sind(p.dip);
    
    x0 = x0 + 0.5*(stepxl+stepxw);
    y0 = y0 + 0.5*(stepyl+stepyw);
    z0 = z0 + 0.5*stepz;
    
    [xx,yy]=meshgrid(1:1:Nl,1:1:Nw);
    
    x = x0 + (xx-1)*stepxl + (yy-1)*stepxw;
    y = y0 + (xx-1)*stepyl + (yy-1)*stepyw;
    z = z0 + (yy-1)*(stepz);
    
    xc = [xc;ext_x];
    yc = [yc;ext_y];
    zc = [zc;ext_z];
    u = [u;ext_u];

    %figure,surf(x,y,z,'EdgeColor','None'),hold on, plot3(xc,yc,zc,'x'),grid on,plot3(ext_x,ext_y,ext_z,'r*')
    %figure,plot3(xc,yc,zc,'o'),grid on
    %x = reshape(x,Nl*Nw,1);
    %y = reshape(y,Nl*Nw,1);
    %z = reshape(z,Nl*Nw,1);
    
    
    v = griddata(xc,yc,zc,u,x,y,z);
    
    %x = reshape(x,Nl,Nw);
    %y = reshape(y,Nl,Nw);
    %z = reshape(z,Nl,Nw);
    %v = reshape(v,Nl,Nw);
    if plot_or_not ~= 0
        figure
        surf(x,y,z,v,'EdgeColor','None');
        colormap(jet),colorbar
        %hold on
        %plot3(xc,yc,zc,'x'),grid on,plot3(ext_x,ext_y,ext_z,'r*')
    end
    
   