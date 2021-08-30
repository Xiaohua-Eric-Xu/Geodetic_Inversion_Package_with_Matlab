function [X,Y,Z,XB,YB,ZB,ll,ww] = generate_grid(mode,patchi,dw,dl,inc,plt)
    % this function is used to generate grids for each plane
    % mode 0 gives back the center of the top edge of the patch (X,Y,Z)
    % mode 1 gives back the edcmp source center(left up corner) of the patch (X,Y,Z)
    
    X0 = patchi.x;
    Y0 = patchi.y;
    Z0 = patchi.z; 
    len = patchi.len;
    wid = patchi.wid; 
    dip = patchi.dip;
    strike = patchi.strike;
    
    if mode ~= 0 && mode ~= 1
        disp('ERROR: Please choose the correct mode for grid_generating');
        return
    end

    dip = dip/180*pi;
    strike = strike/180*pi;
    if inc == 1 
        nw = round(wid/dw);                     % find the number of layers
    else
        nw = round(log(wid/dw*(inc-1)+1)/log(inc));
    end
    
    nl = round(len/dl);                         % top patch number
    
    if nw == 0, nw = 1; end
    if nl == 0, nl = 1; end
    
    if inc==1 
        dW(nw)=wid/nw;                          % top patch width 
    else
        dW(nw)=wid*(inc-1)/(inc^nw-1);
    end
    
    dL(nw)=len/nl;                              % top patch length
    nL(nw)=nl;                

    if nw>1 
        for j=1:nw-1
            dW(j)=dW(nw)*inc^(nw-j);            % patch width
        end
    end
    
    for j=nw-1:-1:1
        if round(len/dL(j+1)/inc) <= 1
            dL(j) = len;
            nL(j) = 1;
        else
            nL(j)=round(len/dL(j+1)/inc);       % patch length
            dL(j)=len/nL(j); 
        end
    end
    
    % generate the grid
    coss  = cos(strike);
    sins  = sin(strike);
    cosd  = cos(dip);
    sind  = sin(dip);
    
    XC = zeros(sum(nL),1);
    YC = zeros(sum(nL),1);
    ZC = zeros(sum(nL),1);
    ll = zeros(sum(nL),1);
    ww = zeros(sum(nL),1);
    XB = zeros(4*sum(nL),1);
    YB = zeros(4*sum(nL),1);
    ZB = zeros(4*sum(nL),1);
    
    for j = 1:nw
        if mod(nL(j),2) == 0  % even # of patches in length
            iC = nL(j)/2;
        else
            iC = fix(nL(j)/2);
        end
        depth0 = -Z0 + sum(dW(nw:-1:j+1))*abs(sind);  % depth0 is depth in vertical direction
        ex = sum(dW(nw:-1:j+1))*cosd;
        for i = 1:nL(j)
            k = sum(nL(1:j))-nL(j)+i;
            if mod(nL(j),2) == 0  % even # of patches in length
                ey = 0.5*(2*(i-iC)-1)*dL(j);
            else
                ey = (i-iC-1)*dL(j);
            end
            XC(k) = ex*coss + ey*sins+X0;
            YC(k) = -ex*sins + ey*coss+Y0;
            ZC(k) = -depth0;
            ll(k) = dL(j);
            ww(k) = dW(j);
            Lx = 0.5*dL(j)*sins;
            Ly = 0.5*dL(j)*coss;
            W = dW(j)*cosd;
            XB(4*(k-1)+1:4*k) = -[Lx-W*coss Lx -Lx -Lx-W*coss] + XC(k);
            YB(4*(k-1)+1:4*k) = -[Ly+W*sins Ly -Ly -Ly+W*sins] + YC(k);
            ZB(4*(k-1)+1:4*k) = [-dW(j)*sind 0 0 -dW(j)*sind] + ZC(k);
        end
    end
    
    if mode == 0
        X = XC; Y = YC; Z = ZC;
    elseif mode == 1
        X = XC - 0.5*ll.*sin(strike);
        Y = YC - 0.5*ll.*cos(strike);
        Z = ZC;
    end
    
    if plt == 1
        %figure
        plot3(X,Y,Z,'o'), hold on
        %pause
        plot3(X0,Y0,Z0,'kx'), hold on
        %pause
        for j=1:length(XC)
            line([XB((j-1)*4+1:(j-1)*4+4);XB((j-1)*4+1)],[YB((j-1)*4+1:(j-1)*4+4);YB((j-1)*4+1)],[ZB((j-1)*4+1:(j-1)*4+4);ZB((j-1)*4+1)],'Color','k'),hold on
            %pause(0.001)
        end
        grid on,axis equal
    end
    %pause
    
    