function GreenO = generate_green_o(XS,YS,ZS,LL,WW,DIP,STRIKE,xO,yO,azO,tpO,nu,fault_type)
    % this function is used to generate greens function for geological data
    % fault_type(e.g. [1,1,0]) 1 for strike, 2 for dip, 3 for normal
    % fault azimuth is consistent with the slip on its left side slip
    % direction. for az:-->, the upper part right is its positive direction
    
    
    GreenO = [];
    if isempty(xO)
        return
    end
    
    STRIKE = STRIKE/180*pi;
    DIP = DIP/180*pi;
    
    XC = XS + 0.5*LL.*sin(STRIKE);
    YC = YS + 0.5*LL.*cos(STRIKE);
    ZC = ZS;
    
    for j = 1:1:length(XC)
        x = xO - XC(j);
        y = yO - YC(j);
    	dx = cosd(180-azO)*10;
        dy = sind(180-azO)*10;
        for k = 1:1:3
            if fault_type(k) ~= 0
                [ux1,uy1,uz1] = calc_okada(1,1,x+dx,y+dy,nu,DIP(j),-ZC(j),LL(j),WW(j),k,STRIKE(j),tpO);
                [ux2,uy2,uz2] = calc_okada(1,1,x-dx,y-dy,nu,DIP(j),-ZC(j),LL(j),WW(j),k,STRIKE(j),tpO);
                GreenO = [GreenO (ux1-ux2).*cosd(90-azO)+(uy1-uy2).*sind(90-azO)];
            end
        end
    end
    %figure;plot(x,y,'ko');hold on;plot(x+dx,y+dy,'rx');plot(x-dx,y-dy,'bx');axis equal