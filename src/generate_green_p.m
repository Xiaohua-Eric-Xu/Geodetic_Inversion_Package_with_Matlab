function GreenP = generate_green_p(XS,YS,ZS,LL,WW,DIP,STRIKE,xP,yP,tpP,look,nu,fault_type)
    % used to generate green's function matrix for phase data
    % fault_type(e.g. [1,1,0]) 1 for strike, 2 for dip, 3 for normal
    
    GreenP = [];
    if isempty(xP)
        return
    end
    
    STRIKE = STRIKE/180*pi;
    DIP = DIP/180*pi;
    
    XC = XS + 0.5*LL.*sin(STRIKE);
    YC = YS + 0.5*LL.*cos(STRIKE);
    ZC = ZS;
    
    
    for j = 1:length(XC)
        x = xP - XC(j);
        y = yP - YC(j);
    	
        for k = 1:1:3
            if fault_type(k) ~= 0
                [ux,uy,uz] = calc_okada(1,1,x,y,nu,DIP(j),-ZC(j),LL(j),WW(j),k,STRIKE(j),tpP);
                GreenP = [GreenP ux.*look(:,1) + uy.*look(:,2) + uz.*look(:,3)];
            end
        end
    end