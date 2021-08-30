function GreenA = generate_green_a(XS,YS,ZS,LL,WW,DIP,STRIKE,xA,yA,tpA,phi,nu,fault_type,dat_az)
    % this function is used for generating the Green's function for
    % azimuth_offset data
    GreenA = [];
    if isempty(xA)
        return
    end
    
    PHI = [];
    phi = phi/180*pi;
    
    for j = 1:1:length(phi)
        PHI = [PHI; zeros(dat_az(j),1)+phi(j)];
    end
     
    STRIKE = STRIKE/180*pi;
    DIP = DIP/180*pi;
    
    XC = XS + 0.5*LL.*sin(STRIKE);
    YC = YS + 0.5*LL.*cos(STRIKE);
    ZC = ZS;
    
    for j = 1:length(XC)
        x = xA - XC(j);
        y = yA - YC(j);
    	
        for k = 1:1:3
            if fault_type(k) ~= 0
                [ux,uy,uz] = calc_okada(1,1,x,y,nu,DIP(j),-ZC(j),LL(j),WW(j),k,STRIKE(j),tpA);
                GreenA = [GreenA ux.*sin(PHI) + uy.*cos(PHI)];
            end
        end
    end