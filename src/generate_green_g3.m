function GreenG = generate_green_g3(XS,YS,ZS,LL,WW,DIP,STRIKE,xG,yG,zG,nu,fault_type,gps_type,dat_gps)
    % this function is used to generate greens function for gps data
    % fault_type(e.g. [1,1,0]) 1 for strike, 2 for dip, 3 for normal
    % gps_type(e.g. [1,1,0;1,0,1]) first row represent for horizontal, 
    %       second row represent for vertical
    GreenG = [];
    if isempty(xG)
        return
    end
    
    STRIKE = STRIKE/180*pi;
    DIP = DIP/180*pi;
    
    XC = XS + 0.5*LL.*sin(STRIKE);
    YC = YS + 0.5*LL.*cos(STRIKE);
    ZC = ZS;
    
    for j = 1:1:length(XC)
        for m = 1:1:3
            l = 1;
            for k = 1:1:length(dat_gps)
                G = [];
                x = xG(l:l-1+dat_gps(k)) - XC(j);
                y = yG(l:l-1+dat_gps(k)) - YC(j);
                z = zG(l:l-1+dat_gps(k));
            
                if fault_type(m) ~= 0
                    [ux,uy,uz] = calc_okada_in(1,x,y,z,nu,DIP(j),STRIKE(j),-ZC(j),LL(j),WW(j),m);
                    if gps_type(1,k) ~= 0
                        G = [G;ux;uy];
                    end
                    if gps_type(2,k) ~= 0
                        G = [G;uz];
                    end
                end
                l = l + ((gps_type(1,k) ~= 0)*2 + (gps_type(2,k) ~= 0))*dat_gps(k);
            end
            GreenG = [GreenG G];
        end
    end
        