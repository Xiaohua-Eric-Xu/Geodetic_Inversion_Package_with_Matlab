function [xxx,yyy,ddd,azOO,tpp,sOO,dat_gov] = proj_gov2trace(patch,xall,yall,dall,azall,tall,sall,dat_gov,interval,closeness)
    % DESCRIPTION:
    % this function is used to project geological observations to the fault
    % trace we wanna use for inversion
    % the patch is a struct that contains {x,y,z,len,wid,dip,strike}, which
    % is the fault trace. Unit is degree and meter
    % x,y,d,az are column vectors
    % CREATED: 04/04/2014 BY: Xiaohua Xu
    
    % figure;scatter(x,y,20,d,'filled');colorbar;
    
    xxx = [];
    yyy = [];
    ddd = [];
    tpp = [];
    azOO = [];
    sOO = [];
    
    for i = 1:1:length(dat_gov)
        x = xall( sum(dat_gov(1:i))-dat_gov(i)+1 : sum(dat_gov(1:i)) );
        y = yall( sum(dat_gov(1:i))-dat_gov(i)+1 : sum(dat_gov(1:i)) );
        d = dall( sum(dat_gov(1:i))-dat_gov(i)+1 : sum(dat_gov(1:i)) );
        az = azall( sum(dat_gov(1:i))-dat_gov(i)+1 : sum(dat_gov(1:i)) );
        t = tall( sum(dat_gov(1:i))-dat_gov(i)+1 : sum(dat_gov(1:i)) );
        s = sall( sum(dat_gov(1:i))-dat_gov(i)+1 : sum(dat_gov(1:i)) );
        figure
        scatter(x,y,10,d,'filled');colorbar;
        
        xx = [];
        yy = [];
        ux = [];
        uy = [];
        dd = [];
        tp = [];
        ss = [];
        num = zeros(length(patch),1);
        dis = zeros(length(patch),1);
        num_num = [0;num];
    
        for j = 1:1:length(patch)
            num(j) = floor(patch(j).len/interval);
            dis(j) = patch(j).len/num(j);
        end
   
        for j = 1:1:length(patch)
            %if j == length(patch) || j == length(patch)-1 closeness = 2000; end
            k = -1/tand(90 - patch(j).strike);
            x0 = patch(j).x - 0.5*patch(j).len*sind(patch(j).strike);
            y0 = patch(j).y - 0.5*patch(j).len*cosd(patch(j).strike);
        
            for l = 1:1:num(j)
                x1 = x0 + dis(j)*sind(patch(j).strike);
                y1 = y0 + dis(j)*cosd(patch(j).strike);
            
                jdl = ((y >= k*(x-x0)+y0 & y <= k*(x-x1)+y1)|(y <= k*(x-x0)+y0 & y >= k*(x-x1)+y1)) & ((y <= tand(90 - patch(j).strike)*(x-x0)+y0 & y > tand(90 - patch(j).strike)*(x-x0)+y0-closeness/sind(patch(j).strike)) | (y >= tand(90 - patch(j).strike)*(x-x0)+y0 & y < tand(90 - patch(j).strike)*(x-x0)+y0-closeness/sind(patch(j).strike)));
                jdr = ((y >= k*(x-x0)+y0 & y <= k*(x-x1)+y1)|(y <= k*(x-x0)+y0 & y >= k*(x-x1)+y1)) & ((y > tand(90 - patch(j).strike)*(x-x0)+y0 & y < tand(90 - patch(j).strike)*(x-x0)+y0+closeness/sind(patch(j).strike)) | (y < tand(90 - patch(j).strike)*(x-x0)+y0 & y > tand(90 - patch(j).strike)*(x-x0)+y0+closeness/sind(patch(j).strike)));
                
                if sum(jdl) > 0 || sum(jdr) > 0
                    xtmp = (x'*jdl + x'*jdr)/(sum(jdl) + sum(jdr));
                    ytmp = (y'*jdl + y'*jdr)/(sum(jdl) + sum(jdr));
                    ttmp = (t'*jdl + t'*jdr)/(sum(jdl) + sum(jdr));
                    stmp = (s'*jdl + s'*jdr)/(sum(jdl) + sum(jdr));
                else
                    xtmp = x0;
                    ytmp = y0;
                    ttmp = 0;
                    stmp = 0;
                end
            
                if sum(jdl) > 0
                    dtl = sum(d.*cosd(az-patch(j).strike).*jdl)/sum(jdl);
                else
                    dtl = 0;
                end
                if sum(jdr) > 0
                    dtr = sum(d.*cosd(az-patch(j).strike).*jdr)/sum(jdr);
                else
                    dtr = 0;
                end

                xt = x0 + ([xtmp,ytmp]-[x0,y0])*[sind(patch(j).strike),cosd(patch(j).strike)]'*sind(patch(j).strike);
                yt = y0 + ([xtmp,ytmp]-[x0,y0])*[sind(patch(j).strike),cosd(patch(j).strike)]'*cosd(patch(j).strike);
                st = patch(j).strike;
                xx = [xx;xt];
                yy = [yy;yt];
                ux = [ux;(dtl+dtr)*sind(st)];
                uy = [uy;(dtl+dtr)*cosd(st)];
                dd = [dd;dtl+dtr];
                tp = [tp;ttmp];
                ss = [ss;stmp];
                x0 = x1;
                y0 = y1;
                
                % add some average
                dd(sum(num_num(1:j))+1:sum(num_num(1:j+1))) = conv(dd(sum(num_num(1:j))+1:sum(num_num(1:j+1))),ones(5,1)/5,'same');
            end
        end
    
    

        djd = (dd == 0);
        xx(djd) = [];
        yy(djd) = [];
        dd(djd) = [];
        ux(djd) = [];
        uy(djd) = [];
        tp(djd) = [];
        ss(djd) = [];
        azO = 90 - atan2(uy,ux)/pi*180;
        
        figure
        scatter(xx,yy,10,dd,'filled');colorbar;
    
        dat_gov(i) = length(xx);
        xxx = [xxx;xx];
        yyy = [yyy;yy];
        ddd = [ddd;dd];
        tpp = [tpp;tp];
        azOO = [azOO;azO];
        sOO = [sOO;ss];
    end
    
            
            
            
            
            
            
            
            