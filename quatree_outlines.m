function [x,y,z,nn,X1,Y1,w,h,wn,hm] = quatree_outlines(a,b,c,r,up,down)
% a is x values
% b is y values
% c is the values to be subsampled
% r is the parameter that controls the sampling rate
% up & down are the thresholds below which the function will not subdivide
% a data patch.
% X1, Y1 are arrays of the coordinates of rectangles representing data points
% w, h are arrays of the widths/heights of rectangles
% wn, hm are arrays of subsampling steps

    x = [];
    y = [];
    z = [];
    nn = [];
    X1 = [];Y1=[];w=[];h=[];wn = [];hm = [];
    [m,n] = size(c);
    if sum(sum(isnan(c))) >= m*n*0.99 
        return %% patches with 99% NaN values would be abandoned
    end
   
    if m <= down || n <=down
        x = [x ; nanmean(reshape(a,numel(a),1))];
        y = [y ; nanmean(reshape(b,numel(b),1))];
        z = [z ; nanmean(reshape(c,numel(c),1))];
        nn = [nn ; sum(sum(~isnan(c)))]; % number of finite points (NaN) averaged in the smallest patch
        dx = (a(1,n)-a(1,1))/(n-1);
        dy = (b(m,1)-b(1,1))/(m-1);
        X1 = [X1;a(1,1)-dx/2];Y1 = [Y1;b(1,1)-dy/2]; 
        w = [w;(a(1,n)-a(1,1))*n/(n-1)];h=[h;(b(m,1)-b(1,1))*m/(m-1)];
        wn = [wn;n];hm = [hm;m];
        return
    end

    dx = (a(1,n)-a(1,1))/(n-1);
    dy = (b(m,1)-b(1,1))/(m-1);

    for j = 1:floor(m/2):floor(m/2)+1
        for k = 1:floor(n/2):floor(n/2)+1
            Step1 = floor(m/2); Step2 = floor(n/2);
            if (j ~= 1) && (mod(m,2) ~= 0) Step1 = ceil(m/2);end
            if (k ~= 1) && (mod(n,2) ~= 0) Step2 = ceil(n/2);end
            q = c(j:j+Step1-1,k:k+Step2-1) - nanmean(reshape(c(j:j+Step1-1,k:k+Step2-1),numel(c(j:j+Step1-1,k:k+Step2-1)),1));
            q(isnan(q)) = 0;
            if rms(rms(q))<= r && m <= up && n <= up  
                pz = c(j:j+Step1-1,k:k+Step2-1);
                px = a(j:j+Step1-1,k:k+Step2-1);
                py = b(j:j+Step1-1,k:k+Step2-1);
                if sum(sum(isnan(pz))) >= numel(pz)*0.75 % when 75% of a patch is NaN then the patch is abandoned
                    continue
                end
                px(isnan(pz)) = NaN;
                py(isnan(pz)) = NaN;
                x = [x ; nanmean(reshape(px,numel(px),1))];
                y = [y ; nanmean(reshape(py,numel(py),1))];
                z = [z ; nanmean(reshape(pz,numel(pz),1))];
                nn = [nn ; sum(sum(~isnan(pz)))];
                X1 = [X1;a(j,k)-dx/2];Y1 = [Y1;b(j,k)-dy/2]; 
                w = [w;(a(j,k+Step2-1)-a(j,k))*Step2/(Step2-1)];h=[h;(b(j+Step1-1,k)-b(j,k))*Step1/(Step1-1)];
                wn = [wn;Step1];hm = [hm;Step2];
            else
                [x0,y0,z0,n0,X10,Y10,w0,h0,wn0,hm0] = quatree_outlines(a(j:j+Step1-1,k:k+Step2-1),b(j:j+Step1-1,k:k+Step2-1),c(j:j+Step1-1,k:k+Step2-1),r,up,down); %when rms(q) > r, subdivide the patch
                x = [x ; x0];
                y = [y ; y0];
                z = [z ; z0];
                nn = [nn ; n0];
                X1 = [X1; X10]; Y1 = [Y1; Y10]; 
                w = [w;w0]; h = [h;h0]; wn=[wn;wn0]; hm = [hm;hm0];
            end
        end
    end
    
    return
