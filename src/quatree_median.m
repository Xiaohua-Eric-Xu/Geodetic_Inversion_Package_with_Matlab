function [x,y,z,nn] = quatree_median(a,b,c,r,up,down)
% a is x values
% b is y values
% c is the values to be subsampled
% r is the parameter that controls the sampling rate

    x = [];
    y = [];
    z = [];
    nn = [];
    
    [m,n] = size(c);
    
    if sum(sum(isnan(c))) >= m*n*0.99 
        return
    end
    
    if m <= down || n <=down
        x = [x ; nanmean(reshape(a,numel(a),1))];
        y = [y ; nanmean(reshape(b,numel(b),1))];
        z = [z ; nanmean(reshape(c,numel(c),1))];
        nn = [nn ; sum(sum(~isnan(c)))];
        return
    end
    
    for j = 1:floor(m/2):floor(m/2)+1
        for k = 1:floor(n/2):floor(n/2)+1
            %q = detrend2(a(j:j+floor(m/2)-1,k:k+floor(n/2)-1),b(j:j+floor(m/2)-1,k:k+floor(n/2)-1),c(j:j+floor(m/2)-1,k:k+floor(n/2)-1));
            %q(isnan(q)) = 0;
            %if rms(rms(q))<= r
            q = c(j:j+floor(m/2)-1,k:k+floor(n/2)-1) - nanmean(reshape(c(j:j+floor(m/2)-1,k:k+floor(n/2)-1),numel(c(j:j+floor(m/2)-1,k:k+floor(n/2)-1)),1));
            q(isnan(q)) = nanmedian(reshape(q,numel(q),1));
            if max(max(q))-min(min(q))<= r && m <= up && n <= up
            %if max(max(c(j:j+floor(m/2)-1,k:k+floor(n/2)-1)))-min(min(c(j:j+floor(m/2)-1,k:k+floor(n/2)-1))) <= r
                pz = c(j:j+floor(m/2)-1,k:k+floor(n/2)-1);
                px = a(j:j+floor(m/2)-1,k:k+floor(n/2)-1);
                py = b(j:j+floor(m/2)-1,k:k+floor(n/2)-1);
                if sum(sum(isnan(pz))) >= numel(pz)*0.5 
                    continue
                end
                
                px(isnan(pz)) = NaN;
                py(isnan(pz)) = NaN;
                pz = reshape(pz,numel(pz),1);
                ztmp = nanmedian(pz);
                [~,iztmp] = min(abs(pz-ztmp));
                px = reshape(px,numel(px),1);
                py = reshape(py,numel(py),1);
                xtmp = px(iztmp);
                ytmp = py(iztmp);
                
                x = [x ; xtmp];
                y = [y ; ytmp];
                z = [z ; ztmp];
                nn = [nn ; sum(~isnan(pz))];
            else
                [x0,y0,z0,n0] = quatree_median(a(j:j+floor(m/2)-1,k:k+floor(n/2)-1),b(j:j+floor(m/2)-1,k:k+floor(n/2)-1),c(j:j+floor(m/2)-1,k:k+floor(n/2)-1),r,up,down);
                x = [x ; x0];
                y = [y ; y0];
                z = [z ; z0];
                nn = [nn ; n0];
            end
        end
    end
    
    return
    
    