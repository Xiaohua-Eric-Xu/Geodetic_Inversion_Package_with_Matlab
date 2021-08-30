function write_inv(str,xs,ys,zs,len,wid,dip,strike,Us,Ud,Un,num_grid)
    % this function is used to write the inversion results to a file named
    
    fid = fopen(str,'w');
    if Us == 0
        Us = zeros(size(xs));
    end
    if Ud == 0
        Ud = zeros(size(xs));
    end
    if Un == 0
        Un = zeros(size(xs));
    end
    
    Layer = 100;
    dpth(1) = Layer;
    for j = 2:1:length(xs)
        if zs(j) == zs(j-1);
            dpth(j) = dpth(j-1);
        elseif zs(j) > zs(j-1);
            dpth(j) = dpth(j-1)+1;
        else
            Layer = Layer + 100;
            dpth(j) = Layer;
        end
    end
    
    N = Layer/100;
    for j = 1:1:N
        ii = dpth >= j*100 & dpth < (j+1)*100;
        dpth(ii) = -(dpth(ii)-max(dpth(ii))-1);
    end
    
    l = 0;
    for j = 1:1:length(num_grid)
        for k = 1:1:num_grid(j)
            l = l+1;
            fprintf(fid,'%d\t%d\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n',...
                l,j,dpth(l),xs(l),ys(l),zs(l),len(l),wid(l),dip(l),strike(l),Us(l),Ud(l),Un(l));
        end
    end
        