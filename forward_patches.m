function [ux,uy,uz] = forward_patches(str,x,y)
    % this function is used to plot the result on ground surface with
    % forwarding the model
    % mode 1 for x, mode 2 for y, mode 3 for z
    % there is a bug !!!!!!!
    % 
    a = load(str);
    xs = a(:,4);
    ys = a(:,5);
    zs = a(:,6);

    L = a(:,7);
    W = a(:,8);
    dip = a(:,9)/180*pi;
    strike = a(:,10)/180*pi;
    Us = a(:,11);
    Ud = a(:,12);
    Un = a(:,13);
    clear a
    
    %xvec = [xmin:xspac:xmax];
    %yvec = [ymin:yspac:ymax];
    
    ux = zeros(length(x),1);
    %ux = reshape(ux,numel(ux),1);
    size(ux)
    uy = ux;
    uz = ux;
    tp = ux;

    str = 'working on patch   1 ...';
    fprintf('\n');
    for j = 1:length(Us)
        fprintf(str);
        %[x,y] = meshgrid(xvec-ys(j),yvec-xs(j));
        %x = reshape(x,numel(x),1);
        %y = reshape(y,numel(y),1);
        if Us(j) ~= 0
            [ux1,uy1,uz1] = ...
                calc_okada(1,Us(j),x-ys(j),y-xs(j),0.25,dip(j),zs(j),L(j),W(j),1,strike(j),tp);
        end
    	if Ud(j) ~= 0
            [ux2,uy2,uz2] = ...
                calc_okada(1,Ud(j),x-ys(j),y-xs(j),0.25,dip(j),zs(j),L(j),W(j),2,strike(j),tp);
        end
        if Un(j) ~= 0
            [ux3,uy3,uz3] = ...
                calc_okada(1,Un(j),x-ys(j),y-xs(j),0.25,dip(j),zs(j),L(j),W(j),3,strike(j),tp);
        end
        ux = ux + ux1 + ux2;
        uy = uy + uy1 + uy2;
        uz = uz + uz1 + uz2;
        tmp_str = repmat('\b',1,4+numel(num2str(j)));
        str = [tmp_str,num2str(j),' ...'];
    end
    fprintf('\n');
    fprintf('Finished computing ...\n');
    %ux = reshape(ux,length(xvec),length(yvec));
    %uy = reshape(uy,length(xvec),length(yvec));
    %uz = reshape(uz,length(xvec),length(yvec));
    %figure,imagesc(flipud(ux')),colormap(jet),colorbar,caxis([-100 100]), title('ux')
    %figure,imagesc(flipud(uy')),colormap(jet),colorbar,caxis([-100 100]), title('uy')
    %figure,imagesc(flipud(uz')),colormap(jet),colorbar,caxis([-100 100]), title('uz')

    %if mode == 1
    %    imagesc(flipud(ux'))
    %elseif mode == 2
    %    imagesc(flipud(uy'))
    %elseif mode == 3
    %    imagesc(flipud(uz'))
    %end
    
    
    
    