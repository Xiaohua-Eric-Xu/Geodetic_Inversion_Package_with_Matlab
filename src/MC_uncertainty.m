close all
NN = 1000;amp = 1;mp_size = 2000;
U_test = zeros(size(U,1),NN);
xxx = -mp_size*100:100:mp_size*100;
yyy = -mp_size*100:100:mp_size*100;

for k = 1:1:NN
    dP_tmp = [];dA_tmp = [];
    str = ['Loading data...Realization ',num2str(k),'...'];
    disp(str)
    % load the descending data
    if num_des ~= 0
        for j = 1:1:num_des
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['des' num2str(j)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            xx = xx - xo;
            yy = yy - yo;
            noise = amp*rednoise(mp_size*2+1,mp_size*2+1,1.33);
            for n = 1:1:length(xx)
                dP_tmp = [dP_tmp;dat(n,7)+noise(round(xx(n)/1000+mp_size),round(yy(n)/1000)+mp_size)];             % displacement (in cm)  (in case of different difinition, this value might need to be reversed.)
            end
        end
    end
    
    if num_asc ~= 0
        for j = 1:1:num_asc
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['asc' num2str(j)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            xx = xx - xo;
            yy = yy - yo;
            noise = amp*rednoise(mp_size*2+1,mp_size*2+1,1.33);
            for n = 1:1:length(xx)
                dP_tmp = [dP_tmp;dat(n,7)+noise(round(xx(n)/1000+mp_size),round(yy(n)/1000)+mp_size)];             % displacement (in cm)  (in case of different difinition, this value might need to be reversed.)
            end             
        end
    end
    
    if num_azi ~= 0
        for j = 1:1:num_azi
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['azi' num2str(j)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            xx = xx - xo;
            yy = yy - yo;
            noise = amp*rednoise(mp_size*2+1,mp_size*2+1,1.33);
            for n = 1:1:length(xx)
                dA_tmp = [dA_tmp;dat(n,3)*100+noise(round(xx(n)/1000+mp_size),round(yy(n)/1000)+mp_size)];             % displacement (in cm)  (in case of different difinition, this value might need to be reversed.)
            end   
        end
    end
    
    % GPS
    dG_tmp = dG + sG.*randn(size(sG))/length(sG)*GW;
    
    % load the geological observation data
    dOO_tmp = dOO + sOO.*randn(size(sOO))*mean(sO)/mean(sOO);
    
    D_new = [SWP*dP_tmp;dA_tmp;dG_tmp;dOO_tmp];
    D_all_new = [W.*D_new;zeros(size(Green_all,1)-size(D_new,1),1)];
    
    rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,100);
    Green_all_new = [Greens rmp; Smooth/mean(max(Smooth,[],2))*SF/size(Smooth,1) zeros(size(Smooth,1),size(rmp,2))];
    for n = 1:1:length(W)
        Green_all_new(n,:) = Green_all_new(n,:)*W(n);
    end
    
    [U,~,~,~] = lsqlin(Green_all_new,D_all_new,[],[],[],[],lb,ub);
    
    rmp0 = rmp;
    Urmp0 = U(end-size(rmp0,2)+1:end);
    
    D_all_new = [W.*(D_new-rmp0*Urmp0); zeros(size(Green_all_new,1)-size(D_new,1),1)];
    rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,1e6);
    Green_all_new = [Greens rmp; Smooth/mean(max(Smooth,[],2))*SF/size(Smooth,1) zeros(size(Smooth,1),size(rmp,2))];
    for j = 1:1:length(W)
        Green_all_new(j,:) = Green_all_new(j,:)*W(j);
    end
    
    [U,~,~,~] = lsqlin(Green_all_new,D_all_new,[],[],[],[],lb,ub);
    
    U_test(:,k) = U;
end

u_std = std(U_test'-repmat(mean(U_test'),NN,1));
Us_std = u_std(1:2:end-size(rmp,2));
Ud_std = u_std(2:2:end-size(rmp,2));
write_inv('std.inv',XS,YS,ZS,LL,WW,DIP,STRIKE,Us_std,Ud_std,0,num_grid);
plot_patches('std.inv',13,0)


