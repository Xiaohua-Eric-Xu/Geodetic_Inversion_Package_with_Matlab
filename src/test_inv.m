% this is used to testify whether the inversion code could recover isolated
% slip patches.
% run main_inv.m first and then run this piece of code.

U_test = zeros(size(U));
for j = 1:1:(length(U)-length(Urmp))/2
    if rand()<0.2
        U_test(2*j-1) = -300;
    end
end
Us_test = U_test(1:2:end-size(rmp,2));
Ud_test = U_test(2:2:end-size(rmp,2));
write_inv('test1.inv',XS,YS,ZS,LL,WW,DIP,STRIKE,Us_test,Ud_test,0,num_grid);
plot_patches('test1.inv',13,6)
%%
D_t = Greens*U_test(1:end-length(Urmp));
Green_test = [Greens rmp; Smooth/mean(max(Smooth,[],2))*0/size(Smooth,1) zeros(size(Smooth,1),size(rmp,2))];

for j = 1:1:length(W)
    Green_test(j,:) = Green_test(j,:)*W(j);
end

D_test = [D_t.*W;zeros(size(Green_all,1)-size(D_t,1),1)];

[U_test_result,~,~,~] = lsqlin(Green_test,D_test,[],[],[],[],lb,ub);

Urmp_test = U_test_result(end-size(rmp,2)+1:end);
Us_test = U_test_result(1:2:end-size(rmp,2));
Ud_test = U_test_result(2:2:end-size(rmp,2));

    ms = ([Greens rmp]*U_test_result-D_t).^2;
    ms0 = D_t.^2;
    
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
        end
    end
    pct = 100*(ms0_i-ms_i)./ms0_i;
    fprintf(' rms misfit (dat., res.) and the reduction percentage = %e %e (%d%%) \n',sum(ms0),sum(ms),round(100*(sum(ms0)-sum(ms))/sum(ms0)));
    fprintf(' reduction percentage for each dataset = %d%% %d%% %d%% %d%%, only %d%% to go!\n',round(pct(1)),round(pct(2)),...
        round(pct(3)),round(pct(4)),400-round(sum(pct)));
    %fprintf(' sqrt(resnorm), mean(residual). = %e %e \n',sqrt(resnorm), mean(residual));     
    fprintf('exitflag = %d\n', exitflag);


write_inv('test.inv',XS,YS,ZS,LL,WW,DIP,STRIKE,Us_test,Ud_test,0,num_grid);
plot_patches('test.inv',13,6)




