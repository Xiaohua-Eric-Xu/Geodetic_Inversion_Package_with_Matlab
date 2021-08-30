% loop over smooth factor
% close all
SF0 = SF;
N = 30;
pct0 = [];
Greens0 = []; Green_all1 = [];
for j = 1:1:length(W)
    Greens0(j,:) = [Greens(j,:) rmp(j,:)]*W(j);
end

for k = 1:1:N
    SF1 = SF0*exp((k-round(N/2))/2);
    Green_all1 = [Greens0 ; Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1) zeros(size(Smooth,1),size(rmp,2))];

    [U,resnorm,residual,exitflag] = lsqlin(Green_all1,D_all,[],[],[],[],lb,ub);
    ms = ([Greens rmp]*U-(D-rmp0*Urmp0)).^2;
    ms0 = D.^2;
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            rms0(k,j) = sqrt(sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))))/num_dat(j+1));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
        end
    end
    pct0(k,:) = 100*(ms0_i-ms_i)./ms0_i;
    
    %Urmp = U(end-size(rmp,2)+1:end);
    %Us = U(1:2:end-size(rmp,2));
    %Ud = U(2:2:end-size(rmp,2));
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,Us,Ud,0,num_gridD);
    %m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    %disp([num2str(k),' Mw is ', num2str(m)]);
    
end
%%

figure
semilogx(SF0*exp(([1:N]-round(N/2))/2),pct0(:,1),'b-x'), hold on,grid on
semilogx(SF0*exp(([1:N]-round(N/2))/2),pct0(:,2),'r-x')
semilogx(SF0*exp(([1:N]-round(N/2))/2),pct0(:,3),'m-x')
semilogx(SF0*exp(([1:N]-round(N/2))/2),pct0(:,4),'g-x')
semilogx(SF0*exp(([1:N]-round(N/2))/2),sum(pct0')/4,'k-o')
semilogx(SF,pct(1),'b*',SF,pct(2),'r*',SF,pct(3),'m*',SF,pct(4),'g*')
figure
semilogx(SF0*exp(([1:N]-round(N/2))/2),rms0(:,1),'b-x'), hold on,grid on
semilogx(SF0*exp(([1:N]-round(N/2))/2),rms0(:,2),'r-x')
semilogx(SF0*exp(([1:N]-round(N/2))/2),rms0(:,3),'m-x')
semilogx(SF0*exp(([1:N]-round(N/2))/2),rms0(:,4),'g-x')
semilogx(SF0*exp(([1:N]-round(N/2))/2),sum(rms0')/4,'k-o')
