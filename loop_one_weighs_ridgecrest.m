close all

%AW = 0.5; GW = 0.2; OW = 0.01;SF = 0.16;
%Green_all0 = [Greens rmp];
%SM = Smooth/mean(max(Smooth,[],2))*SF/size(Smooth,1);
SF1 = SF;
%% azimuth offset
N = 30;pct1= [];
if AW ~= 0
for k = 1:1:N
    
    AW1 = AW*exp((k-round(N/2))/5); 
    GW1 = GW;
    OW1 = OW;
    W1 = [];
    if num_des+num_asc ~= 0 && PW ~= 0
        W1 = [W1; 1./sP];
    elseif num_des+num_asc ~= 0 && PW == 0
        W1 = [W1; sP*0];    
    end
    if num_azi ~= 0 && AW ~= 0
        W1 = [W1; 1./sA/AW*AW1];
    elseif num_azi ~= 0 && AW == 0
        W1 = [W1; sA*0];
    end
    if num_gps ~= 0 && GW ~= 0
        W1 = [W1; 1./sG1/GW*GW1];
        W1 = [W1; 1./sG2/GW*GW1];
    elseif num_gps ~=0 && GW == 0
        W1 = [W1; sG*0];
    end
    if num_gov ~= 0 && OW ~= 0
        W1 = [W1; 1./sOO/OW*OW1];
    elseif num_gov ~= 0 && OW == 0
        W1 = [W1; sOO*0];
    end
    for j = 1:1:length(W1)
        Greens0(j,:) = [Greens(j,:) rmp(j,:)]*W1(j);
    end
    
    %Green_all1 = [Green_all1; SM zeros(size(Smooth,1),size(rmp,2))];
    Green_all1 = [Greens0; [Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1),zeros(size(Smooth));zeros(size(Smooth)),Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1)],zeros(size(Smooth,1)*2,size(rmp,2))];
    %D_all1 = [W.*(D-rmp0*Urmp0); zeros(size(Green_all,1)-size(D,1),1)];
    D_all1 = [(D-rmp0*Urmp0).*W1; zeros(size(Green_all1,1)-size(D,1),1)];
    %[U,resnorm,residual,exitflag] = lsqlin(Green_all1,D_all1,[],[],[],[],lb,ub);
    [U,resnorm,residual,exitflag] = lsqlin([Green_all1;mmt1*dw1;mmt2*dw2],[D_all1;m1*dw1;m2*dw2],[],[],[],[],lb,ub);
    
    ms = ([Greens rmp]*U-(D-rmp0*Urmp0)).^2;
    ms0 = D.^2;
    
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps*6); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            rms1(k,j) = sqrt(sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))))/num_dat(j+1));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
            rms1(k,j) = -1;
        end
    end
    pct1(k,:) = 100*(ms0_i-ms_i)./ms0_i;
    
    %Urmp = U(end-size(rmp,2)+1:end);
    %Us = U(1:2:end-size(rmp,2));
    %Ud = U(2:2:end-size(rmp,2));
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,Us,Ud,0,num_gridD);
    %m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    %disp([num2str(k),' Mw is ', num2str(m)]);
    
end
figure
semilogx(AW*exp(([1:N]-round(N/2))/5),pct1(:,1),'b-x'), hold on,grid on
semilogx(AW*exp(([1:N]-round(N/2))/5),pct1(:,2),'r-o')
semilogx(AW*exp(([1:N]-round(N/2))/5),pct1(:,3),'m-x')
semilogx(AW*exp(([1:N]-round(N/2))/5),pct1(:,4),'g-x')
semilogx(AW,pct(1),'b*',AW,pct(2),'r*',AW,pct(3),'m*',AW,pct(4),'g*')
figure
semilogx(AW*exp(([1:N]-round(N/2))/5),rms1(:,1),'b-x'), hold on,grid on
semilogx(AW*exp(([1:N]-round(N/2))/5),rms1(:,2),'r-o')
semilogx(AW*exp(([1:N]-round(N/2))/5),rms1(:,3),'m-x')
semilogx(AW*exp(([1:N]-round(N/2))/5),rms1(:,4),'g-x')
end
%% GPS
N = 30;pct2 = [];
if GW ~= 0;
for k = 1:1:N
    
    AW1 = AW; 
    GW1 = GW*exp((k-round(N/2))/5);
    OW1 = OW;
    W1 = [];
    if num_des+num_asc ~= 0 && PW ~= 0
        W1 = [W1; 1./sP];
    elseif num_des+num_asc ~= 0 && PW == 0
        W1 = [W1; sP*0];    
    end
    if num_azi ~= 0 && AW ~= 0
        W1 = [W1; 1./sA/AW*AW1];
    elseif num_azi ~= 0 && AW == 0
        W1 = [W1; sA*0];
    end
    if num_gps ~= 0 && GW ~= 0
        W1 = [W1; 1./sG1/GW*GW1];
        W1 = [W1; 1./sG2/GW*GW1];
    elseif num_gps ~=0 && GW == 0
        W1 = [W1; sG*0];
    end
    if num_gov ~= 0 && OW ~= 0
        W1 = [W1; 1./sOO/OW*OW1];
    elseif num_gov ~= 0 && OW == 0
        W1 = [W1; sOO*0];
    end
    for j = 1:1:length(W1)
        Greens0(j,:) = [Greens(j,:) rmp(j,:)]*W1(j);
    end
    
    %Green_all1 = [Green_all1; SM zeros(size(Smooth,1),size(rmp,2))];
    Green_all1 = [Greens0; [Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1),zeros(size(Smooth));zeros(size(Smooth)),Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1)],zeros(size(Smooth,1)*2,size(rmp,2))];
    %D_all1 = [W.*(D-rmp0*Urmp0); zeros(size(Green_all,1)-size(D,1),1)];
    D_all1 = [(D-rmp0*Urmp0).*W1; zeros(size(Green_all1,1)-size(D,1),1)];
    %[U,resnorm,residual,exitflag] = lsqlin(Green_all1,D_all1,[],[],[],[],lb,ub);
    [U,resnorm,residual,exitflag] = lsqlin([Green_all1;mmt1*dw1;mmt2*dw2],[D_all1;m1*dw1;m2*dw2],[],[],[],[],lb,ub);
    
    ms = ([Greens rmp]*U-(D-rmp0*Urmp0)).^2;
    ms0 = D.^2;
    
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps*6); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            rms2(k,j) = sqrt(sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))))/num_dat(j+1));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
            rms2(k,j) = -1;
        end
    end
    pct2(k,:) = 100*(ms0_i-ms_i)./ms0_i;
    
    %Urmp = U(end-size(rmp,2)+1:end);
    %Us = U(1:2:end-size(rmp,2));
    %Ud = U(2:2:end-size(rmp,2));
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,Us,Ud,0,num_gridD);
    %m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    %disp([num2str(k),' Mw is ', num2str(m)]);
    
end
figure
semilogx(GW*exp(([1:N]-round(N/2))/5),pct2(:,1),'b-x'),hold on, grid on
semilogx(GW*exp(([1:N]-round(N/2))/5),pct2(:,2),'r-x')
semilogx(GW*exp(([1:N]-round(N/2))/5),pct2(:,3),'m-o')
semilogx(GW*exp(([1:N]-round(N/2))/5),pct2(:,4),'g-x')
semilogx(GW,pct(1),'b*',GW,pct(2),'r*',GW,pct(3),'m*',GW,pct(4),'g*')
figure
semilogx(GW*exp(([1:N]-round(N/2))/5),rms2(:,1),'b-x'), hold on,grid on
semilogx(GW*exp(([1:N]-round(N/2))/5),rms2(:,2),'r-x')
semilogx(GW*exp(([1:N]-round(N/2))/5),rms2(:,3),'m-o')
semilogx(GW*exp(([1:N]-round(N/2))/5),rms2(:,4),'g-x')
end
%% Geological
N = 30;pct3 = [];
if OW ~= 0
for k = 1:1:N
    
    AW1 = AW; 
    GW1 = GW;
    OW1 = OW*exp((k-round(N/2))/2);
    W1 = [];
    if num_des+num_asc ~= 0 && PW ~= 0
        W1 = [W1; 1./sP];
    elseif num_des+num_asc ~= 0 && PW == 0
        W1 = [W1; sP*0];    
    end
    if num_azi ~= 0 && AW ~= 0
        W1 = [W1; 1./sA/AW*AW1];
    elseif num_azi ~= 0 && AW == 0
        W1 = [W1; sA*0];
    end
    if num_gps ~= 0 && GW ~= 0
        W1 = [W1; 1./sG1/GW*GW1];
        W1 = [W1; 1./sG2/GW*GW1];
    elseif num_gps ~=0 && GW == 0
        W1 = [W1; sG*0];
    end
    if num_gov ~= 0 && OW ~= 0
        W1 = [W1; 1./sOO/OW*OW1];
    elseif num_gov ~= 0 && OW == 0
        W1 = [W1; sOO*0];
    end
    for j = 1:1:length(W1)
        Greens0(j,:) = [Greens(j,:) rmp(j,:)]*W1(j);
    end
    
    %Green_all1 = [Green_all1; SM zeros(size(Smooth,1),size(rmp,2))];
    Green_all1 = [Greens0; [Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1),zeros(size(Smooth));zeros(size(Smooth)),Smooth/mean(max(Smooth,[],2))*SF1/size(Smooth,1)],zeros(size(Smooth,1)*2,size(rmp,2))];
    %D_all1 = [W.*(D-rmp0*Urmp0); zeros(size(Green_all,1)-size(D,1),1)];
    D_all1 = [(D-rmp0*Urmp0).*W1; zeros(size(Green_all1,1)-size(D,1),1)];
    %[U,resnorm,residual,exitflag] = lsqlin(Green_all1,D_all1,[],[],[],[],lb,ub);
    [U,resnorm,residual,exitflag] = lsqlin([Green_all1;mmt1*dw1;mmt2*dw2],[D_all1;m1*dw1;m2*dw2],[],[],[],[],lb,ub);
    
    ms = ([Greens rmp]*U-(D-rmp0*Urmp0)).^2;
    ms0 = D.^2;
    
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps*6); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            rms3(k,j) = sqrt(sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))))/num_dat(j+1));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
            rms3(k,j) = -1;
        end
    end
    pct3(k,:) = 100*(ms0_i-ms_i)./ms0_i;
    
    %Urmp = U(end-size(rmp,2)+1:end);
    %Us = U(1:2:end-size(rmp,2));
    %Ud = U(2:2:end-size(rmp,2));
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,Us,Ud,0,num_gridD);
    %m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    %disp([num2str(k),' Mw is ', num2str(m)]);
    
end
figure
semilogx(OW*exp(([1:N]-round(N/2))/5),pct3(:,1),'b-x'),hold on, grid on
semilogx(OW*exp(([1:N]-round(N/2))/5),pct3(:,2),'r-x')
semilogx(OW*exp(([1:N]-round(N/2))/5),pct3(:,3),'m-x')
semilogx(OW*exp(([1:N]-round(N/2))/5),pct3(:,4),'g-o')
semilogx(OW,pct(1),'b*',OW,pct(2),'r*',OW,pct(3),'m*',OW,pct(4),'g*')
figure
semilogx(OW*exp(([1:N]-round(N/2))/5),rms3(:,1),'b-x'), hold on,grid on
semilogx(OW*exp(([1:N]-round(N/2))/5),rms3(:,2),'r-x')
semilogx(OW*exp(([1:N]-round(N/2))/5),rms3(:,3),'m-x')
semilogx(OW*exp(([1:N]-round(N/2))/5),rms3(:,4),'g-o')
end


