%% loop over depth and dip...one at a time...
close all
Ndep = 31;decdep = 1;Ndip = 21; decdip = 1;
patch1 = patch;
%{
patch1(1).wid = 20000;
patch1(2).wid = 20000;
patch1(3).wid = 20000;
patch1(4).wid = 20000;
patch1(5).wid = 20000;

patch1(1).dip = 123;
patch1(2).dip = 105;
patch1(3).dip = 91;
patch1(4).dip = 91;
patch1(5).dip = 95;
%}
N = length(patch);

%%

for i = 1:1:N
%% depth
%{
p = patch1;
pctdep = [];
%widset = patch1(i).wid + sign([1:Ndep]-round(Ndep/2)).*exp(abs([1:Ndep]-round(Ndep/2))/3.5)*1e3;
zset = patch(i).z + sign([1:Ndep]-round(Ndep/2)).*(exp(abs([1:Ndep]-round(Ndep/2))/3.5)-1)*200;%1000*([1:Ndep]-round(Ndep/2));
%zset(zset>0) = [];
%sset = patch(i).strike + sign([1:Ndip]-round(Ndip/2)).*(exp(abs([1:Ndip]-round(Ndip/2))/4.5)-1);

for dj = 1:1:length(zset)
    p(i).z = zset(dj);
    %p(i).strike = sset(dj);
    
    XSS=[]; YSS=[]; ZSS=[]; XBB=[]; YBB=[]; ZBB=[]; LLL=[]; WWW=[]; DIPP=[]; STRIKEE=[]; num_gridD=[];
    %figure
    for j = 1:1:num_sorc
        [xs,ys,zs,xb,yb,zb,ll,ww] = generate_grid(1, p(j), dw, dl, inc, 0);
        XSS = [XSS; xs];
        YSS = [YSS; ys];
        ZSS = [ZSS; zs];
        XBB = [XBB; xb];
        YBB = [YBB; yb];
        ZBB = [ZBB; zb];
        LLL = [LLL; ll];
        WWW = [WWW; ww];
        DIPP = [DIPP; zeros(size(xs))+patch(j).dip];
        STRIKEE = [STRIKEE; zeros(size(xs))+patch(j).strike];
        num_gridD = [num_gridD; length(xs)];
    end
    clear xs ys zs xb yb zb ll ww
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,0,0,0,num_gridD);
    
    GreenPP = []; GreenAA = []; GreenGG = []; GreenOO = [];
    if num_des+num_asc ~= 0
        GreenPP = generate_green_p(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xP,yP,tpP,look,nu,fault_type);
    end
    if num_azi ~= 0
        GreenAA = generate_green_a(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xA,yA,tpA,phi,nu,fault_type,dat_az);
    end
    if num_gps ~= 0
        GreenGG = generate_green_g(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xG,yG,tpG,nu,fault_type,gps_type,dat_gps);
    end
    if num_gov ~= 0
        %[xOO,yOO,dOO,azOO,tpOO,sOO,dat_gov] = proj_gov2trace(patch,xO,yO,dO,azO,tpO,sO,dat_gov,300,800);
        GreenOO = generate_green_o(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xOO,yOO,azOO,tpOO,nu,fault_type);
    end
    Greenss = [GreenPP;GreenAA;GreenGG;GreenOO];
    
    SmoothH = Smooth;
    lbb = lb;
    ubb = ub;
    %{
    if SF ~= 0
        SmoothH = generate_smoothness(XSS,YSS,ZSS,LLL,WWW,STRIKEE,num_gridD,fault_type,SDF,SSEG,SID,SSID,0);
    end
    %}
    rmpp = [];
    if RMP == 1
        rmpp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,100); 
    end
    
    Green_all3 = [Greenss rmpp; SmoothH/mean(max(SmoothH,[],2))*SF/size(SmoothH,1) zeros(size(SmoothH,1),size(rmpp,2))];
    for j = 1:1:length(W)
        Green_all3(j,:) = Green_all3(j,:)*W(j);
    end
    D_all3 = [W.*D; zeros(size(Green_all3,1)-size(D,1),1)];

    %{
    lbb = zeros(size(Green_all3,2),1)-PMAX;
    ubb = zeros(size(Green_all3,2),1)+PMAX;
    if PSC > 0
        lbb(1:sum(fault_type):end-size(rmp,2)) = 0;
    elseif PSC < 0
        ubb(1:sum(fault_type):end-size(rmp,2)) = 0;
    end
    if PDC > 0
        lbb(2:sum(fault_type):end-size(rmp,2)) = 0;
    elseif PDC < 0
        ubb(2:sum(fault_type):end-size(rmp,2)) = 0;
    end
    %}
    
    [U,resnorm,residual,exitflag] = lsqlin(Green_all3,D_all3,[],[],[],[],lbb,ubb);
    rmp3 = rmpp;
    Urmp3 = U(end-size(rmp3,2)+1:end);
    D_all3 = [W.*(D-rmp3*Urmp3); zeros(size(Green_all3,1)-size(D,1),1)];
    rmpp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,1e6);
    Green_all3 = [Greenss rmpp; SmoothH/mean(max(SmoothH,[],2))*SF/size(SmoothH,1) zeros(size(SmoothH,1),size(rmpp,2))];
    for j = 1:1:length(W)
        Green_all3(j,:) = Green_all3(j,:)*W(j);
    end
    [U,resnorm,residual,exitflag] = lsqlin(Green_all3,D_all3,[],[],[],[],lbb,ubb);
    
    Urmp = U(end-size(rmp,2)+1:end);
    Us = U(1:2:end-size(rmp,2));
    Ud = U(2:2:end-size(rmp,2));
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,Us,Ud,0,num_gridD);
    m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    disp([num2str(dj),' Mw is ', num2str(m)]);
    
    ms = ([Greenss rmpp]*U-(D-rmp3*Urmp3)).^2;
    ms0 = D.^2;
    
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            rms4(dj,j) = sqrt(sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))))/num_dat(j+1));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
            rms4(dj,j) = -1;
        end
    end
    pctdep(dj,:) = 100*(ms0_i-ms_i)./ms0_i;
end

widset = zset;
%widset = sset;
figure
subplot(1,2,1)
plot(widset,pctdep(:,1),'b-x'), hold on,grid on
plot(widset,pctdep(:,2),'r-x')
plot(widset,pctdep(:,3),'m-x')
plot(widset,pctdep(:,4),'g-x')
plot(widset,sum(pctdep'/4),'k-o')
plot(widset,pctdep*[PW; AW; GW; OW],'k-x')

subplot(1,2,2)
plot(widset,rms4(:,1),'b-x'), hold on,grid on
plot(widset,rms4(:,2),'r-x')
plot(widset,rms4(:,3),'m-x')
plot(widset,rms4(:,4),'g-x')
plot(widset,sum(rms4'/4),'k-o')
plot(widset,rms4*[PW; AW; GW; OW],'k-x')
%}
%% dip
p = patch1;
pctdip = [];
dipset = patch1(i).dip + sign([1:Ndip]-round(Ndip/2)).*(exp(abs([1:Ndip]-round(Ndip/2))/3.5)-1);
for dipj = 1:1:Ndip
    p(i).dip = dipset(dipj);
    if p(i).dip == 90
        p(i).dip = 89.9;
    end
    
    XSS=[]; YSS=[]; ZSS=[]; XBB=[]; YBB=[]; ZBB=[]; LLL=[]; WWW=[]; DIPP=[]; STRIKEE=[]; num_gridD=[];
    %figure
    for j = 1:1:num_sorc
        [xs,ys,zs,xb,yb,zb,ll,ww] = generate_grid(1, p(j), dw, dl, inc, 0);
        XSS = [XSS; xs];
        YSS = [YSS; ys];
        ZSS = [ZSS; zs];
        XBB = [XBB; xb];
        YBB = [YBB; yb];
        ZBB = [ZBB; zb];
        LLL = [LLL; ll];
        WWW = [WWW; ww];
        DIPP = [DIPP; zeros(size(xs))+patch(j).dip];
        STRIKEE = [STRIKEE; zeros(size(xs))+patch(j).strike];
        num_gridD = [num_gridD; length(xs)];
    end
    clear xs ys zs xb yb zb ll ww
    
    GreenPP = []; GreenAA = []; GreenGG = []; GreenOO = [];
    if num_des+num_asc ~= 0
        GreenPP = generate_green_p(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xP,yP,tpP,look,nu,fault_type);
    end
    if num_azi ~= 0
        GreenAA = generate_green_a(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xA,yA,tpA,phi,nu,fault_type,dat_az);
    end
    if num_gps ~= 0
        GreenGG = generate_green_g(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xG,yG,tpG,nu,fault_type,gps_type,dat_gps);
    end
    if num_gov ~= 0
        %[xOO,yOO,dOO,azOO,tpOO,sOO,dat_gov] = proj_gov2trace(patch,xO,yO,dO,azO,tpO,sO,dat_gov,300,800);
        GreenOO = generate_green_o(XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,xOO,yOO,azOO,tpOO,nu,fault_type);
    end
    Greenss = [GreenPP;GreenAA;GreenGG;GreenOO];
    
    SmoothH = Smooth;
    %{
    SmoothH = [];
    if SF ~= 0
        SmoothH = generate_smoothness(XSS,YSS,ZSS,LLL,WWW,STRIKEE,num_gridD,fault_type,SDF,SSEG,SID,SSID,0);
    end
    %}
    rmpp = [];
    if RMP == 1
        rmpp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,100); 
    end
    
    Green_all3 = [Greenss rmpp; SmoothH/mean(max(SmoothH,[],2))*SF/size(SmoothH,1) zeros(size(SmoothH,1),size(rmpp,2))];
    for j = 1:1:length(W)
        Green_all3(j,:) = Green_all3(j,:)*W(j);
    end
    D_all3 = [W.*D; zeros(size(Green_all3,1)-size(D,1),1)];
    
    lbb = zeros(size(Green_all3,2),1)-PMAX;
    ubb = zeros(size(Green_all3,2),1)+PMAX;
    if PSC > 0
        lbb(1:sum(fault_type):end-size(rmp,2)) = 0;
    elseif PSC < 0
        ubb(1:sum(fault_type):end-size(rmp,2)) = 0;
    end
    if PDC > 0
        lbb(2:sum(fault_type):end-size(rmp,2)) = 0;
    elseif PDC < 0
        ubb(2:sum(fault_type):end-size(rmp,2)) = 0;
    end
    
    [U,~,~,~] = lsqlin(Green_all3,D_all3,[],[],[],[],lbb,ubb);
    rmp3 = rmpp;
    Urmp3 = U(end-size(rmp3,2)+1:end);
    D_all3 = [W.*(D-rmp3*Urmp3); zeros(size(Green_all3,1)-size(D,1),1)];
    rmpp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,1e6);
    Green_all3 = [Greenss rmpp; SmoothH/mean(max(SmoothH,[],2))*SF/size(SmoothH,1) zeros(size(SmoothH,1),size(rmpp,2))];
    for j = 1:1:length(W)
        Green_all3(j,:) = Green_all3(j,:)*W(j);
    end
    [U,resnorm,residual,exitflag] = lsqlin(Green_all3,D_all3,[],[],[],[],lbb,ubb);
    
    Urmp = U(end-size(rmp,2)+1:end);
    Us = U(1:2:end-size(rmp,2));
    Ud = U(2:2:end-size(rmp,2));
    %write_inv('test.inv',XSS,YSS,ZSS,LLL,WWW,DIPP,STRIKEE,Us,Ud,0,num_gridD);
    m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    disp([num2str(dipj),' Mw is ', num2str(m)]);
    
    ms = ([Greenss rmpp]*U-(D-rmp3*Urmp3)).^2;
    ms0 = D.^2;
    
    num_dat = [0;sum(dat_ph); sum(dat_az); sum(dat_gps); sum(dat_gov)];
    for j = 1:1:length(num_dat)-1
        if num_dat(j+1)~=0
            ms_i(j) = sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            ms0_i(j) = sum(ms0(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))));
            rms5(dipj,j) = sqrt(sum(ms(sum(num_dat(1:j))+1:sum(num_dat(1:j+1))))/num_dat(j+1));
        else
            ms_i(j) = -1;
            ms0_i(j) = -1;
            rms5(dipj,j) = -1;
        end
    end
    pctdip(dipj,:) = 100*(ms0_i-ms_i)./ms0_i;

end
figure
subplot(1,2,1)
plot(dipset,pctdip(:,1),'b-x'), hold on,grid on
plot(dipset,pctdip(:,2),'r-x')
plot(dipset,pctdip(:,3),'m-x')
plot(dipset,pctdip(:,4),'g-x')
plot(dipset,sum(pctdip'/4),'k-o')
plot(dipset,pctdip*[PW; AW; GW; OW],'k-x')
subplot(1,2,2)
plot(dipset,rms5(:,1),'b-x'), hold on,grid on
plot(dipset,rms5(:,2),'r-x')
plot(dipset,rms5(:,3),'m-x')
plot(dipset,rms5(:,4),'g-x')
plot(dipset,sum(rms5'/4),'k-o')
plot(dipset,rms5*[PW; AW; GW; OW],'k-x')

end
                    