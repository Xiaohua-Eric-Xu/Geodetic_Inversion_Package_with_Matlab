    clear all
    close all
%%    
% path to where the data and initial model are:
    earthquake = 'new/Kilauea';
    %earthquake = 'new/Maule';
    input_path = ['/Users/xix016/Documents/Inversion/',earthquake,'/input2/'];
    model_path = ['/Users/xix016/Documents/Inversion/',earthquake,'/model2/'];
%    earthquake = 'Hector_corse';
%    input_path = ['/Users/xix016/Documents/InSAR_Data/Inversion/stress_driven_corse/',earthquake,'/input/'];
%    model_path = ['/Users/xix016/Documents/InSAR_Data/Inversion/stress_driven_corse/',earthquake,'/model/'];
    
% set the imput config file name under the model path
    config_inv = 'config.invBestGPSInSAR';

% adding search path for the inversion package
    addpath /Users/xix016/Documents/Inversion/Inversion

% go to the model path to load the parameters and create inversion model
    gotodir=['cd ' model_path];
    eval(gotodir);

% initiate the model
    xP=[]; yP=[]; dP=[]; sP=[]; tpP=[]; look = [];              % LOS (phase) data 
    xA=[]; yA=[]; dA=[]; sA=[]; tpA=[]; PHI = [];               % AZO (amplitude) data
    xG=[]; yG=[]; dG=[]; sG=[]; tpG=[];                         % GPS data
    xO=[]; yO=[]; dO=[]; sO=[]; uO=[]; azO=[]; tpO=[];          % geological observation data
    dat_ph=[]; dat_az=[]; dat_gps=[]; dat_gov=[];

%% modify the config.inv file under model path to change the parameters
    % read the origin for the model
    xo = cell2mat(inifile(config_inv,'read',{'model','origin','xo','d'}));
    yo = cell2mat(inifile(config_inv,'read',{'model','origin','yo','d'}));
    [xo,yo]=utm2ll(xo,yo,0,1);      % convert to UTM
    X0 = cell2mat(inifile(config_inv,'read',{'model','origin','xo','d'}));
    
    % read the fault patch parameters
    dw = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','top_patch_width','d'}));
    dl = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','top_patch_length','d'}));
    inc = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','patch_increment_factor','d'}));
    ss = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','strike_slip','d'}));
    ds = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','dip_slip','d'}));
    ns = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','normal_slip','d'}));
    fault_type = [ss; ds; ns];
    
    % positivity constraints
    PSC = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','positivity_strike','d'}));
    PDC = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','positivity_dip','d'}));
    PNC = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','positivity_normal','d'}));
    PMAX = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','positivity_max','d'}));
    BC0 = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','bottom_zero_constraint','d'}));
    
    % smoothness constraints
    SF = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','smooth_factor','d'}));
    SSEG = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','smooth_between_segments','d'}));
    SDF = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','smooth_dip_over_strike','d'}));
    
    if SDF < 1
        SF = SF/SDF;
    end
    
    % zero edge constrains
    BOT = cell2mat(inifile(config_inv,'read',{'model','edge_constraints','bot','d'}));
    SIDE = cell2mat(inifile(config_inv,'read',{'model','edge_constraints','side','d'}));
    TOP = cell2mat(inifile(config_inv,'read',{'model','edge_constraints','top','d'}));
    num_side = cell2mat(inifile(config_inv,'read',{'model','edge_constraints','num_side','d'}));
    SIDEID = [];
    if SIDE ~= 0
        for j = 1:1:num_side
            SIDEID = [SIDEID; cell2mat(inifile(config_inv,'read',{'model','edge_constraints',['side' num2str(j)],'d'}))];
        end
    end
    
    % wtighs for the data
    PW = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','weight_phase','d'}));
    AW = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','weight_azi','d'}));
    GW = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','weight_gps','d'}));
    OW = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','weight_gov','d'}));
    SWP = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','switch_phase','d'}));
    
    % read the fault trace
    num_sorc = cell2mat(inifile(config_inv,'read',{'model','model_params','num_of_sources','d'}));
    % build a structure for each source
    patch(1:num_sorc) = struct('x',0,'y',0,'z',0,'len',0,'wid',0,'dip',0,'strike',0);
    
    for j = 1:1:num_sorc
        patch(j).x = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'x','d'}));
        patch(j).y = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'y','d'}));
        patch(j).z = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'z','d'}));
        patch(j).len = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'len','d'}));
        patch(j).wid = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'wid','d'}));
        patch(j).dip = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'dip','d'}));
        patch(j).strike = cell2mat(inifile(config_inv,'read',{'model',['trace' num2str(j)],'strike','d'}));
        if patch(j).strike == 90
            patch(j).strike = 89.9;
        end
    end
    
    % load the model material parameters
    nu = cell2mat(inifile(config_inv,'read',{'model','model_params','poisson_ratio','d'}));
    
    % read the smoothness parametres between segmentselse
    SID = [];
    SSID = [];
    if SSEG ~= 0
        num_seg_smooth = cell2mat(inifile(config_inv,'read',{'model','smooth','num_seg_smooth','d'}));
        for j = 1:1:num_seg_smooth
            SID = [SID; cell2mat(inifile(config_inv,'read',{'model','smooth',['smo' num2str(j)],'d'}))];
        end
        num_inter_smooth = cell2mat(inifile(config_inv,'read',{'model','smooth','num_inter_smooth','d'}));
        for j = 1:1:num_inter_smooth
            SSID = [SSID; cell2mat(inifile(config_inv,'read',{'model','smooth',['smoi' num2str(j)],'d'}))];
        end
    end
    
    % other parameters
    RMP =  cell2mat(inifile(config_inv,'read',{'inversion','inv_params','remove_ramp','d'}));
    TP = cell2mat(inifile(config_inv,'read',{'inversion','inv_params','consider_topography','d'}));
    
%% load the data and plot them
    % load the descending data
    num_des = cell2mat(inifile(config_inv,'read',{'data','data_params','num_des_sources','d'}));
    if num_des ~= 0
        for j = 1:1:num_des
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['des' num2str(j)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            %[xx,yy] = utm2ll(dat(:,1)-X0+3,dat(:,2),0,1);
            dat(:,7) = dat(:,7)/10;
            xx = xx - xo;
            yy = yy - yo;
            ii = (xx>2e5);
            dat(ii,:) = [];
            xx(ii) = [];
            yy(ii) = [];
            xP = [xP;xx];
            yP = [yP;yy];
            dP = [dP;dat(:,7)];             % displacement (in cm)  (in case of different difinition, this value might need to be reversed.)
            if TP~=0 tpP = [tpP;dat(:,3)]; else tpP = [tpP;dat(:,3)*0]; end      % topography (in m)
            sP = [sP;dat(:,8)];             % error
        	look = [look;dat(:,4:6)];       % look vector (direction cosines)
            dat_ph = [dat_ph length(xx)];        
            plot_data(xx,yy,dat(:,7));
            title(['sampled LOS data descending ',num2str(j),', cm'])
            xlabel('m')
            ylabel('m')
        end
    end
    
    % load the ascending data
    num_asc = cell2mat(inifile(config_inv,'read',{'data','data_params','num_asc_sources','d'}));
    if num_asc ~= 0
        for j = 1:1:num_asc
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['asc' num2str(j)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            %[xx,yy] = utm2ll(dat(:,1)-X0+3,dat(:,2),0,1);
            dat(:,7) = dat(:,7)/10;
            xx = xx - xo;
            yy = yy - yo;
            xP = [xP;xx];             % longitude
            yP = [yP;yy];             % latitude
            dP = [dP;dat(:,7)];             % displacement (in cm)  (in case of different difinition, this value might need to be reversed.)
            if TP~=0 tpP = [tpP;dat(:,3)]; else tpP = [tpP;dat(:,3)*0]; end         % topography (in m)
            sP = [sP;dat(:,8)];             % error
            look = [look;dat(:,4:6)];       % look vector (direction cosines)
            dat_ph = [dat_ph length(xx)];
            plot_data(xx,yy,dat(:,7));
            title(['sampled LOS data ascending ',num2str(j),', cm'])
            xlabel('m')
            ylabel('m')
        end
    end
    
    % load the azimuth data
    num_azi = cell2mat(inifile(config_inv,'read',{'data','data_params','num_azi_sources','d'}));
    if num_azi ~= 0
        for j = 1:1:num_azi
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['azi' num2str(j)],'s'}))]);
            %[xx,yy] = utm2ll(dat(:,1)-X0+3,dat(:,2),0,1);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            xx = xx - xo;
            yy = yy - yo;
            xA = [xA;xx];                   % longitude
            yA = [yA;yy];                   % latitude
            dA = [dA;dat(:,3)*100];         % displacement (in cm)  
            tpA = [tpA;zeros(size(dat(:,3)))];    % topography (in m) not considered here, can be added later
            sA = [sA;dat(:,4)];             % error
            dat_az = [dat_az length(xx)];
            plot_data(xx,yy,dat(:,3)*100);
            title(['sampled Azimuth Offset data ',num2str(j),', cm'])
            xlabel('km')
            ylabel('km')
            phi(j) = cell2mat(inifile(config_inv,'read',{'data','data_params',['phi' num2str(j)],'d'}));
        end
    end
    
    % load the GPS data
    num_gps = cell2mat(inifile(config_inv,'read',{'data','data_params','num_gps_sources','d'}));
    if num_gps ~= 0
        figure
        hold on
        title('GPS observation','Fontsize',18)
        set(gca,'box','on','Fontsize',18);
        for j = 1:1:num_gps
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['gps' num2str(1)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1)-X0+3,dat(:,2),0,1);
            ii = dat(:,1)-X0>3;
            xx(ii) = utm2ll(dat(ii,1)-X0+3,dat(ii,2),0,1)+2*utm2ll(dat(ii,1)*0+5.999999999,dat(ii,2),0,1);
            ii = dat(:,1)-X0<-3;
            xx(ii) = utm2ll(dat(ii,1)-X0+3,dat(ii,2),0,1)-2*utm2ll(dat(ii,1)*0+5.999999999,dat(ii,2),0,1);
            %xx = xx - xo;
            yy = yy - yo;
            gps_h(j) = cell2mat(inifile(config_inv,'read',{'data','data_params',['gps' num2str(j) 'h'],'d'}));
            gps_v(j) = cell2mat(inifile(config_inv,'read',{'data','data_params',['gps' num2str(j) 'v'],'d'}));
            dat_gps = [dat_gps length(xx)];
            if gps_h(j) ~= 0
                xG = [xG;xx;xx];                    % longitude
                yG = [yG;yy;yy];                    % latitude
                dG = [dG;[dat(:,3);dat(:,4)]];   % displacement in East North Up directions (in cm) 
                sG = [sG;[dat(:,6);dat(:,7)]];	% error in East North Up directions (in cm)
                if TP~=0 tpG = [tpG;[dat(:,9);dat(:,9)]]; else tpG = [tpG;[dat(:,9)*0;dat(:,9)*0]]; end   % topography (in m)
            end
            if gps_v(j) ~= 0
                xG = [xG;xx];
                yG = [yG;yy];
                dG = [dG;dat(:,5)];
                sG = [sG;dat(:,8)];
                if TP~=0 tpG = [tpG;dat(:,9)]; else tpG = [tpG;dat(:,9)*0]; end
            end
            %{
            gpsE = [gpsE;dat(:,3)];
            gpsN = [gpsN;dat(:,4)];
            gpsU = [gpsU;dat(:,5)];
            %}
            quiver(xx,yy,dat(:,3)/10,dat(:,4)/10,'r');
            quiver(xx,yy,dat(:,5)*0,dat(:,5)/10,'b');
        end
        gps_type = [gps_h; gps_v];
    end
    
    % load the geological observation data
    num_gov = cell2mat(inifile(config_inv,'read',{'data','data_params','num_gov_sources','d'}));
    if num_gov ~= 0
        for j = 1:1:num_gov
            dat = load ([input_path cell2mat(inifile(config_inv,'read',{'data','data_files',['gov' num2str(1)],'s'}))]);
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            xx = xx - xo;
            yy = yy - yo;
            xO = [xO;xx];                   % longitude
            yO = [yO;yy];                   % latitude
            dO = [dO;dat(:,3)*100];                 % displacement magnitude (in cm)
            sO = [sO;dat(:,4)];                 % error in East North Up directions (in cm)
            azO = [azO;dat(:,5)];                   % azimuth of the displacement (in degrees)
            tpO = [tpO;dat(:,6)];           % topography (in m) not considered here, can be added later
            dat_gov = [dat_gov length(xx)];
        end
    figure
    quiver(xO,yO,dO.*cosd(90-azO),dO.*sind(90-azO),'b'); hold on;
    title('Geological observation, projected to trace','Fontsize',18)
    set(gca,'box','on','Fontsize',18);
    end
    sP = sqrt(sqrt(sP));
%% generate the patch model for inversion
    XS=[]; YS=[]; ZS=[]; XB=[]; YB=[]; ZB=[]; LL=[]; WW=[]; DIP=[]; STRIKE=[]; num_grid=[];
    %{
    for j = 1:1:6
        patch(j).wid = 7000;
    end
    %}
    figure
    for j = 1:1:num_sorc
        [xs,ys,zs,xb,yb,zb,ll,ww] = generate_grid(1, patch(j), dw, dl, inc, 1);
        XS = [XS; xs];
        YS = [YS; ys];
        ZS = [ZS; zs];
        XB = [XB; xb];
        YB = [YB; yb];
        ZB = [ZB; zb];
        LL = [LL; ll];
        WW = [WW; ww];
        DIP = [DIP; zeros(size(xs))+patch(j).dip];
        STRIKE = [STRIKE; zeros(size(xs))+patch(j).strike];
        num_grid = [num_grid; length(xs)];
    end
    clear xs ys zs xb yb zb ll ww
    
    save MODEL patch XS YS ZS XB YB ZB LL WW DIP STRIKE num_grid
    
    % save the result
    %write_inv('slip.inv',XS,YS,ZS,LL,WW,DIP,STRIKE,0,0,0,num_grid);
    
    
%% Gnenrate the green's function
    GreenP = []; GreenA = []; GreenG = []; GreenO = [];
    if num_des+num_asc ~= 0
        GreenP = generate_green_p(XS,YS,ZS,LL,WW,DIP,STRIKE,xP,yP,tpP,look,nu,fault_type);
    end
    if num_azi ~= 0
        GreenA = generate_green_a(XS,YS,ZS,LL,WW,DIP,STRIKE,xA,yA,tpA,phi,nu,fault_type,dat_az);
    end
    if num_gps ~= 0
        GreenG = generate_green_g(XS,YS,ZS,LL,WW,DIP,STRIKE,xG,yG,tpG,nu,fault_type,gps_type,dat_gps);
    end
    xOO = [];yOO = [];dOO = [];azOO = [];tpOO = [];sOO = [];
    if num_gov ~= 0
        [xOO,yOO,dOO,azOO,tpOO,sOO,dat_gov] = proj_gov2trace(patch,xO,yO,dO,azO,tpO,sO,dat_gov,150,500);
        GreenO = generate_green_o(XS,YS,ZS,LL,WW,DIP,STRIKE,xOO,yOO,azOO,tpOO,nu,fault_type);
    end
    Greens = [GreenP;GreenA;GreenG;GreenO];
    
%% create the smoothness matrix
    Smooth = [];
    if SF ~= 0
        %load Smooth
        Smooth = generate_smoothness(XS,YS,ZS,LL,WW,STRIKE,num_grid,fault_type,SDF,SSEG,SID,SSID,TOP,BOT,SIDE,SIDEID,1);
    end
    
%% add ramp to the dataset
    rmp = [];
    if RMP == 1
        rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,100); 
        % the last num should be typical value of rmp/ (ratio of cond number 
        % * typical value of GreenP), in order not to weigh the rmp too
        % much
    end
    
%%% prepare the observation vector for inversion
    % weighing inside each dataset while weighing between datasets
    % Hector_corse 
    % AW = 0.4; GW = 0.3; OW = 0;
    % SF = 0.12;
    % Landers
    % AW = 0.3; GW = 0.45; OW = 0.08;
    % SF = 0.2;
    % SF = 0.13;
    % SF = 0.17;
    % SF = 0.07;
    %GW = 0.1;
    %SF = 4;
    %SF = 0.15;
    D = []; W = [];
    if num_des+num_asc ~= 0 && PW ~= 0
        nP = mean(sP);
        sP = sP/mean(sP);
        sP = sP/PW*length(sP);
        D = [D; SWP*dP];
        W = [W; 1./sP];
    elseif num_des+num_asc ~= 0 && PW == 0
        D = [D; SWP*dP];
        W = [W; sP*0];    
    end
    if num_azi ~= 0 && AW ~= 0
        nA = mean(sA);
        sA = sA/mean(sA);
        sA = sA/AW*length(sA);
        D = [D; dA];
        W = [W; 1./sA];
    elseif num_azi ~= 0 && AW == 0
        D = [D; dA];
        W = [W; sA*0];
    end
    if num_gps ~= 0 && GW ~= 0
        nG = mean(sG);
        sG = sG/mean(sG);
        sG = sG/GW*length(sG);
        D = [D; dG];
        W = [W; 1./sG];
    elseif num_gps ~=0 && GW == 0
        D = [D; dG];
        W = [W; sG*0];
    end
    if num_gov ~= 0 && OW ~= 0
        nO = mean(sO);
        sOO = sOO/mean(sOO);
        sOO = sOO/OW*length(sOO);
        D = [D; dOO];
        W = [W; 1./sOO];
    elseif num_gov ~= 0 && OW == 0
        D = [D; dOO];
        W = [W; sOO*0];
    end

    D_all = W.*D;
    %Smooth = [Smooth;diag(ones(sum(num_grid)*2,1))/mean(LL)/3];
% inversion    
%%% prepare the Green's matrix for inversion
    %SF = 1;
    Green_all = [Greens rmp; Smooth/mean(max(Smooth,[],2))*SF/size(Smooth,1) zeros(size(Smooth,1),size(rmp,2))];
    for j = 1:1:length(W)
        Green_all(j,:) = Green_all(j,:)*W(j);
    end
    D_all = [D_all;zeros(size(Green_all,1)-size(D,1),1)];
    
%%% set lsqlin params and bounds
    lb = zeros(size(Green_all,2),1)-PMAX;
    ub = zeros(size(Green_all,2),1)+PMAX;
    if PSC > 0
        lb(1:sum(fault_type):end-size(rmp,2)) = 0;
    elseif PSC < 0
        ub(1:sum(fault_type):end-size(rmp,2)) = 0;
    end
    if PDC > 0
        lb(2:sum(fault_type):end-size(rmp,2)) = 0;
    elseif PDC < 0
        ub(2:sum(fault_type):end-size(rmp,2)) = 0;
    end

    %options = optimset('LargeScale','on','DiffMaxChange',1e-3,'DiffMinChange',1e-9,'TolCon',1e-9,'TolFun',1e-9,'TolPCG',1e-9,'TolX',1e-9,'MaxIter',1e9,'MaxPCGIter',1e9);
    [U,resnorm,residual,exitflag] = lsqlin(Green_all,D_all,[],[],[],[],lb,ub);
    if RMP == 1
        rmp0 = rmp;
        Urmp0 = U(end-size(rmp0,2)+1:end);
        D_all = [W.*(D-rmp0*Urmp0); zeros(size(Green_all,1)-size(D,1),1)];
        rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,1e6);
        Green_all = [Greens rmp; Smooth/mean(max(Smooth,[],2))*SF/size(Smooth,1) zeros(size(Smooth,1),size(rmp,2))];
        for j = 1:1:length(W)
            Green_all(j,:) = Green_all(j,:)*W(j);
        end
    
        [U,resnorm,residual,exitflag] = lsqlin(Green_all,D_all,[],[],[],[],lb,ub);
    else
        rmp0 = 0;
        Urmp0 = 0;
    end
    Urmp = U(end-size(rmp,2)+1:end);
    Us = U(1:2:end-size(rmp,2));
    Ud = U(2:2:end-size(rmp,2));
    write_inv('slip.inv',XS,YS,ZS,LL,WW,DIP,STRIKE,Us,Ud,0,num_grid);
    save ramps rmp0 Urmp0 rmp Urmp

    %plot_patches('slip.inv',13,8);

%%% plot the residuals and results

    ms = ([Greens rmp]*U-(D-rmp0*Urmp0)).^2;
    ms0 = D.^2;
    
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
    fprintf(' reduction percentage for each dataset = %f%% %d%% %d%% %d%%, only %d%% to go!\n',pct(1),round(pct(2)),...
        round(pct(3)),round(pct(4)),400-round(sum(pct)));
    %fprintf(' sqrt(resnorm), mean(residual). = %e %e \n',sqrt(resnorm), mean(residual));     
    fprintf('exitflag = %d\n', exitflag);
    
    mu = 30e9;
 %{   
    dp_list = ZS;%-WW.*sind(DIP);
    mu_list = dp_list*0;
    mu_list(dp_list>-5e3) = 21.4e9;
    mu_list(dp_list<=-5e3 & dp_list>-20e3) = 36.3e9;
    mu_list(dp_list<=-20e3 & dp_list>-35e3) = 43.2e9;
    mu_list(dp_list<=-35e3 & dp_list>-45e3) = 47.5e9;
    mu_list(dp_list<=-45e3 & dp_list>-55e3) = 68.3e9;
    mu_list(dp_list<=-55e3 & dp_list>-90e3) = 68.3e9;
    mu_list(dp_list<-90e3)=75.1e9;
%}    
    
    m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    %m = 2/3*log10(sqrt(sum(mu_list.*Us/100.*LL.*WW)^2+sum(mu_list.*Ud/100.*LL.*WW)^2))-6.07;
    str = ['seismic moment magnitude:', num2str(m)];
    disp(str);
    
    plot_patches('slip.inv',23,15);
    %2/3*log10(mu*(sqrt(sum(Ud/100.*LL.*WW)^2)))-6.07
    %2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2)))-6.07
    %hold on,plot(0,0,'ks','MarkerSize',10)
    
    plot_res(Greens,rmp,U,D,Urmp0,rmp0,dat_ph,dat_az,dat_gps,dat_gov,xP,yP,xA,yA,xG,yG,xOO,yOO,0,11,xo,yo);
    
    %x2 = sum((([Greens rmp]*U + rmp0*Urmp0 - D)./[sP/length(sG)*nG/2;sA;sG/length(sG)*nG;sOO]).^2)/length(D);
    %disp(['X2=',num2str(x2)])





