    clear all
    close all
    
% path to where the data and initial model are:
    earthquake = 'new/Nepal';
    input_path = ['/Users/xix016/Documents/Inversion/',earthquake,'/test_input/'];
    model_path = ['/Users/xix016/Documents/Inversion/',earthquake,'/test_model/'];
%    earthquake = 'Hector_corse';
%    input_path = ['/Users/xix016/Documents/InSAR_Data/Inversion/stress_driven_corse/',earthquake,'/input/'];
%    model_path = ['/Users/xix016/Documents/InSAR_Data/Inversion/stress_driven_corse/',earthquake,'/model/'];
    
% set the imput config file name under the model path
    config_inv = 'config.inv';

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
            xx = xx - xo;
            yy = yy - yo;
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
            [xx,yy] = utm2ll(dat(:,1),dat(:,2),0,1);
            xx = xx - xo;
            yy = yy - yo;
            gps_h(j) = cell2mat(inifile(config_inv,'read',{'data','data_params',['gps' num2str(j) 'h'],'d'}));
            gps_v(j) = cell2mat(inifile(config_inv,'read',{'data','data_params',['gps' num2str(j) 'v'],'d'}));
            dat_gps = [dat_gps length(xx)];
            if gps_h(j) ~= 0
                xG = [xG;xx;xx];                    % longitude
                yG = [yG;yy;yy];                    % latitude
                dG = [dG;[dat(:,3);dat(:,4)]/10];   % displacement in East North  directions (in cm) 
                sG = [sG;[dat(:,6);dat(:,7)]/10];	% error in East North directions (in cm)
                if TP~=0 tpG = [tpG;[dat(:,9);dat(:,9)]]; else tpG = [tpG;[dat(:,9)*0;dat(:,9)*0]]; end   % topography (in m)
            end
            if gps_v(j) ~= 0
                xG = [xG;xx];                       % longitude
                yG = [yG;yy];                       % latitude
                dG = [dG;dat(:,5)/10];              % displacement in up direction (in cm)
                sG = [sG;dat(:,8)/10];              % error in up directions (in cm)
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
    
    %{
    %% solve for the 0 order
        XS=[]; YS=[]; ZS=[]; XB=[]; YB=[]; ZB=[]; LL=[]; WW=[]; DIP=[]; STRIKE=[]; num_grid=[];
    %{
    for j = 1:1:6
        patch(j).wid = 7000;
    end
    %}
    figure
    for j = 1:1:num_sorc
        [xs,ys,zs,xb,yb,zb,ll,ww] = generate_grid(1, patch(j), dw, dl, 1, 1);
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
    
    rmp = [];
    if RMP == 1
        rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,100); 
        % the last num should be typical value of rmp/ (ratio of cond number 
        % * typical value of GreenP), in order not to weigh the rmp too
        % much
    end   
    
    %% create model
    D = []; W = [];
    if num_des+num_asc ~= 0 && PW ~= 0
        sP = sP/mean(sP);
        sP = sP/PW*length(sP);
        D = [D; SWP*dP];
        W = [W; 1./sP];
    elseif num_des+num_asc ~= 0 && PW == 0
        D = [D; SWP*dP];
        W = [W; sP*0];    
    end
    if num_azi ~= 0 && AW ~= 0
        sA = sA/mean(sA);
        sA = sA/AW*length(sA);
        D = [D; dA];
        W = [W; 1./sA];
    elseif num_azi ~= 0 && AW == 0
        D = [D; dA];
        W = [W; sA*0];
    end
    if num_gps ~= 0 && GW ~= 0
        sG = sG/mean(sG);
        sG = sG/GW*length(sG);
        D = [D; dG];
        W = [W; 1./sG];
    elseif num_gps ~=0 && GW == 0
        D = [D; dG];
        W = [W; sG*0];
    end
    if num_gov ~= 0 && OW ~= 0
        sOO = sOO/mean(sOO);
        sOO = sOO/OW*length(sOO);
        D = [D; dOO];
        W = [W; 1./sOO];
    elseif num_gov ~= 0 && OW == 0
        D = [D; dOO];
        W = [W; sOO*0];
    end

    D_all = W.*D;
 
    %% create the smoothness matrix
    Smooth = [];
    if SF ~= 0
        %load Smooth
        Smooth = generate_smoothness(XS,YS,ZS,LL,WW,STRIKE,num_grid,fault_type,SDF,SSEG,SID,SSID,TOP,BOT,SIDE,SIDEID,1);
    end
%% inversion    
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
%% regular inversion
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
    fprintf(' reduction percentage for each dataset = %d%% %d%% %d%% %d%%, only %d%% to go!\n',round(pct(1)),round(pct(2)),...
        round(pct(3)),round(pct(4)),400-round(sum(pct)));
    %fprintf(' sqrt(resnorm), mean(residual). = %e %e \n',sqrt(resnorm), mean(residual));     
    fprintf('exitflag = %d\n', exitflag);
    
    mu = 30e9;
    m = 2/3*log10(mu*(sqrt(sum(Us/100.*LL.*WW)^2+sum(Ud/100.*LL.*WW)^2)))-6.07;
    str = ['seismic moment magnitude:', num2str(m)];
    disp(str);
    
    plot_patches('slip.inv',2,15);
    plot_res(Greens,rmp,U,D,Urmp0,rmp0,dat_ph,dat_az,dat_gps,dat_gov,xP,yP,xA,yA,xG,yG,xOO,yOO,0,11,xo,yo)
    
%}
    
    
%% new data after subtraction
    %D_old = D; D = [Greens rmp]*U-(D-rmp0*Urmp0);
    %Wei = W; W = [];
    %% create model
    D = []; W = [];
    if num_des+num_asc ~= 0 && PW ~= 0
        sP = sP/mean(sP);
        sP = sP/PW*length(sP);
        D = [D; SWP*dP];
        W = [W; 1./sP];
    elseif num_des+num_asc ~= 0 && PW == 0
        D = [D; SWP*dP];
        W = [W; sP*0];    
    end
    if num_azi ~= 0 && AW ~= 0
        sA = sA/mean(sA);
        sA = sA/AW*length(sA);
        D = [D; dA];
        W = [W; 1./sA];
    elseif num_azi ~= 0 && AW == 0
        D = [D; dA];
        W = [W; sA*0];
    end
    if num_gps ~= 0 && GW ~= 0
        sG = sG/mean(sG);
        sG = sG/GW*length(sG);
        D = [D; dG];
        W = [W; 1./sG];
    elseif num_gps ~=0 && GW == 0
        D = [D; dG];
        W = [W; sG*0];
    end
    if num_gov ~= 0 && OW ~= 0
        sOO = sOO/mean(sOO);
        sOO = sOO/OW*length(sOO);
        D = [D; dOO];
        W = [W; 1./sOO];
    elseif num_gov ~= 0 && OW == 0
        D = [D; dOO];
        W = [W; sOO*0];
    end
            
    Wei = W; W = [];
    
    x = [xP;xA;xG;xO]; y = [yP;yA;yG;yO]; vec = [];
    if num_des+num_asc ~= 0
        vec = [vec;look];
    end
    if num_azi ~= 0
        PHI = [];
        for j = 1:1:length(phi)
            PHI = [PHI; zeros(dat_az(j),1)+phi/180*pi];
        end
        vec = [vec;[sin(PHI),cos(PHI),zeros(size(PHI))]];
    end
    if num_gps ~= 0
        for j = 1:1:length(num_gps)
            if gps_h(j) ~= 0
                vec = [vec;repmat([1,0,0],dat_gps(j),1)];
                vec = [vec;repmat([0,1,0],dat_gps(j),1)];
            end
            if gps_v(j) ~= 0
                vec = [vec;repmat([0,0,1],dat_gps(j),1)];
            end
        end
    end
    if num_gov ~=0
        for j = 1:1:length(num_gps)
            % compute the greens' function differently and the vector
            % should be different. see generate_green_o.m
        end
    end
        
 %% compute the green's function
    %plot_trace(patch)
    N = 1000;
    nn = 0;
    ld = 0;
    if ld == 0
        !mkdir green
        for j = 1:1:length(patch)
            L = patch(j).len; W = patch(j).wid;
            sd = sind(patch(j).dip); cd = cosd(patch(j).dip);
            ss = sind(90-patch(j).strike); cs = cosd(90-patch(j).strike);
            dl = L/N;
            dw = W/N;
            a = 0.5;
            b = dw*dl/4/pi;
            X = patch(j).x; Y = patch(j).y; Z = patch(j).z;
        
            %plot3([X - 0.5*L*cs, X + 0.5*L*cs, X + 0.5*L*cs + W*cd*ss, X - 0.5*L*cs + W*cd*ss, X - 0.5*L*cs],...
                    %[Y - 0.5*L*ss, Y + 0.5*L*ss, Y + 0.5*L*ss - W*cd*cs, Y - 0.5*L*ss - W*cd*cs, Y - 0.5*L*ss],...
                    %[Z, Z, Z - W*sd, Z - W*sd, Z],'k-o')
            %grid on, axis equal, hold on
        
            xob = x - (X - 0.5*L*cs + W*cd*ss);
            yob = y - (Y - 0.5*L*ss - W*cd*cs);
            %xob = x - (X - W*cd*ss);
            %yob = y - (Y + W*cd*cs);
            depth = -Z + W*sd;
        
            %xob = xob - W*cd*ss;
            %yob = yob + W*cd*cs;
            %depth = -ZS + W*sd;
            % rotate the coordinates
            rotx =  xob*cs + yob*ss;
            roty = -xob*ss + yob*cs;
        
            greensx = zeros(N);
            greensy = greensx;
            greensz = greensx;
            greendx = greensx;
            greendy = greensx;
            greendz = greensx;
            greens = greensx;
            greend = greensx;
            greenxx = greensx;
            greenyy = greensx;
            
            for k = 1:1:length(rotx)
                
                str = repmat('\b',1,nn);
                fprintf(str);

                str = ['computing green function...',num2str(k)];
                fprintf(str);
                nn = length(str);

                for m = 1:1:1000
                    for n = 1:1:1000
                        l = (m-1)*dl;
                        w = (n-1)*dw;
                        xx = rotx(k)-l; yy = roty(k)-w*cd; dd = depth-w*sd;
                        p = yy*cd + dd*sd;
                        q = yy*sd - dd*cd;
                        R = sqrt(xx^2+yy^2+dd^2);
                        I1 = a*yy*(1/R/(R+dd)^2-xx^2*(3*R+dd)/R^3/(R+dd)^3);
                        I2 = a*xx*(1/R/(R+dd)^2-yy^2*(3*R+dd)/R^3/(R+dd)^3);
                        I3 = a*xx/R^3-I2;
                        I4 = a*(-xx*yy*(2*R+dd)/R^3/(R+dd)^2);
                        I5 = a*(1/R/(R+dd)-xx^2*(2*R+dd)/R^3/(R+dd)^2);

                        greensx(m,n) = -2*b*(3*xx^2*q/R^5+I1*sd);
                        greensy(m,n) = -2*b*(3*xx*yy*q/R^5+I2*sd);
                        greensz(m,n) = -2*b*(3*xx*dd*q/R^5+I4*sd);

                        greendx(m,n) = -2*b*(3*xx*p*q/R^5-I3*sd*cd);
                        greendy(m,n) = -2*b*(3*yy*p*q/R^5-I1*sd*cd);
                        greendz(m,n) = -2*b*(3*dd*p*q/R^5-I5*sd*cd);

                    end
                end
                %rotate back to the coordinates
                greenxx = greensx*cs-greensy*ss;
                greenyy = greensx*ss+greensy*cs;
                greens = greenxx*vec(k,1)+greenyy*vec(k,2)+greensz*vec(k,3);

                greenxx = greendx*cs-greendy*ss;
                greenyy = greendx*ss+greendy*cs;
                greend = greenxx*vec(k,1)+greenyy*vec(k,2)+greendz*vec(k,3);

                str = ['green/p',[num2str(j)],'_green',num2str(k)];
                save(str,'greens','greend');
                
            end
        end
        fprintf('\n');
    else
        j = 1;
        L = patch(j).len; W = patch(j).wid;
        sd = sind(patch(j).dip); cd = cosd(patch(j).dip);
        ss = sind(90-patch(j).strike); cs = cosd(90-patch(j).strike);
        dl = L/N;
        dw = W/N;
    end
%% compute gamma    
    gamma = zeros(length(x));
    
    ld2 = 0;
    if ld2 == 0
        for k = 1:1:length(x)

            str = ['green/p',[num2str(1)],'_green',num2str(k)];
            load(str);
            green1s = greens*Wei(k);
            green1d = greend*Wei(k);

            for j = k:1:length(x)

                str1 = repmat('\b',1,nn);
                fprintf(str1);
                
                str1 = ['computing gamma, line ',num2str(k),' row ',num2str(j),' ...'];
                fprintf(str1);
                nn = length(str1);

                str = ['green/p',[num2str(1)],'_green',num2str(j)];
                load(str);
                green2s = greens*Wei(j);
                green2d = greend*Wei(j);
                gamma(k,j) = sum(sum(green1s.*green2s+green1d.*green2d))/dl/dw;
                % use the symmatry
                gamma(j,k) = gamma(k,j);

            end

        end
        fprintf('\n');
        save Gamma gamma
    else
        load Gamma.mat
    end
    gamma2 = gamma;
 %% Remove ramp
    % currently use 1 for all the data
    %W = W*0+1;
    tmp_str = '_save';
    rmv_rmp = 0;
    if rmv_rmp ~= 0;
    ld3 = 1;
    if ld3 == 0
        load(['Ud',tmp_str]);
        load(['Us',tmp_str]);

        D_tmp = [];
        for j = 1:1:length(x)
            str1 = repmat('\b',1,nn);
            fprintf(str1);

            str1 = ['computing forward model on data ',num2str(j),' ...'];
            fprintf(str1);
            nn = length(str1);

            str = ['green/p',[num2str(1)],'_green',num2str(j)];
            load(str);
            D_tmp = [D_tmp; sum(sum(fliplr(ud').*greend+fliplr(us').*greens))];
        end
        fprintf('\n');
        save Dtmp D_tmp
    else
        load(['Ud',tmp_str]);
        load(['Us',tmp_str]);
        load(['Dtmp',tmp_str]);
    end
    load(['Urmp',tmp_str]);
    
    rmp0 = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,100); 
    rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,1e6);
    D_rmp = rmp*Urmp + rmp0*Urmp0;
    D_tmp = D_tmp + D_rmp;
    end
    
%%    
    Wei = Wei*0+1;
    D_all = Wei.*D;
    jj = cumsum([dat_ph,dat_az,dat_gps,dat_gov]);
    jj = [0,jj];
    %for j = 1:1:length(jj)-1
    %    plot_data(x(jj(j)+1:jj(j+1)),y(jj(j)+1:jj(j+1)),D(jj(j)+1:jj(j+1))-D_tmp(jj(j)+1:jj(j+1)));
    %end
    %% spectral analysis
    gamma = gamma2.*(Wei*Wei');
    
    odr = 10;
    n = N^2*2;
    m = length(x);
    O = zeros(m, odr);
    A = zeros(odr, odr);
    psi = zeros(odr, n);
    
    [O,A] = eigs(gamma, odr);
    R = diag(A);
    
    nn = 0;
    %%
    ld4 = 0;
    if ld4 == 0
        for k = 1:1:length(x)


            str1 = repmat('\b',1,nn);
            fprintf(str1);

            str1 = ['computing psi, data point ',num2str(k),' ...'];
            fprintf(str1);
            nn = length(str1);

            str = ['green/p',[num2str(1)],'_green',num2str(k)];
            load(str);

            green1s = reshape(greens,1,N^2)*Wei(k)/dl/dw;
            green1d = reshape(greend,1,N^2)*Wei(k)/dl/dw;
            tmp1 = repmat(R.^(-0.5) .* O(k,:)',1,n);
            tmp2 = repmat([green1s,green1d],odr,1);

            psi =  psi+ tmp1.*tmp2;    
        end
        !mkdir psi
        fprintf('\n');
        for j = 1:1:odr
            str = ['psi/psi',num2str(j)];
            tmp_psi = psi(j,:);
            save(str,'tmp_psi');
        end
    else
        for j = 1:1:odr
            str = ['psi/psi',num2str(j)];
            load(str);
            psi(j,:) = tmp_psi;
        end
    end
    %%
    a = (R.^(-0.5)).*(D_all' * O)';
    mm = (a' * psi )';
    sa = R.^(-0.5);
    Us = mm(1:N^2); Us = reshape(Us,N,N);
    Ud = mm(N^2+1:2*N^2); Ud = reshape(Ud,N,N);
    mm_err = sa'*abs(psi);
    clear psi;
    Us_err = mm_err(1:N^2); Us_err = reshape(Us_err,N,N);
    Ud_err = mm_err(N^2+1:2*N^2); Ud_err = reshape(Ud_err,N,N);
    
    mu = 30e9;
    Us_tol = Us+fliplr(us');
    Ud_tol = Ud+fliplr(ud');
    m = 2/3*log10(mu*(sqrt(sum(sum(Us_tol/100*dl*dw))^2+sum(sum((Ud_tol>0).*Ud_tol/100*dl*dw))^2)))-6.07;
    str = ['Mw = ',num2str(m)];
    disp(str)
    
    
    %%
    figure
    subplot(3,2,1)
    imagesc(Ud_tol),colorbar,title('Ud')
    subplot(3,2,2)
    imagesc(Us_tol),colorbar,title('Us')
    subplot(3,2,3)
    imagesc(Ud),colorbar,title('Ud')
    subplot(3,2,4)
    imagesc(Us),colorbar,title('Us')
    subplot(3,2,5)
    imagesc(Ud_err),colorbar,title('Ud_e')
    subplot(3,2,6)
    imagesc(Us_err),colorbar,title('Us_e')
    
    
    figure
    zl = (1:1:N)*patch.wid/N*sind(patch.dip);
    plot(sum(Ud_tol)*dl/1000/100,zl,'b','LineWidth',3),hold on, grid on
    plot(sum(Ud_tol)*dl/1000/100-sum(Ud_err)*dl/1000/100,zl,'r','LineWidth',2)
    plot(sum(Ud_tol)*dl/1000/100+sum(Ud_err)*dl/1000/100,zl,'r','LineWidth',2)
    
    %% plot the model
    x0 = patch.x - 0.5*patch.len*sind(patch.strike);
    y0 = patch.y - 0.5*patch.len*cosd(patch.strike);
    z0 = patch.z;
    
    dl = patch.len/N;
    dw = patch.wid/N;
    
    stepxl = dl*sind(patch.strike);
    stepxw = dw*cosd(patch.dip)*cosd(patch.strike);
    stepyl = dl*cosd(patch.strike);
    stepyw = -dw*cosd(patch.dip)*sind(patch.strike);
    stepz = -dw*sind(patch.dip);
    
    x0 = x0 + 0.5*(stepxl+stepxw);
    y0 = y0 + 0.5*(stepyl+stepyw);
    z0 = z0 + 0.5*stepz;
    
    [xx,yy]=meshgrid(1:1:N,1:1:N);
    
    x = x0 + (xx-1)*stepxl + (yy-1)*stepxw;
    y = y0 + (xx-1)*stepyl + (yy-1)*stepyw;
    z = z0 + (yy-1)*(stepz);
    
    figure
    subplot(3,2,1)
    surf(x,y,z,ud,'EdgeColor','None'),colorbar,axis equal,view(0,90)
    
    subplot(3,2,2)
    surf(x,y,z,us,'EdgeColor','None'),colorbar,axis equal,view(0,90)
    
    subplot(3,2,3)
    surf(x,y,z,fliplr(Ud_tol)','EdgeColor','None'),colorbar,axis equal,view(0,90)
    
    subplot(3,2,4)
    surf(x,y,z,fliplr(Us_tol)','EdgeColor','None'),colorbar,axis equal,view(0,90)
    
    subplot(3,2,5)
    surf(x,y,z,fliplr(Ud_err)','EdgeColor','None'),colorbar,axis equal,view(0,90)
    
    subplot(3,2,6)
    surf(x,y,z,fliplr(Us_err)','EdgeColor','None'),colorbar,axis equal,view(0,90)
    
    
    %% res
    x = [xP;xA;xG;xO]; y = [yP;yA;yG;yO];
    D_for = [];
    for j = 1:1:length(x)
        str1 = repmat('\b',1,nn);
        fprintf(str1);

        str1 = ['computing forward model on data ',num2str(j),' ...'];
        fprintf(str1);
        nn = length(str1);

        str = ['green/p',[num2str(1)],'_green',num2str(j)];
        load(str);
        D_for = [D_for; sum(sum(Ud_tol.*greend+Us_tol.*greens))];
    end
    %%
    jj = cumsum([dat_ph,dat_az,dat_gps,dat_gov]);
    jj = [0,jj];
    for j = 1:1:length(jj)-1
        plot_data(x(jj(j)+1:jj(j+1)),y(jj(j)+1:jj(j+1)),D(jj(j)+1:jj(j+1))-D_for(jj(j)+1:jj(j+1)));
    end
    
    
    
    
    
    
    
    