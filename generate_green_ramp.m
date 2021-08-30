function rmp = generate_green_ramp(xP,yP,xA,yA,dat_ph,dat_az,dat_gps,dat_gov,ratio)
    % adding plane ramps to inversion: ax+by+c
    % ramps are only for insar data, phase and azimuth_offset
    
    rmp = [];
    
    n_rmp = length(dat_ph) + length(dat_az);
    dat = [dat_ph dat_az];
    x = [xP;xA];
    y = [yP;yA];
    
    lnth = sum(dat_ph) + sum(dat_az) + sum(dat_gps)*3 + sum(dat_gov);
    for j = 1:1:n_rmp
        add_col = zeros(lnth,3);
        ii = sum(dat(1:j))-dat(j)+1 : sum(dat(1:j));
        add_col(ii,1) = x(ii)/ratio;
        add_col(ii,2) = y(ii)/ratio;
        add_col(ii,3) = 1;
        rmp = [rmp add_col];
        %{
        figure
        plot3(add_col(:,1),add_col(:,2),add_col(:,3),'.')
        grid on
        %}
    end
        
        
        