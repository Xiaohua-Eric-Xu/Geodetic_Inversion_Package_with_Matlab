function plot_res(Greens,rmp,U,D,Urmp0,rmp0,dat_ph,dat_az,dat_gps,dat_gov,xP,yP,xA,yA,xG,yG,xO,yO,output,zone,xo,yo)
    % plot residuals for the inversion
    s = 5500;
    forwards = [Greens rmp]*U + rmp0*Urmp0;
    r = forwards - D;
    ii = abs(r)>30;
    % mark some large residues
    %save ii ii
    
    %r = [Greens rmp]*U + rmp0*Urmp0;
    x = [xP;xA;xG;xO];
    y = [yP;yA;yG;yO];
    num_dat = [0 dat_ph dat_az dat_gps*3 dat_gov];
    for j = 1:1:length(dat_ph)+length(dat_az)+length(dat_gps)+length(dat_gov)
        if j <= length(dat_ph)+length(dat_az) || j > length(dat_ph)+length(dat_az)+length(dat_gps)
            if num_dat(j+1)~=0
                xx = x(sum(num_dat(1:j))+1:sum(num_dat(1:j+1)));
                yy = y(sum(num_dat(1:j))+1:sum(num_dat(1:j+1)));
                rr = r(sum(num_dat(1:j))+1:sum(num_dat(1:j+1)));
                plot_data(xx,yy,rr)
                [xx,yy] = utm2ll(xx+xo,yy+yo,zone,2);
                if output ~= 0
                    fid = fopen(['dataset',num2str(j),'.residuals'],'w');
                    for k = 1:1:length(xx)
                        fprintf(fid,'%.6f\t%.6f\t%.6f\n',xx(k),yy(k),rr(k));
                    end
                    fclose(fid);
                    ff = forwards(sum(num_dat(1:j))+1:sum(num_dat(1:j+1)));
                    fid = fopen(['dataset',num2str(j),'.forwards'],'w');
                    for k = 1:1:length(xx)
                        fprintf(fid,'%.6f\t%.6f\t%.6f\n',xx(k),yy(k),ff(k));
                    end
                    fclose(fid);
                end
            end
        else
            if num_dat(j+1)~=0
                gps_e = D(sum(num_dat(1:j))+1:sum(num_dat(1:j))+num_dat(j+1)/3);
                gps_n = D(sum(num_dat(1:j))+num_dat(j+1)/3+1:sum(num_dat(1:j))+2*num_dat(j+1)/3);
                gps_u = D(sum(num_dat(1:j))+2*num_dat(j+1)/3+1:sum(num_dat(1:j))+3*num_dat(j+1)/3);
                gps_er = r(sum(num_dat(1:j))+1:sum(num_dat(1:j))+num_dat(j+1)/3);
                gps_nr = r(sum(num_dat(1:j))+num_dat(j+1)/3+1:sum(num_dat(1:j))+2*num_dat(j+1)/3);
                gps_ur = r(sum(num_dat(1:j))+2*num_dat(j+1)/3+1:sum(num_dat(1:j))+3*num_dat(j+1)/3);
                %{
                figure, hold on
                quiver(xG(1:end/6),yG(1:end/6),gps_e(1:end/2)*s,gps_n(1:end/2)*s,0,'r')
                quiver(xG(1:end/6),yG(1:end/6),(gps_e(1:end/2)+gps_er(1:end/2))*s,(gps_n(1:end/2)+gps_nr(1:end/2))*s,0,'m')
                quiver(xG(1:end/6),yG(1:end/6),0*xG(1:end/6),gps_u(1:end/2)*s,0,'b')
                quiver(xG(1:end/6),yG(1:end/6),0*xG(1:end/6),(gps_u(1:end/2)+gps_ur(1:end/2))*s,0,'c')
                figure, hold on
                quiver(xG(1:end/6),yG(1:end/6),gps_e(end/2+1:end)*s,gps_n(end/2+1:end)*s,0,'r')
                quiver(xG(1:end/6),yG(1:end/6),(gps_e(end/2+1:end)+gps_er(end/2+1:end))*s,(gps_n(end/2+1:end)+gps_nr(end/2+1:end))*s,0,'m')
                quiver(xG(1:end/6),yG(1:end/6),0*xG(1:end/6),gps_u(end/2+1:end)*s,0,'b')
                quiver(xG(1:end/6),yG(1:end/6),0*xG(1:end/6),(gps_u(end/2+1:end)+gps_ur(end/2+1:end))*s,0,'c')
                %}
                figure, hold on
                quiver(xG(1:end/3),yG(1:end/3),gps_e*s,gps_n*s,0,'b')
                quiver(xG(1:end/3),yG(1:end/3),(gps_e-gps_er)*s,(gps_n-gps_nr)*s,0,'c')
                quiver(xG(1:end/3),yG(1:end/3),0*xG(1:end/3),gps_u*s,0,'r')
                quiver(xG(1:end/3),yG(1:end/3),0*xG(1:end/3),(gps_u-gps_ur)*s,0,'m')
               
                figure, hold on
                plot(gps_e,gps_e-gps_er,'go')
                plot(gps_n,gps_n-gps_nr,'b*')
                plot(gps_u,gps_u-gps_ur,'rx')
                a = min([gps_e;gps_n;gps_u]);
                b = max([gps_e;gps_n;gps_u]);
                plot([a,b],[a,b],'k')
                if output ~= 0
                    [xx,yy] = utm2ll(xG(1:end/3)+xo,yG(1:end/3)+yo,zone,2);
                    gps_ef = forwards(sum(num_dat(1:j))+1:sum(num_dat(1:j))+num_dat(j+1)/3);
                    gps_nf = forwards(sum(num_dat(1:j))+num_dat(j+1)/3+1:sum(num_dat(1:j))+2*num_dat(j+1)/3);
                    gps_uf = forwards(sum(num_dat(1:j))+2*num_dat(j+1)/3+1:sum(num_dat(1:j))+3*num_dat(j+1)/3);
                    fid = fopen(['gpsset',num2str(j),'.residuals'],'w');
                    for k = 1:1:length(xx)
                        fprintf(fid,'%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',xx(k),yy(k),gps_er(k),gps_nr(k),gps_ur(k));
                    end
                    fclose(fid);
                    fid = fopen(['gpsset',num2str(j),'.forwards'],'w');
                    for k = 1:1:length(xx)
                        fprintf(fid,'%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',xx(k),yy(k),gps_ef(k),gps_nf(k),gps_uf(k));
                    end
                    fclose(fid);
                end
            end
        end
    end
                
                