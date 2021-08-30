    % interatively fit the data with dip and strike slip kernals.
    x = [xP;xA;xG;xO]; y = [yP;yA;yG;yO];D_for = [];
    clearvars -except D_all x y N noise mu dat_ph dat_az dat_gps dat_gov dl dw odr sG
    
    D_save = D_all;
    N = 1000;
    n=2*N^2;
    %noise = [ones(dat_ph(1),1)*2.3;ones(dat_ph(2),1)*5.4;ones(dat_ph(3),1)*4.1;sG];
    noise = [ones(dat_ph(1),1)*5;ones(dat_ph(2),1)*5;sG];
    psi = zeros(odr, n);
    load OOAA.mat
    di = 10;
    imax = 400;
    
    Mw=[]; X2 = [];
    
    for odr1 = di:di:imax
        O = OO(:,1:odr1);
        A = AA(1:odr1,1:odr1);
        
        R = diag(A);
        
        for j = odr1-di+1:1:odr1
            str = ['psi/psi',num2str(j)];
            load(str);
            psi(j,:) = tmp_psi;%/sqrt(tmp_psi*tmp_psi'*dl*dw);
        end
        % compute expansian coefficients
        a = (R.^(-0.5)).*(D_all' * O)';
        mm = (a' * psi(1:odr1,:) )';

        Us = mm(1:N^2); Us = reshape(Us,N,N);
        Ud = mm(N^2+1:2*N^2); Ud = reshape(Ud,N,N);
        
        m = 2/3*log10(mu*(sqrt(sum(sum(Us/100*dl*dw))^2+sum(sum(Ud/100*dl*dw))^2)))-6.07;
        Mw = [Mw;m];
        
        
        %x = [xP;xA;xG;xO]; y = [yP;yA;yG;yO];
        D_for = [];
        nn = 0;
        for j = 1:1:length(x)
            str1 = repmat('\b',1,nn);
            fprintf(str1);

            str1 = ['computing forward model on data ',num2str(j),' ...'];
            fprintf(str1);
            nn = length(str1);

            str = ['green/p',[num2str(1)],'_green',num2str(j)];
            load(str);
            D_for = [D_for; sum(sum(Ud.*greend+Us.*greens))];
            %D_for = [D_for; sum(sum(Ud_tol.*greend))];
        end
        x2 = sum((D_for-D_save).^2./noise./noise/length(x));
        X2 = [X2;x2];
        str = ['Mw = ',num2str(m),', X2 = ',num2str(x2)];
        fprintf('\n');
        disp(str)
        fprintf('\n');
    end
    figure,plot([di:di:imax],Mw,'b-.')
    figure,plot([di:di:imax],X2,'b-.')
    save fitting2 Mw X2
    
    
    