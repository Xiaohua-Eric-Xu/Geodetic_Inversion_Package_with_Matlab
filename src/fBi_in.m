function [f1a,f2a,f3a,f1b,f2b,f3b,f1c,f2c,f3c] = fBi_in(epsn,eta,a,dip,fault_type,q,z)
    % fault type: 1.strike-slip     2.dip-slip      3.tensile
    % setup some parameters
    tol = 1e-15;
    cosd  = cos(dip);
    sind  = sin(dip);
    %tand  = tan(delta);
    cosd2 = cosd*cosd;
    sind2 = sind*sind;
    %cssnd = cos(delta)*sin(delta);
    
    % reversed sign of z somehow
    %d=c-z;
    
    %{
    strike2 = -strike+pi/2;
    coss  = cos(strike2);
    sins  = sin(strike2);
    %rot = [coss -sins ; sins coss];
    rotx =  x*coss+y*sins;
    roty = -x*sins+y*coss;
    p = roty*cosd + d*sind;
    q = roty*sind - d*cosd;
    %}
    R     = sqrt(epsn.^2 + eta.^2 + q.^2);
    ytil  = eta*cosd + q*sind;		
    dtil  = eta*sind - q*cosd;
    X     = sqrt(epsn.^2+q.^2);
    if abs(q) < tol
        Theta = 0;
    else
        Theta = atan((epsn.*eta)./(q.*R));
    end
    ctil = dtil + z;
    
    % compute I-s
    if (abs(cosd) > 1e-15)
        if abs(R+eta) < tol
            I3 = 1./cosd.*ytil./(R+dtil) - 1./cosd2.*(-log(R-eta) - sind.*log(R+dtil));
        else
            I3 = 1./cosd.*ytil./(R+dtil) - 1./cosd2.*(log(R+eta) - sind.*log(R+dtil));
        end
        if abs(epsn) < tol
            I4 = 0;
        else
            I4 = sind/cosd*epsn./(R+dtil) + 2/cosd2.*atan( (eta.*(X+q.*cosd) + X.*(R+X).*sind) / (epsn.*(R+X).*cosd) );
        end
    else
        if abs(R+eta) < tol
            I3 = 0.5*(eta./(R+dtil) + ytil.*q./(R+dtil).^2 + log(R-eta));
        else
            I3 = 0.5*(eta./(R+dtil) + ytil.*q./(R+dtil).^2 - log(R+eta));
        end
        if abs(epsn) < tol
            I4 = 0;
        else
            I4 = 0.5*epsn.*ytil./(R+dtil).^2;
        end
    end
    I1 = -epsn/(R+dtil)*cosd-I4*sind;
    I2 = log(R+dtil) + I3.*sind;
    
    % compute h, X-s, Y-s and Z-s
    if abs(R+epsn) < tol
        X11 = 0;
        X32 = 0;
        X53 = 0;
    else
        X11 = 1 ./ (R.*(R+epsn));
        X32 = (2*R+epsn) ./ (R.^3 .* (R+epsn).^2);
        X53 = (8*R.^2 + 9*R.*epsn + 3*epsn.^2) ./ (R.^5 .* (R+epsn).^3);
    end
    if abs(R+eta) < tol
        Y11 = 0;
        Y32 = 0;
        Y53 = 0;
    else
        Y11 = 1 ./ (R.*(R+eta));
        Y32 = (2*R+eta) ./ (R.^3 .* (R+eta).^2);
        Y53 = (8*R.^2 + 9*R.*eta + 3*eta.^2) ./ (R.^5 .* (R+eta).^3);
    end
    
   
    h = q.*cosd - z;
    Z32 = sind./R.^3 - h.*Y32;
    Z53 = 3*sind./R.^5 - h.*Y53;
    Y0 = Y11 - epsn.^2.*Y32;
    Z0 = Z32 - epsn.^2.*Z53;
    
    if fault_type == 1
        f1a = Theta/2               +   a/2.*epsn.*q.*Y11; 
        f2a =                           a/2*q./R;
        if abs(R+eta) < tol
            f3a = (1-a)/2*(-log(R-eta)) -   a/2*q.^2.*Y11;
        else
            f3a = (1-a)/2*log(R+eta)    -   a/2*q.^2.*Y11;
        end
        
        f1b = -epsn.*q.*Y11-Theta   -   (1-a)/a*I1.*sind; 
        f2b = -q./R                 +   (1-a)/a*ytil./(R+dtil).*sind;
        f3b = q.^2.*Y11             -   (1-a)/a*I2.*sind;
        
        f1c = (1-a)*epsn.*Y11.*cosd          -  a*epsn.*q.*Z32; 
        f2c = (1-a)*(cosd./R+2*q.*Y11.*sind) -  a*ctil.*q./R.^3; 
        f3c = (1-a)*q.*Y11.*cosd             -  a*(ctil.*eta./R.^3 - z.*Y11 + epsn.^2.*Z32);
    elseif fault_type == 2
        f1a = a/2*q./R; 
        f2a = Theta/2               +   a/2*eta.*q.*X11;
        if abs(R+epsn) < tol
            f3a = (1-a)/2*(-log(R-epsn))   -   a/2*q.^2.*X11;
        else
            f3a = (1-a)/2*log(R+epsn)   -   a/2*q.^2.*X11;
        end
        
        f1b = -q/R                  +   (1-a)/a*I3.*sind.*cosd;
        f2b = -eta.*q.*X11-Theta    -   (1-a)/a*epsn./(R+dtil).*sind.*cosd;
        f3b = q.^2.*X11             +   (1-a)/a*I4.*sind.*cosd;
        
        f1c = (1-a)*cosd./R         -   q.*Y11.*sind-a.*ctil.*q./R.^3;
        f2c = (1-a).*ytil.*X11      -   a*ctil.*eta.*q.*X32;
        f3c = -dtil.*X11            -   epsn.*Y11.*sind-a*ctil.*(X11-q.^2.*X32);
    elseif fault_type == 3
        if abs(R+eta) < tol
            f1a = -(1-a)/2*(-log(R-eta))-a/2*q.^2.*Y11;
        else
            f1a = -(1-a)/2*log(R+eta)-a/2*q.^2.*Y11;
        end
        if abs(R+epsn) < tol
            f2a = -(1-a)/2*(-log(R-epsn))-a/2*q.^2.*X11;
        else
            f2a = -(1-a)/2*log(R+epsn)-a/2*q.^2.*X11;
        end
        f3a = Theta/2-a/2*q.*(eta.*X11+epsn.*Y11);
        
        f1b = q.^2.*Y11     -   (1-a)/a.*I3.*sind2;
        f2b = q.^2.*X11     +   (1-a)/a*epsn./(R+dtil).*sind2;
        f3b = q.*(eta.*X11+epsn.*Y11)-Theta-(1-a)/a.*I4.*sind2;
        
        f1c = -(1-a)*(sind./R+q.*Y11.*cosd)-a*(z.*Y11-q.^2.*Z32);
        f2c = (1-a)*2*epsn.*Y11.*sind+dtil.*X11-a*ctil.*(X11-q.^2.*X32);
        f3c = (1-a)*(ytil.*X11+epsn.*Y11.*cosd)+a*q.*(ctil.*eta.*X32+epsn.*Z32);
    end
    
    
    
    
    
end