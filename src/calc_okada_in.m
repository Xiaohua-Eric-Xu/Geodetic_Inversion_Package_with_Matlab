function [ux,uy,uz] = calc_okada_in(U,x,y,z,nu,dip,strike,c,len,W,fault_type)

    cosd  = cos(dip);
    sind  = sin(dip);
    
    c=c+W*sind;
    x=x-W*cosd*cos(strike);
    y=y+W*cosd*sin(strike);
    
    % reversed sign of z somehow
    d=c-z; 
    strike2 = -strike+pi/2;
    coss  = cos(strike2);
    sins  = sin(strike2);
    %rot = [coss -sins ; sins coss];
    rotx =  x*coss+y*sins;
    roty = -x*sins+y*coss;
    
    L     = len/2;
    Cnst = U/(2*pi);

    p = roty*cosd + d*sind;	%a matrix eqn. (30)
    q = roty*sind - d*cosd;	%a matrix eqn. (30)
    a = 0.5/(1-nu);		% (lambda+mu)/(lambda+2mu) = 0.5/(1-poisson's ratio)
    
    [f11a,f21a,f31a,f11b,f21b,f31b,f11c,f21c,f31c] = fBi_in(rotx+L,p,a,dip,fault_type,q,z);
    [f12a,f22a,f32a,f12b,f22b,f32b,f12c,f22c,f32c] = fBi_in(rotx+L,p-W,a,dip,fault_type,q,z);
    [f13a,f23a,f33a,f13b,f23b,f33b,f13c,f23c,f33c] = fBi_in(rotx-L,p,a,dip,fault_type,q,z);
    [f14a,f24a,f34a,f14b,f24b,f34b,f14c,f24c,f34c] = fBi_in(rotx-L,p-W,a,dip,fault_type,q,z);
    u1a = f11a - f12a - f13a + f14a;
    u2a = f21a - f22a - f23a + f24a;
    u3a = f31a - f32a - f33a + f34a;
    u1b = f11b - f12b - f13b + f14b;
    u2b = f21b - f22b - f23b + f24b;
    u3b = f31b - f32b - f33b + f34b;
    u1c = f11c - f12c - f13c + f14c;
    u2c = f21c - f22c - f23c + f24c;
    u3c = f31c - f32c - f33c + f34c;
    
    d=c+z;
    p = roty*cosd + d*sind;	%a matrix eqn. (30)
    q = roty*sind - d*cosd;
    
    [f11a,f21a,f31a,~,~,~,~,~,~] = fBi_in(rotx+L,p,a,dip,fault_type,q,-z);
    [f12a,f22a,f32a,~,~,~,~,~,~] = fBi_in(rotx+L,p-W,a,dip,fault_type,q,-z);
    [f13a,f23a,f33a,~,~,~,~,~,~] = fBi_in(rotx-L,p,a,dip,fault_type,q,-z);
    [f14a,f24a,f34a,~,~,~,~,~,~] = fBi_in(rotx-L,p-W,a,dip,fault_type,q,-z);
    u1acir = f11a - f12a - f13a + f14a;
    u2acir = f21a - f22a - f23a + f24a;
    u3acir = f31a - f32a - f33a + f34a;
    
    uxj = Cnst*(u1a-u1acir+u1b+z.*u1c);
    uyj = Cnst*((u2a-u2acir+u2b+z.*u2c).*cosd - (u3a-u3acir+u3b+z.*u3c).*sind);
    uz = Cnst*((u2a-u2acir+u2b-z.*u2c).*sind + (u3a-u3acir+u3b-z.*u3c).*cosd);
    
    ux = (uxj*coss-uyj*sins);
    uy = (uxj*sins+uyj*coss);
    
    
end


