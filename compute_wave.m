N = 150;
dl = 6.735721933941810e+05/1000;
dw = 2.5e+05/1000;
%fdir = '/net/radar/Volumes/radar2/xxu/Inversion/new/test/model/';
fdir = '/radar/xix016/Documents/Inversion/new/test3/model/';
dir = '/Users/xix016/Desktop/test_dir/test/';
str = ['cd ',dir];
eval(str);
!rmdir psi/wl
!mkdir psi/wl
str = '';
for i = 1:1:N
    fprintf(str);
    str = ['computing for psi',num2str(i),'...'];
    nn = length(str);
    fprintf(str);
    str = ['psi/psi',num2str(i),'.mat'];
    % copy the file over
    file = [fdir,str];
    file2 = [dir,'tmp.mat'];
    str = ['!cp ',file,' ',file2];
    eval(str);
    
    load(file2);
    us = reshape(tmp_psi(1:1000000),1000,1000);
    ud = reshape(tmp_psi(1000001:2000000),1000,1000);
    ws = get_wavelength(us,dl,dw,0);
    wd = get_wavelength(ud,dl,dw,0);
    %str = ['psi/wl/wl',num2str(i)];
    %save(str,'ws','wd');
    save('tmp2.mat','ws','wd');
    str = ['!cp tmp2.mat ',fdir,'psi/wl/wl',num2str(i),'.mat'];
    eval(str);
    str = repmat('\b',1,nn);
end
fprintf('\n');
