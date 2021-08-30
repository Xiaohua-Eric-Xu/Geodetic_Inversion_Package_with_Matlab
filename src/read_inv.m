function [xs,ys,zs,len,wid,dip,strike,Us,Ud,Un,num_grid] = read_inv(str)
    % this function is used to read the inversion result into Matlab
    a = load(str);
    xs = a(:,4);
    ys = a(:,5);
    zs = a(:,6);
    len = a(:,7);
    wid = a(:,8);
    dip = a(:,9);
    strike = a(:,10);
    Us = a(:,11);
    Ud = a(:,12);
    Un = a(:,13);
    
    num_grid = zeros(max(a(:,2)),1);
    num_grid(1) = find(a(:,2) == 1,1,'last');
    
    for j = 2:1:length(num_grid)
        num_grid(j) = find(a(:,2) == j,1,'last') - sum(num_grid(1:j-1));
    end
        