% red noise simulator
%      Xiaohua Xu                                                
% modified from : M.Sc. Eng. Hristo Zhivomirov       

function y = rednoise(N1,N2,r)

	if rem(N1,2)
    	M1 = N1+1;
	else
    	M1 = N1;
	end

	if rem(N2,2)
    	M2 = N2+1;
	else
    	M2 = N2;
	end
	x = randn(M1, M2);
	X = fftshift(fft2(x));
	n1 = [[-M1/2:1:-1],[1:1:M1/2]]';
	n2 = [[-M2/2:1:-1],[1:1:M2/2]];

	kk = sqrt(repmat(n1,1,M2).^2+repmat(n2,M1,1).^2).^(r);
	X = X./kk;

	y = ifft2(fftshift(X));

	y = real(y(1:N1, 1:N2));

	y = y - mean(mean(y));
	yrms = sqrt(mean(mean(y.^2)));
	y = y/yrms;

end