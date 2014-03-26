function k = fft2wavenumber(N,d)
% FFT2WAVENUMBER(N,d) computes wavenumber k for every location in fft matrix of size N^d
%  N: size of matrix dimension, powers of 2
%  d: dimensions (1 or 2)
%  k: matrix of wavenumbers  \sqrt{ \sum_i k_i^2 }


kivec=zeros(1,N);
for i=1:N
	kivec(i)=fft2k(N,i);  % k_1
end

if d==1
	k=kivec;
elseif d==2
	k=zeros(N,N);
	
	%k = sqrt(bsxfun(@plus,kivec.^2,kivec'.^2));
	for j=1:N
		Kvec=sqrt( kivec.^2+kivec(j)^2 );
		k(j,:)=Kvec;
	end
	
else
	disp('wrong d');
end



end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
