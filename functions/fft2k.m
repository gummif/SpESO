function k = fft2k(N,i)
% FFT2K(N,i) computes wavenumber k from location in fft matrix (any dimensional)
%  N: size of matrix N=[N1,N2,...], powers of 2
%  i: position in matrix i=[i1,i2,...]
%  k: returns components of k vector k=[k1,k2,...]

k=zeros(1:length(N));

for d=1:length(N)
	
	if i(d) <= N(d)/2+1
		k(d) = i(d)-1;
	else
		k(d) = N(d)-i(d)+1;
	end
	
end



end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
