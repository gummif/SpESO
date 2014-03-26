function y=CS_ranpha_2d(fh,M,N,h,T,fwt,iwt,ind)
% fh: square matrix Ns*Ns=N
% CS_mat size M*N, N a power of 2 
% h: filter (2D square) of length Ns*Ns=N: complex phase randomizer
% T: 0 or 1 for transpose
% fwt and iwt: wavelet transforms (orthonormal) in 2D
% ind: vector of length M , M=Ms^2
% y = CS_ranpha*fh (square matrix)


if T==0  % original
	fh=iwt(fh);
	
	% randomize phase
	fh=fft2(fh);
	fh=real(ifft2(fh.*h));
	
	% subsampling at ind
	y=fh(ind);
	y=vec2square(y);

	
elseif T==1 % transposed

	% De-subsampling at ind
	Ns=sqrt(N);
	y=zeros(Ns,Ns);
	y(ind)=fh;
	
	
	% De-randomize phase
	y=fft2(y);
	y=real(ifft2(y./h));
	
	y=fwt(y);
	
else
	disp('T==0 or T==1');
end



end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
