function y=CS_ranpha(fh,M,N,h,T,fwt,iwt,ind)
% fh: column vector
% CS_mat size M*N, N a power of 2 
% h: filter of length N: complex phase randomizer
% T: 0 or 1 for transpose
% fwt and iwt: wavelet transforms (orthonormal)
% ind: vector of length M (vector of length N with {-1,0,1} values for modulation)
% y = CS_ranpha*fh



if T==0  % original
	fh=iwt(fh);
	
	% randomize phase
	fh=fft(fh);
	fh=ifft(fh.*h);
	
	% subsampling at ind
	y=fh(ind);
	
	% modulated sum / subsampling
% 	fh=fh.*ind;
% 	y=sum(reshape(fh',N/M,M));
% 	y=sqrt(M/N)*y';
% 	%y=y';

	
elseif T==1 % transposed

	% De-subsampling at ind
	y=zeros(N,1);
	y(ind)=fh;
	
	% De-modulated sum / De-subsampling
% 	y=bsxfun(@times, reshape(ind,N/M,M), fh');
% 	%y=y(:)/(N/M);
% 	y=y(:)/sqrt(M/N);
	
	% De-randomize phase
	y=fft(y);
	y=ifft(y./h);
	
	y=fwt(y);
	
else
	disp('T==0 or T==1');
end



end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
