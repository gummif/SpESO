function y=CS_rancon_db(fh,M,N,h,T,fwt,iwt)
% fh: column vector
% CS_mat size M*N, N a power of 2, any M (optimally M-K+3 a power of 2)
% h: filter
% T: 0 or 1 for transpose
% fwt and iwt: wavelet transforms (orthonormal)
% y = CS_rancon*fh


K=length(h);
Nfilt=N+K-1;
do1=ceil((Nfilt-2)/M);  % 1/delta, subsampling rate
ind=2:do1:Nfilt-1;
while M>length(ind)
	do1=do1-1;
	ind=2:do1:Nfilt-1;
end

if T==0  % original
	fh=iwt(fh);
	
	y=conv(fh,reverse(h)); %length(y), Nfilt
	y=y(ind(1:M));		% subsample

	
elseif T==1 % transposed

	fh2=zeros(N+K,1);
	fh2(ind(1:M))=fh;
	
	y=conv(fh2,h);
	y=fwt(y(K:N+K-1));
	

else
	disp('T==0 or T==1');
end



end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
