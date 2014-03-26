function y=CS_rancon_db_2d(fh,M,N,h,T,fwt,iwt)
% fh: square matrix Ns*Ns=N
% CS_mat size M*N, N and M a squared (Ns a power of 2, any Ms (optimally M-K+3 a power of 2)) 
% h: filter (2D square)
% T: 0 or 1 for transpose
% fwt and iwt: wavelet transforms (orthonormal) in 2D
% y = CS_rancon*fh

Ns=sqrt(N);
Ms=sqrt(M);

K=length(h);
Nfilt=Ns+K-1;
do1=ceil((Nfilt-2)/Ms);  % 1/delta, subsampling rate
ind=2:do1:Nfilt-1;
while Ms>length(ind)
	do1=do1-1;
	ind=2:do1:Nfilt-1;
end



if T==0  % original
	fh=iwt(fh);
	
	y=conv2(fh,rot90(h,2));
	%y=y(K:Ns/Ms:end,K+1-sh:Ns/Ms:end);
	y=y(ind(1:Ms),ind(1:Ms));		% subsample
	
elseif T==1 % transposed

	fh2=zeros(Ns+K,Ns+K);
	%fh2(1:Ns/Ms:Ns,1:Ns/Ms:Ns)=fh;
	fh2(ind(1:Ms),ind(1:Ms))=fh;
	
	y=conv2(fh2,h);
	%y=fwt(y(1:Ns,1:Ns));
	y=fwt(y(K:Ns+K-1,K:Ns+K-1));
	
	
else
	disp('T==0 or T==1');
end



end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
