function varargout = synthetic_signal(N,dim,param)
% Create periodic synthetic signal with power-law and exponential spectrum
% N:	signal length (N*1 in 1d or N*N in 2d)
% dim:	'1d' or '2d'
% param: struct with fields (e.g. param.k0=5)
% k0:	frequency location 0
% p0:	power-law k^p0 in k<k0
% k1:	frequency location 1
% p1:	power-law k^p1 in k0<=k<k1
% Ce:	exponential constant k^p1*exp(-Ce*(k-k1)^pe) in k1<=k
% pe:	exponential power    k^p1*exp(-Ce*(k-k1)^pe) in k1<=k
% C1:	constant for power-law in k0<=k<k1 (C0 chosen s.t. E is continuous)
% n:	noise deviation from perfect spectrum ( E=E*(1+n*randn(.,.)) )
% nd:	(optional) noise variance n decreases to n/nd with k
% fwt:	forward wavelet transform operator
% iwt:	inverse wavelet transform

% f:	synthetical signal with spectrum E(k) (if n=0)
% E(k)=
%		C0*k^p0						k<=k0
%		C1*k^p1						k0<k<=k1
%		C1*k^p1*exp(-Ce*(k-k1)^pe)	k1<k
% kvec: vector of frequencies
% Has options for output:
% {1}=f;
% {1}=f; {2}=E;
% {1}=f; {2}=kvec; {3}=E;

% examples:
%synp=struct('k0',4,'p0',1,'k1',70,'p1',-5/3,'Ce',0.002,'pe',1.2,'C1',1e6,'n',0.1);
%[fs,kvecs,Es] = synthetic_signal(1024,'1d',synp);
%figure,plot(fs)
%figure,loglog(kvecs,Es)
%synp=struct('k0',50,'p0',1,'k1',70*10*8,'p1',-5/3,'Ce',0.002/10/8,'pe',1.2,'C1',1e6,'n',0.1);
%[fs,kvecs,Es] = synthetic_signal(1024*64,'1d',synp);

% synp=struct('k0',4,'p0',1,'k1',70,'p1',-5/3,'Ce',0.002,'pe',1.2,'C1',1e6,'n',0.1);
% synp.nd=3;
% [fs,kvecs,Es] = synthetic_signal(1024,'2d',synp);
% figure,imagesc(fs)
% figure,loglog(kvecs,Es)
% figure,imagesc(abs(fftshift(fft2(fs))).^(1/6))

k0=param.k0; p0=param.p0;
k1=param.k1; p1=param.p1;
Ce=param.Ce; pe=param.pe;
C1=param.C1; n=param.n;
try
    nd = param.nd;
catch ME
    nd = 1;
end


C0=C1*k0^(p1-p0);

kvec=0:N/2;

E= [C0*kvec(1:k0).^p0, ...
	C1*kvec(k0+1:k1).^(p1), ...
	C1*kvec(k1+1:end).^(p1).*exp(-Ce*(kvec(k1+1:end)-k1).^pe)]';
% add noise
noise=n*randn(length(kvec),1);
% more smooth with increasing k
noise=noise./linspace(1,nd,length(noise))';

E=E.*(1+noise);
E(E<0)=0;



if strcmp('1d',dim)==1
	
	
	% random phase
	fhat=sqrt(E).*exp(1i*rand(N/2+1,1)*2*pi);
	% conjugate correction
	fhat(N)=0;
	fhat(1)=real(fhat(1));
	fhat(N/2+1)=real(fhat(N/2+1));
	fhat(N/2+2:end)=conj(fhat(N/2:-1:2));
	
	f=real(ifft(fhat));
	
	
	
elseif  strcmp('2d',dim)==1
	
	
	fhat=zeros(N,N);
	
	kivec=zeros(1,N);
	for i=1:N
		kivec(i)=fft2k(N,i);  % k_1 , actual frequencies
	end
	
	Kmat=zeros(N,N);
	for j=1:N
		Kmat(:,j)=round(sqrt( kivec.^2+kivec(j)^2 ))+1; % freq. index for E
	end
	
	% random phase and amplitude
	ind=find(Kmat<=N/2);
	fhat(ind)=randn(length(ind),1).*exp(1i*rand(length(ind),1)*2*pi);
	

	
	
	
	
	
% 	%fhat(ind)=laplace_dist(length(ind),1,0,1).*exp(1i*rand(length(ind),1)*2*pi);
% fhat(ind)=laplace_dist(length(ind),1,0,1).*exp(1i*laplace_dist(length(ind),1,0,1)*2*pi);
% 	%fhat(1:2:end,1:2:end)=2+randn(512,512)*0.4;
% 	
% 	evo = xyindex2x(1:2:N,1:2:N,N,N); %to change
% 	%evo = randperm(N*N);
% 	%evo = evo(1:N/8);
% 	omeg = intersect(evo,find(Kmat>N/16)); 
% 	omi = randperm(length(omeg));
% 	omeg=omeg(omi(1:length(omi)/2));
% 	%fhat(omeg)=0.8*fhat(omeg)+0.2*(2+randn(length(omeg),1)'*1.3);
% 	
% 	evo = xyindex2x(2:3:N,2:3:N,N,N); %to change
% 	omeg = intersect(evo,find(Kmat>N/16)); 
% 	omi = randperm(length(omeg));
% 	omeg=omeg(omi(1:length(omi)/2));
% 	fhat(omeg)=0.5*fhat(omeg)+0.5*(0.2+randn(length(omeg),1)'*0.8);	
% 	
% 		evo = xyindex2x(2:2:N,1:4:N,N,N); %to change
% 	omeg = intersect(evo,find(Kmat>N/16)); 
% 	omi = randperm(length(omeg));
% 	omeg=omeg(omi(1:length(omi)/2));
% 	fhat(omeg)=0.8*fhat(omeg)+0.2*(-2.5+randn(length(omeg),1)'*0.8);	

	% conjugate correction
	fhat([1,N/2+1],[1,N/2+1])=abs(fhat([1,N/2+1],[1,N/2+1])); %4 reals
	fhat(N/2+2:N,1)=conj(fhat(N/2:-1:2,1)); % left col
	fhat(1,N/2+2:N)=conj(fhat(1,N/2:-1:2)); % upper row
	fhat(N/2+2:N,N/2+1)=conj(fhat(N/2:-1:2,N/2+1)); % middle col
	fhat(2:N,N/2+2:N)=conj(rot90(fhat(2:N,2:N/2),2)); %other

	% scale to fit E
	for kk=1:N/2
		ind=find(Kmat==kk);
		fhat(ind)=fhat(ind)/sqrt(sum(abs(fhat(ind)).^2))*sqrt(E(kk));
	end

	f=real(ifft2(fhat));


	%%break
% 	
% 	% wavelet stuff
% 	
% 	%figure,imagesc((abs(param.fwt(f))).^(1/6))
% 	ph_old=angle(fft2(f));
% 	abs_old=abs(fft2(f));
% 	f=param.fwt(f);
% 	
% 	%Z=f(N/16+1:N/8,N/16+1:N/8);
% 	Z=abs(randn(N/16/2,N/16/2)).^3;
% 	x=linspace(1,N/16/2,N/16);
% 		[X,Y]=meshgrid(x,x);
% 		Z = interp2(Z,X,Y,'spline')/mean(mean(abs(Z)));
% 		
% 	im=1;
% 	for d=[8,4,2]
% 		im=im*2;
% 		x=linspace(1,N/d/2,N/d);
% 		[X,Y]=meshgrid(x,x);
% 		Z = interp2(Z,X,Y,'spline')/mean(mean(abs(Z)));
% 		f(N/d+1:N/d*2,N/d+1:N/d*2)=f(N/d+1:N/d*2,N/d+1:N/d*2).*Z/im;
% 		f(N/d+1:N/d*2,1:N/d)=f(N/d+1:N/d*2,1:N/d).*Z/im;
% 		f(1:N/d,N/d+1:N/d*2)=f(1:N/d,N/d+1:N/d*2).*Z/im;
% 		%figure,imagesc((abs(f)).^(1/6))
% 	end
% 	
% 	%f = sign(f).*f.^(1.2);
% 	%f = f.*(1+randn(N,N)*0.3);
% 	%f = f.*laplace_dist(N,N,1,0.2);
% 	%figure,imagesc((abs(f)).^(1/6))
% 	f=param.iwt(f);
% 	f=abs_old.*exp(1i*angle(fft2(f)));
% 	f(1:N/8,1:N/8)=abs_old(1:N/8,1:N/8).*exp(1i*ph_old(1:N/8,1:N/8));
% 	f=real(ifft2(f));
	
	
	
	
	

else
	disp('error: dim 1d or 2d')
	return;
end



if nargout==1
	varargout{1}=f;
elseif nargout==2
	varargout{1}=f;
	varargout{2}=E;
elseif nargout==3
	varargout{1}=f;
	varargout{2}=kvec;
	varargout{3}=E;
end


end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
