function varargout=spectrum_2d(f)
% calculates spectrum of 2D signal f
% assumes periodicity

% f:	2D signal (N*N)

% kvec: vector of wavenumbers
% E:	vector of specrum values at kvec 

% has options for output:
% {1}=E; 
% {1}=kvec; {2}=E;

N=length(f);
u=fft2(f);


% spectrum

E=zeros(N/2+2,1); % 1 for zero, rest for rounding

kivec=zeros(1,N);
for i=1:N
	kivec(i)=fft2k(N,i);  % k_1
end

for j=1:N
	Kvec=round(sqrt( kivec.^2+kivec(j)^2 ))+1; %abs(k) values + 1 for index in E
	Kvecfilt=Kvec(Kvec<=N/2);
	
	
	temp=abs(u(Kvec<=N/2,j)).^2;
	for kk=1:length(temp)
		E(Kvecfilt(kk)) = E(Kvecfilt(kk)) + temp(kk);
	end
	
	%==========
	
	% r([1,2,1])=r([1,2,1]) + [1;3;5]; gives wrong results
	%E(Kvecfilt) = E(Kvecfilt) + abs(u(Kvec<=N/2,j)).^2; %[1:N](Kvec<=N/2)=Kvec<=N/2

end

E=E/N^2;
kvec=0:length(E)-1;
if E(end-1)==0
	E(end-1)=E(end-2)/1e8;
end
if E(end)==0
	E(end)=E(end-2)/1e8;
end

if nargout<=1
	varargout{1}=E;
elseif nargout==2
	varargout{1}=kvec;
	varargout{2}=E;
end

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
