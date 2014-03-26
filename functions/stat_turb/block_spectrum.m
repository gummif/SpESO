function varargout=block_spectrum(f,K,iN,windowt,logav)
% calculates spectrum with block averaging
% assumes non-periodicity (when iN>1) and zero mean

% f: signal vector (of length N)
% K: number of blocks, block length = N/K (actual number of blocks is iN*K-(iN-1) )
% iN: 1 for no over lap, 2 for half overlap, ...
% windowt: window type, ('hanning','hamming','box','welch')
% logav: 1 for average in log space, 2 median, 0 regular

% kvec: vector of frequencies (in interval 1 to N/2)
% spvec: vector of specrum values at kvec
% aven(1): the mean energy at k=0
% has options for output:
% {1}=spvec;
% {1}=kvec; {2}=spvec;
% {1}=kvec; {2}=spvec; {3}=aven(1);

if nargin<5
	logav=0;
end

N=length(f);
N1=floor(N/K); % block length


	if logav==2
		aven=[];
	else
		aven=zeros(N1,1);
	end
W=window_spec(N1,windowt);
jj=0; % num. blocks
for i=0:1/iN:K-1
	f2=f(1+i*N1 : (i+1)*N1).*W;
	f2en=abs(fft(f2)).^2/length(f2);
	if logav==1
		aven=aven+log(f2en);
	elseif logav==2
		aven=[aven,log(f2en)];
	else
		aven=aven+f2en;
	end
	jj=jj+1;
end


% normalization
Wss=sum(W.^2);  % devide by Wss to normalize (N/K for 'box')
if logav==1
	aven=0.5*exp(aven/jj)/Wss;
elseif  logav==2
	aven=0.5*exp(median(aven,2))/Wss;
else
	aven=0.5*aven/Wss/jj;
end


%aven=aven/K/iN;
%aven=aven*mean(abs(fhat(1:K:end)))/mean(abs(aven));


% finalize:
kvec=K*(1:N1/2);
spvec=aven(2:N1/2+1);
spvec(1:end-1)=spvec(1:end-1) + aven(end:-1:N1/2+2);
%varargout{1}=aven(1); % the mean energy at k=0

if nargout==1
	varargout{1}=spvec;
elseif nargout==2
	varargout{1}=kvec;
	varargout{2}=spvec;
elseif nargout==3
	varargout{1}=kvec;
	varargout{2}=spvec;
	varargout{3}=aven(1);
end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
