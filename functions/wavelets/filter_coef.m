function f = filter_coef(f,n,varargin)
% f is a matrix or colomn vector
% put all to zero except the n with the largest absolute value
% leaving f(1:K,1:K) alone if specified
%         f(1:K)
%varargin		(1) K
%				(2) scale, multiply coefs by 2^scale when going to courser
%					scales -> effective epsilon scale dependent.


[N,M]=size(f);



K=0;
if ~isempty(varargin) && varargin{1}>0
	K=varargin{1};
	if M==1
		temp=f(1:K);
	else
		temp=f(1:K,1:K);
	end
end

fkp=f;
if K>1
	if M==1
		fkp(1:K)=zeros(K,1);
	else
		fkp(1:K,1:K)=zeros(K,K);
	end
end

%scale=0;
if length(varargin)==2
	scale=varargin{2};
	
	%make the larger wavelets coefs bigger/smaller
	for j=log2(N-1):-1:log2(max(K,1))
		ind=round(j);
		if M==1
			fkp(1:2^ind)=fkp(1:2^ind)*2^(scale);
		else
			fkp(1:2^ind,1:2^ind)=fkp(1:2^ind,1:2^ind)*2^(scale);
		end
	end
end


%fkp=fkp(:);


[~,i] = sort(abs(fkp(:)));
f(i(1:end-n)) = 0;			% hard threshhold


%f=steinmuThresh(f,fs(n),1.5);

if K>0
	if M==1
		f(1:K)=temp;
	else
		f(1:K,1:K)=temp;
	end
end

end

function ft=steinmuThresh(ft,Tc,mu)
ft= ft.*max(1-(Tc^mu)./(abs(ft).^mu), 0);
end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
