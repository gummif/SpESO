function f = cdf97_1d(f,DIR,varargin)
% size of f: 2^J (column vector) (or a columnwise transform of a matrix)
% discrete wavelet transform, CDF 9/7
% update stencil: 9, prediction stencil: 7
% vanishing moments: 4, 4
% 2nd gen. biorthogonal lifting, periodic or symmetric
% DIR: direction of transform: 1=forward, -1=backward
% varargin:
%			(1)K: number of transform steps (default=J (Inf))
%			(2)persym: periodic 'per' or symmetric 'sym' transform (default='per')
%			(3)adelta: translation of constant a (default=0)

% written by: Gudmundur Adalsteinsson 2012

% http://www.mathworks.com/matlabcentral/fileexchange/5502-filter-coefficients-to-popular-wavelets
% http://www.getreuer.info/home/waveletcdf97

[N,~]=size(f);
J=round(log2(N));

K=J;
persym='per';
adelta=0;


if length(varargin)>=1
	K=min(J,varargin{1}); % use fewer steps if specified
end
if length(varargin)>=2
	persym=varargin{2};
end
if length(varargin)>=3
	adelta=varargin{3};
end

if DIR~=1 && DIR~=-1
	disp('error: DIR = 1 or -1');
	return;
end
if ~strcmp(persym,'per') && ~strcmp(persym,'sym')
	disp('error: persym: per or sym');
	return;
end



% [Do Quan & Yo-Sung Ho. Optimized median lifting scheme for lossy image compression.]
a=-1.5861343420604-adelta;
b1=-1/(4*(1+2*a)^2);
a2=-(1+2*a)^2/(1+4*a);
b2=1/16*(4+(1-8*a)/(1+2*a)^2-2/(1+2*a)^3);
ka=2*(1+2*a)/(1+4*a)*sqrt(2);					% f(1)/sqrt(N) is mean
kb =  1/ka;
a1=-a;
a2=-a2;


% wavelab has reversed roles of dwmf and qmf, or h and h^tilda (see (7.162)
% in wavelet tour) for 'CDF'
% equal to 'Villasenor' par=1 in wavelab

if DIR==1  % FORWARD
	
	if strcmp(persym,'per')
		for i=0:K-1
			N = 2^(J-i);
			
			% restriction and
			% prediction 1 (periodic)
			f(1:N,:)=[f(1:2:N-1,:);
				(f(2:2:N,:)-a1*(f(1:2:N-1,:)+f([3:2:N-1,1],:)))];
			% update 1 (periodic)
			f(1:N/2,:)=f(1:N/2,:)+b1*(f(N/2+1:N,:)+f([N,N/2+1:N-1],:));
			
			% prediction 2 (periodic)
			f(N/2+1:N,:)=f(N/2+1:N,:)-a2*(f(1:N/2,:)+f([2:N/2,1],:));
			% update 2 (periodic)
			f(1:N/2,:)=f(1:N/2,:)+b2*(f(N/2+1:N,:)+f([N,N/2+1:N-1],:));
			
			% normalize
			f(1:N/2,:)=f(1:N/2,:)*ka;
			f(N/2+1:N,:)=f(N/2+1:N,:)*kb;
		end
	else
		for i=0:K-1
			N = 2^(J-i);
			
			% restriction and
			% prediction 1 (symmetric)
			f(1:N,:)=[f(1:2:N-1,:);
				(f(2:2:N,:)-a1*(f(1:2:N-1,:)+f([3:2:N-1,N-1],:)))];
			% update 1 (symmetric)
			f(1:N/2,:)=f(1:N/2,:)+b1*(f(N/2+1:N,:)+f([N/2+1,N/2+1:N-1],:));
			
			% prediction 2 (symmetric)
			f(N/2+1:N,:)=f(N/2+1:N,:)-a2*(f(1:N/2,:)+f([2:N/2,N/2],:));
			% update 2 (symmetric)
			f(1:N/2,:)=f(1:N/2,:)+b2*(f(N/2+1:N,:)+f([N/2+1,N/2+1:N-1],:));
			
			% normalize
			f(1:N/2,:)=f(1:N/2,:)*ka;
			f(N/2+1:N,:)=f(N/2+1:N,:)*kb;
		end
	end
	
else           % BACKWARD
	if strcmp(persym,'per')
		for i=K-1:-1:0
			N = 2^(J-i);
			
			% normalize
			f(1:N/2,:)=f(1:N/2,:)/ka;
			f(N/2+1:N,:)=f(N/2+1:N,:)/kb;
			
			% update 2 (periodic)
			f(1:N/2,:)=f(1:N/2,:)-b2*(f(N/2+1:N,:)+f([N,N/2+1:N-1],:));
			% prediction 2 (periodic)
			f(N/2+1:N,:)=f(N/2+1:N,:)+a2*(f(1:N/2,:)+f([2:N/2,1],:));
			
			% update 1 (periodic)
			f(1:N/2,:)=f(1:N/2,:)-b1*(f(N/2+1:N,:)+f([N,N/2+1:N-1],:));
			
			% restriction and
			% prediction 1 (periodic)
			temp=f(1:N/2,:);
			f(2:2:N,:)=f(N/2+1:N,:)+a1*(temp+temp([2:end,1],:));
			f(1:2:N-1,:)=temp;
		end
	else
		for i=K-1:-1:0
			N = 2^(J-i);
			
			% normalize
			f(1:N/2,:)=f(1:N/2,:)/ka;
			f(N/2+1:N,:)=f(N/2+1:N,:)/kb;
			
			% update 2 (symmetric)
			f(1:N/2,:)=f(1:N/2,:)-b2*(f(N/2+1:N,:)+f([N/2+1,N/2+1:N-1],:));
			% prediction 2 (symmetric)
			f(N/2+1:N,:)=f(N/2+1:N,:)+a2*(f(1:N/2,:)+f([2:N/2,N/2],:));
			
			% update 1 (symmetric)
			f(1:N/2,:)=f(1:N/2,:)-b1*(f(N/2+1:N,:)+f([N/2+1,N/2+1:N-1],:));
			
			% restriction and
			% prediction 1 (symmetric)
			temp=f(1:N/2,:);
			f(2:2:N,:)=f(N/2+1:N,:)+a1*(temp+temp([2:end,end],:));
			f(1:2:N-1,:)=temp;
		end
	end
end
end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
