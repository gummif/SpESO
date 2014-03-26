function f = cdf97_2d(f,DIR,varargin)
% size of f: 2^J x 2^J (matrix)
% requires the cdf97_1d function
% discrete wavelet transform, CDF 9/7
% update stencil: 9, prediction stencil: 7
% vanishing moments: 4, 4
% 2nd gen. biorthogonal lifting, periodic or symmetric
% DIR: direction of transform: 1=forward, -1=backward
% varargin: 
%			(1)K: number of transform steps (default=J (Inf))
%			(2)type2d: 'isotropic' or 'separable' transform (default='isotropic')
%			(3)persym: periodic 'per' or symmetric 'sym' transform (default='per')
%			(4)adelta: translation of constant a (default=0)

% written by: Gudmundur Adalsteinsson 2012

% http://www.mathworks.com/matlabcentral/fileexchange/5502-filter-coefficients-to-popular-wavelets
% http://www.getreuer.info/home/waveletcdf97

N=length(f);
J=log2(N);

K=J;
type2d='isotropic';
persym='per';
adelta=0;

if length(varargin)>=1
	K=min(J,varargin{1}); % use fewer steps if specified
end
if length(varargin)>=2
	type2d=varargin{2};
end
if length(varargin)>=3
	persym=varargin{3};
end
if length(varargin)>=4
	adelta=varargin{4};
end

if DIR~=1 && DIR~=-1
	disp('error: DIR = 1 or -1');
	return;
end
if ~strcmp(persym,'per') && ~strcmp(persym,'sym')
	disp('error: persym: per or sym');
	return;
end


if strcmp(type2d,'isotropic')
	
	if DIR==+1
		stvec=0:K-1;
	else
		stvec=K-1:-1:0;
	end
	
	for st=stvec
		f(1:2^(J-st),1:2^(J-st))=cdf97_1d(f(1:2^(J-st),1:2^(J-st))',DIR,1,persym,adelta)';
		f(1:2^(J-st),1:2^(J-st))=cdf97_1d(f(1:2^(J-st),1:2^(J-st)) ,DIR,1,persym,adelta);
	end
	
elseif strcmp(type2d,'separable')
	
	f=cdf97_1d(f',DIR,K,persym,adelta)';
	f=cdf97_1d(f ,DIR,K,persym,adelta);
	
else
	disp('error: type2d: isotropic or separable')
end


end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
