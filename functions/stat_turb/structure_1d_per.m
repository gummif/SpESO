function S = structure_1d_per(x,p,rmax)
% calculates structure function S_p(r) for periodic 1D x
% columns are measurements in space
% input:
%			x: vector of values
%			p: order (or vector of orders)
% output: 
%			S: column vector of size N-1 * length(p).

[N,~]=size(x);
if rmax>=N
	rmax=N-1;
end
S=zeros(rmax,length(p));

for r=1:rmax
	%S(r) = mean( abs(circshift(x,-r) - x).^p );
	diff=x([1+r:end,1:r]) - x;
	for i=1:length(p)
		S(r,i) = mean( abs(diff).^p(i) );
	end
end

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
