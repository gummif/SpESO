function f = moving_average( f , logsc, ns)
%MOVING_AVERAGE of f. Leaves the 2 elements on boundaries alone.
% if logsc is defined, use log scale average.
% if ns is defined, only average ns from boundary

N=length(f);
ind=3:N-2;
if nargin==1
	f(ind)=0.15*f(ind-2)+0.2*f(ind-1)+0.3*f(ind)+0.2*f(ind+1)+0.15*f(ind+2);
	%ind=4:N-3;
	%f(ind)=0.2*f(ind-3)+0.1*f(ind-1)+0.4*f(ind)+0.1*f(ind+1)+0.2*f(ind+3);
else
	if nargin==2
		ns=0;
	end
	ind=3+ns:N-2-ns;
	f=log(f);
	f(ind)=0.10*f(ind-2)+0.2*f(ind-1)+0.4*f(ind)+0.2*f(ind+1)+0.10*f(ind+2);
	f=exp(f);
end

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
