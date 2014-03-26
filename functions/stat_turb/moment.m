function m = moment(x,n,c)
% the nth sample moment
% if c is something: calculates the central moment

if nargin <= 2
	m=mean(x(:).^n);
else
	m=moment(x-moment(x,1),n);
end

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
