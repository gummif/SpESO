function I = support(f,epsilon)
% support of f, supp(f)
% returns vector I such that f(I) is non-zero, and f(setdiff(1:N,I)) are
% zero
% f: vector
% epsilon: optional, support above a threshold epsilon

if nargin<2
	epsilon=0;
end


I=1:length(f);
I=I(abs(f)>epsilon);


end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
