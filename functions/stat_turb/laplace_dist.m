function Z=laplace_dist(n1,n2,mu,b)
% Z: distributed with Laplace(mu,b) of size n1*n2

U=rand(n1,n2)-0.5;

Z=mu-b*sign(U).*log(1-2*abs(U));

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
