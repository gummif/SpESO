function Z = cauchy(n1,n2,x0,gamma)
% Z: Cauchy distibuted with location x0 and scale gamma>0
%    of size n1*n2

X=randn(n1,n2);
Y=randn(n1,n2);

Z=X./Y;
Z=x0+gamma.*Z;

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
