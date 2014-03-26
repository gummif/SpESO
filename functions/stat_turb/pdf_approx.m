function [xout,n] = pdf_approx(X,bins)

[n,xout]=hist(X,bins);

n=n/length(X)/(xout(2)-xout(1));

%figure;bar(xout,n,1,'w');


end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
