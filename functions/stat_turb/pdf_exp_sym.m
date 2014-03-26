function f=pdf_exp_sym(x,lambda)
% calculates the PDF of Exponential(lambda) at x

f=lambda*exp(-lambda*x).*(x>=0)/2;
f=f+lambda*exp(lambda*x).*(x<0)/2;

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
