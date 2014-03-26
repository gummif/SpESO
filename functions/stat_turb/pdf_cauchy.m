function f=pdf_cauchy(x,x0,gamma)
% calculates the PDF of Cauchy(x0,gamma) at x


f=1/pi*(gamma./((x-x0).^2+gamma^2));

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
