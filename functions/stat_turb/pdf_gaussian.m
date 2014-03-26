function f=pdf_gaussian(x,mu,std)
% calculates the PDF of Gaussian(mu,std) at x


f=1/(std*sqrt(2*pi))*exp(-0.5*((x-mu)./std).^2);

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
