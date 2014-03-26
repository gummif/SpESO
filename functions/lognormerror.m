function e = lognormerror(E,Ec,k,ind)
% relative log norm error on log-log scale 
% of E compared to Ec (w(k)=1/k weighted)
% e = || log(E)-log(Ec)||_w/||log(Ec)||_w
% E,Ec:		vectors
% k:		x-axis coordinates of E,Ec
% ind:		index set (optional)
% e:		the error

if nargin<4
	ind=1:length(k);
end
E=log(E(ind));
Ec=log(Ec(ind));
k=k(ind);
E=E(:);
Ec=Ec(:);
k=k(:);

e=norm( (E-Ec)./sqrt(k) )/norm( (Ec)./sqrt(k) );

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
