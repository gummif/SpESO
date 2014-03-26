function plot_coef(f,x,e)
% plot coefs >= e in f at locations x
% use with figure;hold on;

N=length(f);
J=log2(N);

ms=8;

%figure;hold on;

for j=0:J-1
	ivec=1+2^j:2^(j+1);
	ind=abs(f(ivec))>=e;
	indx=2^(J-1-j):2^(J-j):N;
	if sum(ind)>0
		plot(x(indx(ind)),j,'.k','MarkerSize',ms);
	end

end
%ylabel('level');
%xlabel x;
axis([x(1),x(end),-0.5,J-0.5]);

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
