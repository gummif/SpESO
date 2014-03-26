function plot_coef_im(f,x)
% image plot of coefs in f at locations x
%figure;hold on;  from J-1 to 0
%figure;          from 0 to J-1

N=length(f);
J=log2(N);

A=zeros(J,N);

%figure;hold on;

for j=0:J-1
	ivec=1+2^j:2^(j+1);
% 	ind=abs(f(ivec))>=e;
% 	indx=2^(J-1-j):2^(J-j):N;
% 	if sum(ind)>0
% 		plot(x(indx(ind)),j,'.k','MarkerSize',ms);
% 	end

	int=2^(J-j);
	for i=1:2^j
		A(j+1,1+(i-1)*int:i*int)=f(ivec(i));
	end

end

%imagesc(x,0:J-1,sqrt(abs(A)).*sign(A))
imagesc(x,0:J-1,A)
colormap(gray)

%ylabel('level');
%xlabel x;
axis([x(1),x(end),-0.5,J-0.5]);
set(gca,'YTick',0:J-1)

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
