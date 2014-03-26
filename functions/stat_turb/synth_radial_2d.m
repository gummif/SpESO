function f=synth_radial_2d(c,N,L,x,pos)
% A 2D signal f=c(r)
% isotropic and periodic (if L chosen < radius of cusp)
% c: function handle
% N: size of f
% L: length of f for dx
% x: centering of signal at f(x(1),x(2))
% pos: if pos=1 make negative part zero

f=zeros(N,N);

dx=L/N;

Nr=3*N;			% for periodicity
R=zeros(Nr,Nr);
x=x+N;
jd=((1:Nr)-x(2));
for i=1:Nr
	vec=sqrt( (i-x(1))^2+jd.^2 );
	R(i,:)=vec;
end
R=c(R*dx);
if pos==1
	R(R<0)=0;
end

% add up to make periodic
ix=1:N;
iy=ix;
for i=1:3
	for j=1:3
		f=f+R(ix+(i-1)*N,iy+(j-1)*N);
	end
end

end




% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
