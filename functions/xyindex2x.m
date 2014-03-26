function x1=xyindex2x(xv,yv,Nx,Ny)
% takes xv and yv index vector for A(xv,yv) and makes x1 such that A(x1)
% index the same values. A(xv,yv) means the array of values by any
% combination of indexes from xv and yv. A(1:2,1:2) is therefore a 2*2
% cube, with length(x1)=4.
% xv,yv: 1D index vectors
% N: size of A N*N
% x1: 1d index vector combined of xv and yv

Lx=length(xv);
Ly=length(yv);

x1=zeros(1,Lx*Ly);

%for i=1:length(xv)
	for j=1:Ly
		x1((j-1)*Lx+1 : j*Lx)=sub2ind([Nx Ny], xv, ones(1,Lx)*yv(j));
	end
%end



end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
