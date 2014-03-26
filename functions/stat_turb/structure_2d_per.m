function S = structure_2d_per(x,p,rmax)
% calculates structure function S_p(r) for periodic x
% with r in direction (1,0)

% input:
%			x: 2d matrix N*N
%			p: order 
%			rmax: length of S (at most N-1)
% output: 
%			S: column vector of size rmax

[N,~]=size(x);
if rmax>=N
	rmax=N-1;
end
S=zeros(rmax,1);

for r=1:rmax
	%Smat=circshift(x,-r) - x;
	Smat=x([1+r:N,1:r],:) - x;
	Smat=abs(Smat);
	Smat=Smat.^p;
	%Smat=realpow(Smat,Pmat);
	
	
	S(r) = mean( Smat(:) );
	
	%S(r) = mean(mean( abs(circshift(x,-r) - x).^p ));
end

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
