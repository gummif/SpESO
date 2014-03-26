function A = func2mat(Afunc,N,M)
% func2mat(Afunc,N,M) converts matrix function Afunc to M*N matrix

A=zeros(M,N);
t=zeros(N,1);
for i=1:N
	t(i)=1;
	A(:,i)=Afunc(t);
	t(i)=0;
end


end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
