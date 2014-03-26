function f=waveletmat2vec(x,s)
% convert matrix form to 2d wavelet toolbox vector
% x: square matrix
% s: information matrix
% f: vector

f=zeros(numel(x),1);
J=size(s,1)-1;

nj=s(1);
i1=1; i2=nj^2;
h=x(1:nj,1:nj);
f(i1:i2)=h(:);

for j=2:J
	nj=s(j); Lj=sum(s(1:j-1,2));

	i1=i2+1; i2=i1+nj^2-1;
    h=x(1:nj,(1:nj)+Lj); 
    f(i1:i2)=h(:);
	
    i1=i2+1; i2=i1+nj^2-1;
    h=x((1:nj)+Lj,1:nj); 
    f(i1:i2)=h(:);
	
    i1=i2+1; i2=i1+nj^2-1;
    h=x((1:nj)+Lj,(1:nj)+Lj); 
    f(i1:i2)=h(:);
	
end




end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
