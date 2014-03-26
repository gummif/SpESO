function f=waveletvec2mat(x,s)
% convert 2d wavelet toolbox vector to matrix form
% x: vector
% s: information matrix
% f: square matrix

f=zeros(s(end,:));
J=size(s,1)-1;

nj=s(1);
i1=1; i2=nj^2;
f(1:nj,1:nj)=reshape(x(i1:i2),[nj,nj]);

for j=2:J
	nj=s(j); Lj=sum(s(1:j-1,2));

	i1=i2+1; i2=i1+nj^2-1;
    h=x(i1:i2); 
    f(1:nj,(1:nj)+Lj)=reshape(h,[nj,nj]);
	
    i1=i2+1; i2=i1+nj^2-1;
    h=x(i1:i2); 
    f((1:nj)+Lj,1:nj)=reshape(h,[nj,nj]);
	
    i1=i2+1; i2=i1+nj^2-1;
    h=x(i1:i2); 
    f((1:nj)+Lj,(1:nj)+Lj)=reshape(h,[nj,nj]);
	
end




end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
