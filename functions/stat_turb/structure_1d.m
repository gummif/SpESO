function S = structure_1d(x,p)
% calculates structure function S_p(r) for x of size M*1
% at position x(1)
% columns are measurements in space
% input:
%			x: vector of values
%			r: r>0, separation x(i+r)-x(i)
%			p: order 
% output: 
%			S: column vector of size M-1. The function S_p(r) at x(1).

[M,~]=size(x);
S=zeros(M-1,1);

S=(x(1+1:M)-x(1)).^p;



end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
