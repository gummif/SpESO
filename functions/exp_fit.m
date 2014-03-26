function sc = exp_fit(f,nvec,scale,x0,nsc)
% f: vector to fit to exponential curve to
% nvec: vector of locations for f (or k for freq.)
% scale: function of scale e.g. x, log(x), sqrt(x)
% x0: inital guess vector
% nsc: scaling of nvec, 'log' for semixlog scaling
% Computes exponent and constant for fit c*n^(-s)
% sc(1): s (exponent)
% sc(2): c (constant)
f=f(:);

if nargin<4
	x0=[1,f(1)];
end
if nargin<5
	nsc='';
end

N=length(f);
%nvec=(round(q*N)+1:round(N*(1-q)))';
f=scale(f);

if strcmp(nsc,'log')==1
	func=@(scx)comperrlog(f,nvec,scx(1),scx(2),scale);
else
	func=@(scx)comperr(f,nvec,scx(1),scx(2),scale);
end
op.MaxFunEvals=1100;
sc = fminsearch(func,x0,op);

end

function e=comperr(fhat,nvec,s,c,scale)
e=norm(fhat-scale(c*nvec.^(-s)));
end

function e=comperrlog(fhat,nvec,s,c,scale)
e=norm( (fhat-scale(c*nvec.^(-s)))./sqrt(nvec) );
end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
