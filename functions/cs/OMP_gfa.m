function x = OMP_gfa(y,Psi,Psit,N,B,L,tol,tolf)
% Orthogonal Matching Pursuit (OMP) to solve a CS problem
% y=A*f=Psi*fhat, CS measurements
% Psi: CS matrix function, Psit the transpose
% N: length of f
% B: sparsity of fhat (number of iterations * L)
% L: number of lumped coefficients in each estimation (5-20 a good number)
% tol: tolerence for the least squares
% tolf: tolerence for the final low-tol least squares
% x: approximate solution to y=Psi*fhat

r=y;
psyt=Psit(y);
x=zeros(N,1);
ind=[];
if (nargin < 6), L=1; end
if (nargin < 7), tol=1e-4; end  % default is 1e-6 for lsqr
if (nargin < 8), tolf=1e-10; end % final iteration

for t=1:round(B/L)
	% index set
	[~,i] = sort(abs(Psit(r)));
	i=i(end-L+1:end);
	ind=union(ind,i);
	
	% iterative method (with pseudo-inverse LS) no transpose 
	% (bicgstab, bicgstabl, cgs, gmres, minres, pcg, symmlq, tfqmr)
	x0=x(ind);
	b=psyt(ind);

	[xt,~]  = symmlq(@(xx) itfunc(xx,Psi,Psit,ind,N),b,tol,200,[],[],x0);
%[xt,~] 
	
	% update
	x(1:end)=0;
	x(ind)=xt;
	r=y-Psi(x);
	
%plot(x),xlim([-5,100])
end

% one final low tolerance iteration
[xt,flag] = symmlq(@(xx) itfunc(xx,Psi,Psit,ind,N),b,tolf,200,[],[],x0);
x(1:end)=0;
x(ind)=xt;

end

% without transposes, val=A'*A*x, where A=Psi(:,ind)
function val=itfunc(x,Psi,Psit,ind,N)

x2=zeros(N,1);
x2(ind)=x;

val=Psit(Psi(x2));
val=val(ind);

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
