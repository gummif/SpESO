function varargout = QOMOMPt_gfa(y,Psi,Psit,N,J0,L,beta,alph,tol,tolf)
% QOM Orthogonal Matching Pursuit (OMP) to solve a CS problem (tree)
% y=A*f=Psi*fhat, CS measurements
% Psi: CS matrix function, Psit the transpose
% N: length of f
% J0: directly approximate coefs 1:2^J0
% L: number of coefs to approximate, L(j), at level j
% beta: tree structure multiplicative coef (>= 1)
% alph: function determining the beta multiplication, alph: R^|omeg| -> R
% tol: tolerence for the least squares
% tolf: tolerence for the final low-tol least squares
%%% d:	diagonal for preconditioning
% varargout:
% {1} x: approximate solution to y=Psi*fhat
% {2} outp: information about solution convergence

r=y;
psyt=Psit(y);
x=zeros(N,1);
ind=1:2^J0;
%if (nargin < 6), L=1; end
if (nargin < 8+1), tol=1e-4; end  % default is 1e-6 for lsqr
if (nargin < 9+1), tolf=1e-9; end % final iteration
L(L<0)=0; % negative Lj = 0
L=round(L);
for j=1:log2(N)-1 % larger than 2^j Lj = 2^j
	if L(j)>2^j
		L(j)=2^j;
	end
end

x0=x(ind);
b=psyt(ind);


%b=sign(b).*(abs(b)).^(0.8);
[xt,~] = tfqmr(@(xx) itfunc(xx,Psi,Psit,ind,N),b,tol,200,[],[],x0);
%xt = nsoli(x0,@(xx) nlfunc(xx,Psi,Psit,b,ind,N),[tol,0]);


% solution of oracle
x(ind)=xt;
r=y-Psi(x);


for j=J0:log2(N)-1
	% index set
	omeg=intersect(ind,2^(j-1)+1:2^(j)); % Omega_(j-1)
	omegs=find(abs(x)>alph(x(omeg)));
	omegs=intersect(omegs,omeg);% Omega_(j-1)^star
	omeg=[omegs*2,omegs*2-1]; % Omega_(j)
	
	
	temp=abs(Psit(r)); %figure;plot(x),figure;plot(temp)
	temp(omeg)=temp(omeg)*beta(j);
	[~,i] = sort(temp(2^j+1:2^(j+1)));
	i=i(end-L(j)+1:end) + 2^j;
	ind=union(ind,i);
	
	% figure,plot_coef_im(abs(x)>0,1:length(x))
	
	% iterative method (with pseudo-inverse LS) no transpose
	% (bicgstab, bicgstabl, cgs, gmres, minres, pcg, symmlq, tfqmr)
	x0=x(ind);
	b=psyt(ind);
	
	
	%b=sign(b).*(abs(b)).^(0.8);
	[xt,~]  = tfqmr(@(xx) itfunc(xx,Psi,Psit,ind,N),b,tol,200,[],[],x0);
	%xt = nsoli(x0,@(xx) nlfunc(xx,Psi,Psit,b,ind,N),[tol,0]);
	
	% update
	x(1:end)=0;
	x(ind)=xt;
	r=y-Psi(x);
	
	%plot(x),xlim([-5,100])
end

% one final low tolerance iteration
%[xt,~] = tfqmr(@(xx) itfunc(xx,Psi,Psit,ind,N),b,tolf,200,[],[],x0);
[xt,flag,relres,iter] = tfqmr(@(xx) itfunc(xx,Psi,Psit,ind,N),b,tolf,200,[],[],x0);
%xt = nsoli(x0,@(xx) nlfunc(xx,Psi,Psit,b,ind,N),[tolf,0]);

%fprintf('flag=%d, iter=%d\n',flag,iter);


x(1:end)=0;
x(ind)=xt;

if nargout==1
	varargout{1}=x;
elseif nargout==2
	varargout{1}=x;
	outp.flag=flag;
	outp.relres=relres;
	outp.iter=iter;
	varargout{2}=outp;
end
end

% without transposes, val=A'*A*x, where A=Psi(:,ind)
function val=itfunc(x,Psi,Psit,ind,N)

x2=zeros(N,1);
x2(ind)=x;

val=Psit(Psi(x2));
val=val(ind);

%val=sign(val).*(abs(val)).^(0.8);

end

% dd
function val=nlfunc(x,Psi,Psit,b,ind,N)

x2=zeros(N,1);
x2(ind)=x;

val=Psit(Psi(x2));
val=val(ind)-b;

%val=val(ind);
%val=sign(val).*log(abs(val))-sign(b).*log(abs(b));
val=sign(val).*(abs(val)).^(0.5);

%val(abs(val)>1e9)=-1e9;

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
