function varargout = QOMOMPt_2d_gfa(y,Psi,Psit,N,J0,L,q,beta,alph,tol,tolf)
% QOM Orthogonal Matching Pursuit (OMP) 2D to solve a CS problem
% y:	A*f=Psi*fhat, CS measurements (2D matrix)
% Psi:	CS matrix function, Psit the transpose
% N:	number of elements in f Ns^2=N
% J0:	directly approximate coefs (1:2^J0,1:2^J0)
% L:	number of coefs to approximate, L(j), at level j
% q:	q(j) ratio of coefs d^3/d^total for scale j, negative q regular
% beta: tree structure multiplicative coef (>= 1)
% alph: function determining the beta multiplication, alph: R^|omeg| -> R
% tol:	tolerence for the least squares
% tolf: tolerence for the final low-tol least squares
% varargout:
% {1} x: approximate solution to y=Psi*fhat (2D matrix)
% {2} outp: information about solution convergence

Ns=sqrt(N);

r=y;
psyt=Psit(y);
x=zeros(Ns,Ns);
ind=1:pow2(J0);
ind=xyindex2x(ind,ind,Ns,Ns);
if (nargin < 7), tol=1e-4; end  % default is 1e-6 for lsqr
if (nargin < 8), tolf=1e-9; end % final iteration
L(L<0)=0; % negative Lj = 0
for j=1:log2(Ns)-1 % larger than maximum
	if L(j)>pow2(j*2)*3
		L(j)=pow2(j*2)*3;
	end
end
	
x0=x(ind)';
b=psyt(ind); b=b(:);
[xt,flag] = symmlq(@(xx) itfunc(xx,Psi,Psit,ind,Ns),b,tol,200,[],[],x0);

% solution of oracle
x(ind)=xt;
r=y-Psi(x);


for j=J0:log2(Ns)-1

	itdiff = xyindex2x(1:2^(j-1),1:2^(j-1),Ns,Ns); %to remove
	omeg = setdiff(ind,itdiff,'stable'); % Omega_(j-1)
	omegs=find(abs(x)>alph(x(omeg))); 
	omegs=intersect(omegs,omeg);% Omega_(j-1)^star
	[ix,iy] = ind2sub([Ns Ns],omegs); % change to sub
	ix=[ix*2,ix*2-1,ix*2,  ix*2-1];
	iy=[iy*2,iy*2,  iy*2-1,iy*2-1];
	omeg = sub2ind([Ns Ns],ix,iy); % Omega_(j)
	
	% index set
	temp=abs(Psit(r)); 
	temp(omeg)=temp(omeg)*beta(j);
	temp2=temp(1:2^(j+1),1:2^(j+1));

	
	%============

	
	[~,i] = sort(temp2(:));
	itdiff = xyindex2x(1:2^j,1:2^j,pow2(j+1),pow2(j+1)); %to remove
	i = setdiff(i,itdiff,'stable'); %remove itdiff
	if q(j)<0
		[ix,iy] = ind2sub([pow2(j+1) pow2(j+1)],i); % change to sub
		i = sub2ind([Ns Ns],ix,iy); % change to index in full matrix
		i=unique(i,'stable');
		
		i=i(end-L(j)+1:end);
		ind=union(ind,i);
	else
		itdiff = xyindex2x(2^j+1:2^(j+1),2^j+1:2^(j+1),pow2(j+1),pow2(j+1)); %to remove
		i12 = setdiff(i,itdiff,'stable'); %remove itdiff
		i3 = setdiff(i,i12,'stable'); %remove i12
		
		[ix12,iy12] = ind2sub([pow2(j+1) pow2(j+1)],i12); % change to sub
		i12 = sub2ind([Ns Ns],ix12,iy12); % change to index in full matrix
		i12=unique(i12,'stable');
		
		if round(L(j)*(1-q(j)))>length(i12) 
			q(j)=1-length(i12)/L(j);  % Lj overrules q
		end
		i12=i12(end-round(L(j)*(1-q(j)))+1:end);
		ind=union(ind,i12);
		
		[ix3,iy3] = ind2sub([pow2(j+1) pow2(j+1)],i3); % change to sub
		i3 = sub2ind([Ns Ns],ix3,iy3); % change to index in full matrix
		i3=unique(i3,'stable');
		
		i3=i3(end-round(L(j)*q(j))+1:end);
		ind=union(ind,i3);
	end
	ind=sort(ind);
	
	% iterative method (with pseudo-inverse LS) no transpose
	% (bicgstab, bicgstabl, cgs, gmres, minres, pcg, symmlq, tfqmr)
	x0=x(ind); x0=x0(:);
	b=psyt(ind); b=b(:);
	
	[xt,~]  = symmlq(@(xx) itfunc(xx,Psi,Psit,ind,Ns),b,tol,200,[],[],x0);

	% update
	x(:)=0;
	x(ind)=xt;
	r=y-Psi(x);
	

end

% one final low tolerance iteration
[xt,flag,relres,iter] = symmlq(@(xx) itfunc(xx,Psi,Psit,ind,Ns),b,tolf,200,[],[],x0);
x(:)=0;
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
function val=itfunc(x,Psi,Psit,ind,Ns)

x2=zeros(Ns,Ns);
x2(ind)=x;

val=Psit(Psi(x2));
val=val(ind);  val=val(:);

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
