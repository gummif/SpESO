function p = synthetic_param_w(N,dim,type,extra)
% Parameters synthetic signal synthetic_signal(N,dim,param)

if strcmp('1d',dim)==1
	if type==1
		synpr=extra(1);
		p=struct('rand',@(m,n) eta_dist1(m,n,0.125*synpr,pow2(-1/2),pow2(-5/6)));
	elseif type==2
		e0=extra(1);
		e1=extra(2);
		% choose e0<5/6 and e1>5/6 for -5/3 spectrum
		x0=2^(-e0);
		x1=2^(-e1);
		espec=extra(3);
		p0=(pow2(-espec)-pow2(-e1*2))/(pow2(-e0*2)-pow2(-e1*2));
		p=struct('rand',@(m,n) eta_dist1(m,n,p0,x0,x1));
		
	end
	
elseif strcmp('2d',dim)==1
	if type==1
		pvec=2:4;
		mom=pow2(-pvec*(1/3+2/2));
		
		
		
		
		p0=1/4;
		valop=fminsearch(@(x) norm(mom-mom_dist1(pvec,p0,x(1),x(2))),[0.4,0.3]);
		%p0=valop(1);
		x0=valop(1);
		x1=valop(2);
		
		momop=mom_dist1(pvec,p0,x0,x1);
		%-log2(momop)-2/2*pvec
		%pvec/3

		p=struct('rand',@(m,n) eta_dist1(m,n,p0,x0,x1));
	elseif type==2
		e0=extra(1);
		e1=extra(2);
		% choose e0<5/6+1/2 and e1>5/6+1/2 for -5/3 spectrum
		espec=extra(3);
		x0=2^(-e0);
		x1=2^(-e1);
		p0=(pow2(-espec-1)-pow2(-e1*2))/(pow2(-e0*2)-pow2(-e1*2));
		p=struct('rand',@(m,n) eta_dist1(m,n,p0,x0,x1));
		
	elseif type==3
		e0=4/6+1/2;
		e1=5/6+1/2;
		e2=7/6+1/2;
		x0=2^(-e0);
		x1=2^(-e1);
		x2=2^(-e2);
		%p0=(pow2(-5/6*2-1)-pow2(-e1*2))/(pow2(-e0*2)-pow2(-e1*2));
		p1=1/3;
		p0=1/3*1.1;
		p=struct('rand',@(m,n) eta_dist2(m,n,p0,p1,x0,x1,x2));
		
	end
	
end

end

function eta=eta_dist1(m,n,p0,n0,n1)

X=rand(m,n);
eta=zeros(m,n);
eta(X<p0)=n0;
eta(X>=p0)=n1;

end
function m=mom_dist1(n,p0,n0,n1)
m = n0.^n*p0 + n1.^n*(1-p0);
end


function eta=eta_dist2(m,n,p0,p1,n0,n1,n2)

X=rand(m,n);
eta=zeros(m,n);
eta(X<p0)=n0;
eta(X>=p0 & X<p1+p0)=n1;
eta(X>=p1+p0)=n2;

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
