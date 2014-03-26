function xmin = spectrum_gradopt(x0,Psi,At,D,E,g,param)
% Minimizes spectrum error in log-space by varying L in QOMOMP
% x0:	initial guess (in transformed space(wavelet))
% A:	measurement matrix operator
% At:	measurement matrix operator transpose
% D:	QOMOMP solver function of L, f_app=D(L)
% E:	spectrum function of a signal, [kvec,Evec]=E(f)
% g:	the measurement vector g=Psi \hat f  (or g=Af)
% param:	parameters struct for optimization
%		tolf: final tolerance for QOMOMP
%		d: linspace bounds ( exp(-d,...,d) )
%		df: final variance (sigma^2 for log-normal distribution)
%		iter: number of iterations
%		jvec: vector of scales to optimize
%		res: resolution of evaluations for each scale
%		resf: resolution of final average
% Lopt:	optimized L vector
% Eopt:	optimized spectrum

tol=param.tol;
mu=param.mu;
h=param.h;
tmax=param.tmax;
B=param.B;
N=length(x0);

Eu=E(At(g));
xt=x0*99999999;
errmin=norm(log(Eu./E(At(Psi(x0)))));
xmin=x0;
gradf=x0*0;
f=@(x) norm(log(Eu./E(At(Psi(x)))));

t=1;
while norm(x0-xt)>tol && t<tmax
	
	xih=xt;
	for i=1:N
		xih(i)=xih(i)+h;
		gradf(i)=f(xih);
		xih(i)=xih(i)-2*h;
		gradf(i)=(gradf(i)-f(xih))/(2*h);
		xih(i)=xt(i);
		i
	end
	xt=xt-mu*gradf;
	xt=filter_coef(xt,B);
	
	if errt<errmin
		xmin=xt;
		errmin=errt;
	end

	t=t+1;
	disp(t)
end
disp(t)




end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
