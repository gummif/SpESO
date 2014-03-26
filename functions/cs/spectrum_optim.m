function [Lopt,Eopt] = spectrum_optim(L0,A,At,D,E,g,param)
% Minimizes spectrum error in log-space by varying L in QOMOMP
% L0:	initial guess
% A:	measurement matrix operator
% At:	measurement matrix operator transpose
% D:	QOMOMP solver function of L, f_app=D(L)
% E:	spectrum function of a signal, [kvec,Evec]=E(f)
% g:	the measurement vector g=Af
% param:	parameters struct for optimization
%		tolf: final tolerance for QOMOMP
%		d: linspace bounds ( exp(-d,...,d) )
%		df: final variance (sigma^2 for log-normal distribution)
%		iter: number of iterations
%		jvec: vector of scales to optimize
%		res: resolution of evaluations for each scale
%		resf: resolution of final average
%		dimsc: restrict L(j)<L(j-1)*dimsc (default==1)
% Lopt:	optimized L vector
% Eopt:	optimized spectrum

try dimsc=param.dimsc; 
catch ME, dimsc=1; end

L=L0;
[kvec,Eg]=E(At(g));

jvec=param.jvec;
tolf=param.tolf;
d=param.d;
Res=param.res - mod(param.res,2)+1; %make sure it's odd
Ldist=exp(linspace(-d,d,Res));




% for i=1:param.resf
% 	Li=L;
% 	Li(jvec)=Li(jvec).*exp(randn(length(jvec),1)*param.df);
% 	Li=round(Li);
% 	% QOMOMP
% 	u=D(Li,tolf);
% 	% spectrum
% 	[~,Ea]=E(u);
%
% 	Eopt=Eopt+log(Ea);
%
% end


tolf=tolf*10;

Ljend=L(jvec(end));
L(jvec(end))=0;  % was not

for s=1:param.iter			% iterations
	if s==2
		tolf=tolf/10;
	end
	if s>=2
		d=d*0.6;
		Ldist=exp(linspace(-d,d,Res));
	end
	%figure
	for j=jvec		% scales
		wj=find(kvec>2^(j-1) & kvec<=2^j);
		if s==3
			wj=find(kvec>2^(jvec(1)-1) & kvec<=2^jvec(end)); %all !!!
		end
		Lcoef=zeros(Res,1);
		e=Lcoef;
		if s==1 && j==jvec(end) 
			L(j)=Ljend;
		end
		
		logspec=log(Eg(wj));
		%logspec=(Eg(wj));  % linear
		for i=1:Res	% resolution
			Li=L;
			Li(j)=Li(j)*Ldist(i);
			Li=round(Li);
			Lcoef(i)=Li(j);
			% QOMOMP
			[u,outp]=D(Li,tolf);
			% spectrum
			[~,Ea]=E(At(A(u)));
			% error in log space
			flag = outp.flag ~= 0;  % 0 if 0, 1 if nonzero
			%confact=((1-flag)+flag*outp.relres/tolf); % convergence factor
			confact=1;  % do nothing
			e(i)=norm(log(Ea(wj))-logspec)*confact;
			%e(i)=norm(Ea(wj)-logspec)*confact;  % linear
			
			flag(i)=outp.flag;
			
			fprintf('i');
			
		end
		%plot(Lcoef,e,'-*'),hold all
		
		if L0(j)<L0(j-1)
			eind=Lcoef<=L(j-1)*dimsc;  % restrict L(j)<=L(j-1)*dimsc
		else
			eind=Lcoef<=inf;
		end
		imin=find(e(eind)==min(e(eind)),1,'first');
		L(j)=Lcoef(imin);
		
		fprintf(' j %d ,',L(j));
		disp(flag(:)')
	end
	fprintf('\n');
end
%figure,plot(Lcoef,e,'-*')

disp('final optim')
Eopt=Eg*0;
% average spectra
for i=1:param.resf
	Li=L;
	%Li(jvec)=Li(jvec).*exp(randn(length(jvec),1)*param.df);  % randomize?
	Li=round(Li);
	% QOMOMP
	%u=D(Li,tolf);
	[u,outp]=D(Li,tolf);
	disp(outp)
	% spectrum
	[~,Ea]=E(u);
	
	Eopt=Eopt+log(Ea);
	
end

Eopt=exp(Eopt/param.resf);
Lopt=round(L);






end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
