function varargout = Dspec_function(Dfunc,iwt,y,CSPsi,CSPsit,N,J0,Lj,beta,alph,tol,tolf,q)
ya=size(y,2);

if ya==1
	if nargout==1
		x=Dfunc(y,CSPsi,CSPsit,N,J0,Lj,beta,alph,tol,tolf);
		varargout{1}=iwt(x);
	elseif nargout==2
		[x,outp]=Dfunc(y,CSPsi,CSPsit,N,J0,Lj,beta,alph,tol,tolf);
		varargout{1}=iwt(x);
		varargout{2}=outp;
	end
else
	if nargout==1
		x=Dfunc(y,CSPsi,CSPsit,N,J0,Lj,q,beta,alph,tol,tolf);
		varargout{1}=iwt(x);
	elseif nargout==2
		[x,outp]=Dfunc(y,CSPsi,CSPsit,N,J0,Lj,q,beta,alph,tol,tolf);
		varargout{1}=iwt(x);
		varargout{2}=outp;
	end
end

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
