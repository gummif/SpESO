function L = Lvec_QOMOMP(N,M,dim,type,extra)
% creata an L for QOM Orthogonal Matching Pursuit (OMP)
% for dim='2d' N and M are side lengths
% L: number of coefs to approximate, L(j), at level j

J=log2(N);
J_M=round(floor(log2(M/2)));

if strcmp('1d',dim)==1
	L = pow2(1:J-1)';
	
	if type==1
		mvec(J_M-1:J_M+4)=[0.97,0.38,2/32,2.4/128,2.0/512,0.7/4096];
		mvec(J_M+5:J-1)=1/4096/4;
		
		for j=J_M-1:J-1
			L(j)=L(j)*mvec(j);
		end
		
		if nargin<5
			extra=1;
		end
		% extra scaling
		for j=1:J-1
			a=L(j)/2^j;
			C=a/(1-a);
			a2=extra(1)*a/(1-a*(1-extra(1)));
			L(j)=2^j*a2;
		end
	elseif type==2
		mvec(J_M-1:J_M+4)=[0.95,0.4,1.8/32,4/128,1.5/512,0.5/4096];
		mvec(J_M+5:J-1)=1/4096/4;
		
		for j=J_M-1:J-1
			L(j)=L(j)*mvec(j);
		end
		
		if nargin<5
			extra=1;
		end
		% extra scaling
		for j=1:J-1
			a=L(j)/2^j;
			C=a/(1-a);
			a2=extra(1)*a/(1-a*(1-extra(1)));
			L(j)=2^j*a2;
		end
		
	end
	
	L(L<0)=0; % negative Lj = 0
	for j=1:J-1 % larger than 2^j Lj = 2^j
		if L(j)>2^j
			L(j)=2^j;
		end
	end
	
elseif  strcmp('2d',dim)==1
	L = zeros(J-1,1);
	L(1:J-1)=pow2(2*(1:J-1))*3;
	
	if type==1
		L(J_M) = pow2(2*(J_M)-0.6)*3;
		L(J_M+1) = pow2(2*(J_M)-2)*3;
		L(J_M+2) = pow2(2*(J_M)-3.6)*3;
		L(J_M+3) = pow2(2*(J_M)-6)*3;
		%L(J_M+2:end) = 0;
	end
	
	L(L<0)=0; % negative Lj = 0
	for j=1:J-1 % larger than max
		if L(j)>pow2(2*j)*3
			L(j)=pow2(2*j)*3;
		end
	end
end


L=round(L);

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
