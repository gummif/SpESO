function f = synthetic_signal_w(N,dim,param)
% Create periodic synthetic signal with multiaffine structure
% N:	signal length (N*1 in 1d or N*N in 2d)
% dim:	'1d' or '2d'
% param: struct with fields (e.g. param.k0=5)
% rand:	param.rand: (m,n) -> iid random matrix of size m*n
%		moments of rand control structure function exponent
% A random process for the construction of  multiaffine fields
% R. Benzi, L. Biferale, A. Crisanti, G. Paladin, M. Vergassola, A. Vulpiani
% Physica D  65,  352, 1993.

% f:	wavelet transform of synthetical signal (apply iwt to get signal)



if strcmp('1d',dim)==1
	
	f=zeros(N,1);
	J=log2(N);
	
	
	f(2)=1;
	for j=1:J-1
		ind=2^(j)+1:2^(j+1);
		ni=length(ind);
		fjm1=zeros(ni,1);
		fjm1(1:2:ni)=f(2^(j-1)+1:2^(j)); % parent scale j-1
		fjm1(2:2:ni)=f(2^(j-1)+1:2^(j));
		f(ind)=fjm1.*param.rand(ni,1).*sign(rand(ni,1)-0.5);
	end
	
	
elseif  strcmp('2d',dim)==1
	
	
	f=zeros(N,N);
	J=log2(N);
	
	f(1,1)=0;
	f(1,2)=1;
	f(2,1)=1;
	f(2,2)=1;
	for j=1:J-1
		
		for ww=1:3
			if ww==1
				ind = xyindex2x(1:2^j,2^(j)+1:2^(j+1),N,N);
			elseif ww==2
				ind = xyindex2x(2^(j)+1:2^(j+1),1:2^j,N,N);
			elseif ww==3
				ind = xyindex2x(2^(j)+1:2^(j+1),2^(j)+1:2^(j+1),N,N);
			end
			[ix,iy] = ind2sub([N N],ind);
			ixp=ceil(ix/2);
			iyp=ceil(iy/2);
			indp = sub2ind([N N],ixp,iyp);
			
			f(ind)=f(indp)'.*param.rand(2^(2*j),1).*sign(rand(2^(2*j),1)-0.5);
			if j>=J-1-3 && j<J-2 && ww==3
				f(ind)=f(ind)/1.15;
			end
		end
		
	end
	
	
else
	disp('error: dim 1d or 2d')
	return;
end




end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
