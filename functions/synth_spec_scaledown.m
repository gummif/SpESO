function f=synth_spec_scaledown(f,k0,c)
% scale 2d fourier coefs of f from |k|=k0 to end by amount c on a log scale

[m,n]=size(f);

if m==1 || n==1  %1d
	
	fhat=fft(f);
	N=length(f);
	
	ind=k0:N/2+1;
	amp=logspace(0,log10(c),length(ind));
	fhat(ind)=fhat(ind)./amp';
	fhat(N-ind+2)=fhat(N-ind+2)./amp';
	
	f=real(ifft(fhat));
	
else  %2d
	
	fhat=fft2(f);
	N=length(f);
	
	kivec=zeros(1,N);
	for i=1:N
		kivec(i)=fft2k(N,i);  % k_1 , actual frequencies
	end
	
	Kmat=zeros(N,N);
	for j=1:N
		Kmat(:,j)=round(sqrt( kivec.^2+kivec(j)^2 ))+1; % freq. index for E
	end
	
	k_ind=k0:floor(N/sqrt(2));
	amp=logspace(0,log10(c),length(k_ind));
	for i=1:length(k_ind)
		kk=k_ind(i);
		% 	if kk>N/2*0.95;
		% 		amp(i)=10000;
		% 	end
		ind=find(Kmat==kk);
		fhat(ind)=fhat(ind)/amp(i);
	end
	f=real(ifft2(fhat));
	
end

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
