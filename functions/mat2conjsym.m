function h = mat2conjsym(h)
% change a random matrix to complex conjugate symmetric matrix (for fft2)
% h has elements of form exp(i*x) only
% h:	complex square matrix (powers of 2 size)

N=size(h,1);

h([1,N/2+1],[1,N/2+1])=sign(imag(h([1,N/2+1],[1,N/2+1]))); %4 reals
h(N/2+2:N,1)=conj(h(N/2:-1:2,1)); % left col
h(1,N/2+2:N)=conj(h(1,N/2:-1:2)); % upper row
h(N/2+2:N,N/2+1)=conj(h(N/2:-1:2,N/2+1)); % middle col
h(2:N,N/2+2:N)=conj(rot90(h(2:N,2:N/2),2)); %other

end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
