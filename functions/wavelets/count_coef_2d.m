function [coefnum,nzr] = count_coef_2d(fM)
% count number of non-zero coefficients of each scale of the 2d wavelet transform
% f: is a matrix
% coefnum: vector of coefs at scale j
% nzr: non-zero ratio at scale j between coefnum and wavelet type 3 coefs
% 
% example: fM = filter_coef(fwt(f),M);

[N,M]=size(fM);
J=log2(N);

for j=1:J-1;
	qw1=fM(1:2^j,2^j+1:2^(j+1));
	qw2=fM(2^j+1:2^(j+1),1:2^j);
	qw3=fM(2^j+1:2^(j+1),2^j+1:2^(j+1));
	nnz(qw1);
	nnz(qw2);
	nnz(qw3);
	coefnum(j)=nnz(qw3)+nnz(qw2)+nnz(qw1);
	nzr(j)=nnz(qw3)/coefnum(j);
end

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
