% This code takes a 2D signal f of specified type, applies a sampling 
% or approximation strategy, and computes the energy spectrum. The goal
% is to estimate the energy spectrum well on a logarithmic scale.
% Depending on DECODER, the following methods are used. 
% (32)      -a fixed rate sampling and fourier interpolation
% (31,312)  -best M/2 or M-term wavelet approximation by thresholding
%           -compressive sampling/compressed sensing with linear sampling y=A*f, 
%            and reconstruction by
% (82)         -Spectrum Estimation by Sparse Optimization (SpESO) using
%               Quasi-Oracle Multilevel OMP (QOMOMP)
% The results can be saved into .mat files, and plotted with the plot scripts.
% 
% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
% See file README.md for more information
%

timerID=tic;

SAVEDLOC='./';
path(path,genpath('functions'))

%LOADVAR=0; % load saved variables?
SAVEVAR=1; % save variables?

DWTTYPE=1;  % choose wavelet package
			% 1:Wavelab, 2:MATLAB's Wavelet toolbox, 0:cdf9/7

JH=0;
SYNTH=0;

% =============data type===============
% JH data
% JH=1;  mstr='jh';
% Ns=1024;

% Synthetic data
Ns=1024;
SYNTH=1;  mstr='synth';
synthtype=1;	 % 1:FFT,		2:DWT
syntslope=1;	 % 1:(5/3),		2:(?)
% =====================================
VORT=0;  % use vorticity


N=Ns^2;
% undersampling ratio for each dimension (squared for total)
do1vec=4;	% for SpESO (DECODER=82) use half:  82:(11.31,5.66), 1:(8,4)
SIMN=length(do1vec);
SIMN0=2;		% number of experiments for each case

EXTRA='';		% extra text for special runs, '' for nothing
SEED=1;			% for random numbers

for MATRIX=2      %  2=conv., not used: 1=filter,3=full
	if MATRIX==1
		K=284;
		Msvec = ceil((Ns+K-3)./do1vec);   % db
		Mvec = Msvec.^2;
	elseif MATRIX==2
		K=1;
		Mvec = ceil(N./do1vec.^2);   % db 2ds
		Msvec = round(sqrt(Mvec));
		Mvec = Msvec.^2;
	end
	
	Msmin=Msvec(end);
	Mmin=Mvec(end);
	
	% [31,312,32] , 82
	for DECODER=[31,312,32]  %  31=M/2Best, 312=MBest, 32=shannon, 82=SpESO/QOMOMP
		
		FIGQUAL=sprintf('_%d_%d_%d_%s_%d_%d_%d_%d_%s',N,Mmin,SIMN0,mstr,K,MATRIX,DECODER,SEED,EXTRA);
		SAVEFILE = [SAVEDLOC,'saved_2d',FIGQUAL,'.mat']; % for variables
		
		disp(' ')
		disp(FIGQUAL)
		
		spectra_error_2d_calc;  % Calculation
		
	end
end

toc(timerID)

disp('done')



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
