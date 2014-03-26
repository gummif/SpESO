% spectrum errors (slope and log norm)

timerID=tic;

SAVEDLOC='./';
path(path,genpath('functions'))

%LOADVAR=0; % load saved variables?
SAVEVAR=1; % save variables?

DWTTYPE=1;  % choose wavelet package
			% 1:Wavelab, 2:MATLAB's Wavelet toolbox, 0:cdf9/7

JH=0;
RHHDA=0;
SYNTH=0;

% =============data type===============
% JH data
%JH=1;  mstr='jh';
%N=1024;

% RHHDA hot-wire data
% RHHDA=1;  mstr='rhh';
% N = 1024*32;

% SYNTH data
SYNTH=1;  mstr='synth';
N = 1024*32;
synthtype=1;	 % 1:FFT,     2:DWT
syntslope=1;	 % 1:(3,5/3), 2:(5/3,3)
% =====================================


do1vec=8; 	% for SpESO (DECODER=82) use half:  82:(32,16,8), 1:(16,8,4)
SIMN=length(do1vec);
SIMN0=4;		% number of experiments for each case

EXTRA='';		% extra text for special runs, '' for nothing
K=284;			% filter length

Mvec=ceil((N+K-3)./do1vec);   % any M opt db
Mmin=Mvec(end);

% [1,31,312,32] , 82
for MATRIX=1      % 1=filter,  not used: 2=conv., 3=full
	for DECODER=82  %  1=LOMP, 31=M/2Best, 312=MBest, 32=shannon, 82=SpESO/QOMOMP
		
		FIGQUAL=sprintf('_%d_%d_%d_%s_%d_%d_%d_%s',N,Mmin,SIMN0,mstr,K,MATRIX,DECODER,EXTRA);
		SAVEFILE = [SAVEDLOC,'saved_1d',FIGQUAL,'.mat']; % for variables
		
		disp(' ')
		disp(FIGQUAL)
		
		spectra_error_1d_calc;  % Calculation
		
	end
end


toc(timerID)

disp('done')




% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
