% spectrum errors (slope and log norm)

timerID=tic;

SAVEDLOC='./saved_var/';
path(path,genpath('functions'))

LOADVAR=0; % load saved variables?
SAVEVAR=1; % save variables?

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
do1vec=5.66; %increasing Ns/Ms.    4/5.66, 8/11.31
SIMN=length(do1vec);


EXTRA='f1'; % extra text for special runs, '' for nothing
% w1
% vor,vorticity


SIMN0=8;		% number of experiments for each case
SEED=1;			% for random numbers
% ============================


for MATRIX=2      % 1=filter, 2=conv., 3=full
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
	for DECODER=82  % 0=BP, (1=LOMP), 31=M/2Best,312=MBest, 32=shannon, 82=QOMspec
		
		FIGQUAL=sprintf('_%d_%d_%d_%s_%d_%d_%d_%d_%s',N,Mmin,SIMN0,mstr,K,MATRIX,DECODER,SEED,EXTRA);
		SAVEFILE = [SAVEDLOC,'spectra_error_2d_saved',FIGQUAL,'.mat']; % for variables
		
		disp(' ')
		disp(FIGQUAL)
		
		spectra_error_2d_calc;  % Calculation
		
		%slope_fun_m2d_plot1;  % images
		%slope_fun_m2d_plot2;  % spectra
		
	end
end

toc(timerID)

disp('EXITING')
exit



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
