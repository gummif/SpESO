function w = window_spec(N,type)
% returns a window w for signal processing
% N: lenght of window
% type: type of window (hanning,hamming,box,welch,blackman)
% smoothness order: blackman>hanning>welch>hamming>box

if strcmp(type,'hanning')
	
	w=0.5-0.5*cos(2*pi*(0:N-1)./(N-1));
	
elseif strcmp(type,'hamming')
		
	w=0.54-0.46*cos(2*pi*(0:N-1)./(N-1));
		
elseif strcmp(type,'box')
		
	w=ones(1,N);
	
elseif strcmp(type,'welch') % recommended by NumRec
		
	w=1-(((0:N-1)-0.5*N)./(0.5*N)).^2;
	
elseif strcmp(type,'blackman') % alpha=0.16
	a=0.16;	
	w=(1-a)/2-1/2*cos(2*pi*(0:N-1)./(N-1))+a/2*cos(4*pi*(0:N-1)./(N-1));
	
end

w=w';

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
