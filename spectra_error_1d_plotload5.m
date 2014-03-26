
%========== signals =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-rgb';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
FIGQUAL='';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;

COMPTYPEv=[1];

for COMPTYPE=COMPTYPEv
	
	EXTRASTR='';
	if COMPTYPE==1					 % '64', [32,16,8,4], [1,2,3,4]
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='64';  %8,16,64
		Mratio=4; %16,8,4
		RUNTYPEv=[1,2,3,4,10];
	end

	
	% 1=w353, 2=w533 wavelet
	% 3=f353, 4=f533 fourier
	% 10=rhh
	%RUNTYPEv=10;		%1,2,3,4,5
	for RUNTYPE = RUNTYPEv
		
		if Mratio==4
			AVM='_4';
			M1=num2str(8263);
			M2=num2str(4132);
		end
		if RUNTYPE==1
			RUNSTR=['_w353',EXTRASTR];
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31',RUNSTR,'.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1',RUNSTR,'.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);
		elseif RUNTYPE==2
			RUNSTR=['_w533',EXTRASTR];
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31',RUNSTR,'.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1',RUNSTR,'.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);
		elseif RUNTYPE==3
			RUNSTR=['_f353',EXTRASTR];
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31',RUNSTR,'.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1',RUNSTR,'.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);  %d2av2d10
		elseif RUNTYPE==4
			RUNSTR=['_f533',EXTRASTR];  %s4 (av2),s2
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31',RUNSTR,'.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1',RUNSTR,'.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);  %d2av2d10
		elseif RUNTYPE==10
			RUNSTR=['_',EXTRASTR];
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_rhh_284_1_31',RUNSTR,'.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_rhh_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_rhh_284_1_32',RUNSTR,'.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_rhh_284_1_1',RUNSTR,'.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_rhh_284_1_82',RUNSTR,'d2if.mat']);
			RUNSTR=['_rhh',EXTRASTR];
		end
		
		if RUNTYPE==1
			ylimit=[8e-16 4e-10];
			C1=1;C2=1;p1=-3;p2=-5/3;
		elseif RUNTYPE==2
			ylimit=[8e-16 4e-10]*1000;
			C1=1;C2=1e8;p1=-5/3;p2=-3;
		elseif RUNTYPE==3
			ylimit=[8e-16 4e-10]*100;
			C1=60;C2=60;p1=-3;p2=-5/3;
		elseif RUNTYPE==4
			ylimit=[8e-16 4e-10]*100000;
			C1=60;C2=60e8;p1=-5/3;p2=-3;
		elseif RUNTYPE==10
			ylimit=[8e-12 4e-6]*100000;
			C1=60;C2=60e8;p1=-5/3;p2=-3;
		end
		
		N=length(var1.f_os{1});
		J=log2(N);
		kvec=var1.k_os{1};
		[SIMN0,~]=size(var1.f_os);
		
		f=var1.f_os{1};
		
		ind=N-2200+1:N-100;
		
		figure;
		set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);

		hh=plot(0:length(ind)-1,f(ind),'-','Color',[0 0 0],'LineWidth',0.8);
		set(gca,'YTick',[])

		axis tight
		% \hspace{-15mm} $ \hspace{10mm}
% 		if ~(COMPTYPE==2 || COMPTYPE==3)
% 			ylabel('\hspace{-2mm} $ E(k)$','FontSize', 11);
% 		end
% 		if ~((COMPTYPE==1||COMPTYPE==2||COMPTYPE==3) && (RUNTYPE==1||RUNTYPE==2||RUNTYPE==3))
% 			xlabel('$k$','FontSize', 11);
% 			set(gca,'XTick',[10^2,10^3,10^4])
% 		end
		%set(gca,'Outerposition',[-0.08,0,1.1,1])  %posx,posy,widx,widy
		%set(gca,'XTick',jvec)
		
%		ylim(ylimit)
	%	xlim([50 N/2])
		
		
		%legend([hh0,hh1,hh1b,hh2,hh3,hh4],{'original','$M$-best','$M/2$-best','Shannon','LOMP',CSNAME}...
	%		,'interpreter', 'latex','Location','southwest');
		
		
		if SAVEFIG==1
			export_fig([FIGLOC,'u1',FIGQUAL,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'u1',FIGQUAL,RUNSTR], '-eps', '-transparent', FIGCOLOR)
			disp([FIGLOC,'u1',FIGQUAL,RUNSTR])
		end
		
		
	end
	
end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
