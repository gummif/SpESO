
%========== structure function =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-rgb';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
FIGQUAL='';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;
tic;

COMPTYPEv=1;

for COMPTYPE=COMPTYPEv
	
	if COMPTYPE==1			% N/M=16		 % '8'*2, [4,8], [1,5] (if)
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='8';
		SEEDR=1:2;
		Mratio=4;
		RUNTYPEv=[1,5]; %1
	end
	
	
	% 1=w1, 2=w2, 3=w3 wavelet
	% 5=f1, 6=f1, 7=f3 fourier
	% 10=r1			   radial
	% 15=jh 16=jhvor   JH
	%RUNTYPEv = 7
	for RUNTYPE = RUNTYPEv
		EE=cell(6,1);
		if Mratio==4
			AVM='_4';
			M1=num2str(65536);
			M2=num2str(32761);
		elseif Mratio==8
			AVM='_8';
			M1=num2str(16384);
			M2=num2str(8281);
		end
		mstr='synth';
		if RUNTYPE==1
			RUNSTR='_w1';
			ylimit=[2e-5 0.2];
			C1=1e40;C2=4e7;p1=-3;p2=-5/3;
		elseif RUNTYPE==2
			RUNSTR='_w2';
			ylimit=[6e-8 0.01];
			C1=3*2e0;C2=5*1e7;p1=-5/3;p2=-3;
		elseif RUNTYPE==3
			RUNSTR='_w3';
			ylimit=[2e-5 0.2];
			C1=3e6;C2=4e7;p1=-3;p2=-5/3;
		elseif RUNTYPE==5
			RUNSTR='_f1';
			ylimit=[8e-6 8e-2];
			C1=1e40;C2=1e7;p1=-3;p2=-5/3;
		elseif RUNTYPE==6
			RUNSTR='_f2';
			ylimit=[2e-8 3e-3];
			C1=2e0;C2=1e7;p1=-5/3;p2=-3;
		elseif RUNTYPE==7
			RUNSTR='_f3';
			ylimit=[2e-6 6e-2];
			C1=1e6;C2=1e7;p1=-3;p2=-5/3;
		elseif RUNTYPE==10
			RUNSTR='_r1';  %r1, r1_s
			ylimit=[0.3 2e3];
			C1=1e40;C2=8e11;p1=-3;p2=-2;
		elseif RUNTYPE==15
			mstr='jh';
			RUNSTR='_'; %jh  Le1,Le2
			ylimit=[2e-2 1e4];
			C1=1;C2=2e12;p1=-3;p2=-5/3;
		elseif RUNTYPE==16  %fix
			mstr='jh';
			RUNSTR='_vor'; %jh
			ylimit=[4e-1 1e2];
			C1=1;C2=10e7;p1=-3;p2=-5/3+2;
		end
		RUNSTR4=RUNSTR;
		SIMN0t4=SIMN0t;
		if RUNTYPE==15
			RUNSTR4='_Le3'; %empty or Le1,2,3
			RUNSTR='_';
			%SIMN0t4='2';
		end
		
		fv=cell(1,1);
		ics=[1,2,4,5,6];
		for seedn=SEEDR
			seedstr=['_',num2str(seedn)];
			var1t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M1,'_',SIMN0t,'_',mstr,'_1_2_31',seedstr,RUNSTR,'.mat']);
			% 			var2t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M1,'_',SIMN0t,'_',mstr,'_1_2_312',seedstr,RUNSTR,'.mat']);
			% 			var5t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M1,'_',SIMN0t,'_',mstr,'_1_2_32',seedstr,RUNSTR,'.mat']);
			% 			var4t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M2,'_',SIMN0t4,'_',mstr,'_1_2_82',seedstr,RUNSTR4,'.mat']);
			% 			EEt{1}=var1t.E_cs_res1;
			% 			EEt{2}=var2t.E_cs_res1;
			% 			EEt{4}=var4t.E_cs_res1;
			% 			EEt{5}=var5t.E_cs_res1;
			% 			EEt{6}=var1t.E_os;
			% 			for i=ics
			% 				EE{i}=[EE{i};EEt{i}];
			% 			end
			if seedn==1
				fv=var1t.f_os;
			else
				fv=[fv;var1t.f_os];
			end
		end
		if RUNTYPE==15
			RUNSTR='_jh'; %jh
		end
		
		% ===========================
		
		Ns=size(var1t.f_os{1},1);
		J=log2(Ns);
		%kvec=var1t.k_os{1};
		[SIMN0,~]=size(fv);
		
		for ii=1:SIMN0
			f=fv{ii};
			%ind=1:2:Ns;
			
			
			% structure function
			rmax=81; % p/3 law valid for r in inertial range
			
			% save('./spectra_var/structure2d', 'S');
			for p=1:12
				S{p} = structure_2d_per(f,p,rmax); %direction (1,0)
				S{p} = S{p}/2 + structure_2d_per(f',p,rmax)/2; %direction (0,1)
			end
			
			prange=[1,2,4:12];
			for p=prange
				
				Ssc=abs(S{p})./abs(S{3});
				
				C0=0.6*(1/8).^(p/3-1)*Ssc(8);
				e0=(p/3-1);
				
% 				figure
% 				loglog(Ssc)
% 				hold all
% 				loglog(8:rmax, C0*(8:rmax).^e0,'k--')
% 				xlim([1 rmax])
				
				inguess=[e0,C0];
				ine=8:80;
				if RUNTYPE==5
					ine=5:50; %F
				end
				rvec=(1:length(Ssc))';
				fitscale=@(t) log(t);
				sco=exp_fit(Ssc(ine),rvec(ine),fitscale,inguess,'log');
				expon(p)=sco(1);
% 				loglog(ine,sco(2)*ine.^(-sco(1)),'r-')
			end
			expon=-expon;
			expon(3)=0;
			expon=expon+1;  % /S{3}
			prange=1:12;
			
			exc{RUNTYPE,ii}=expon;
			
		end
		
	end
	
	excm{1}=exc{1,ii}*0;
	excm{5}=excm{1};
	%excv{1}=excm{1};
	for ii=1:SIMN0
		excm{1}=excm{1}+exc{1,ii};
	%	excv{1}=[excv{1};exc{1,ii}];
		excm{5}=excm{5}+exc{5,ii};
	end
	excm{1}=excm{1}./SIMN0;
	excm{5}=excm{5}./SIMN0;
	%excv{1}=var(excv{1},0,1);
	
	figure;
	set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
	
	hhp=plot(prange,prange/3,':','Color',[0.4 0.4 0.4],'LineWidth',0.8);
	hold all
	hh=plot(prange,excm{1}(prange),'+-','Color',[0 0 0],'LineWidth',0.8);
	hhf=plot(prange,excm{5}(prange),'o-','Color',[0 0 0],'LineWidth',0.8);
	
	axis tight
	xlabel('$p$','FontSize', 11);
	ylabel('$\zeta_p$','FontSize', 11);
	
	%set(gca,'Outerposition',[-0.08,0,1.1,1])  %posx,posy,widx,widy
	%set(gca,'XTick',jvec)
	%set(gca,'YTick',[])
	ylim([1/3-0.1 4.3])
	xlim([0.85 12.15])
	
	legend([hh,hhf,hhp],{'Wavelet $(5/3)$','Fourier $(5/3)$','$p/3$'}...
		,'interpreter', 'latex','Location','southeast');
	
	if SAVEFIG==1
		export_fig([FIGLOC,'s2',FIGQUAL,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
		export_fig([FIGLOC,'s2',FIGQUAL,RUNSTR], '-eps', '-transparent', FIGCOLOR)
		disp([FIGLOC,'s2',FIGQUAL,RUNSTR])
	end
	
end
toc


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
