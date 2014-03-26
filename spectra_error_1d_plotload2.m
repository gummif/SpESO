
%========== spectrum cases =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-rgb';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
FIGQUAL='';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;    




COMPTYPEv=[1,10,11];
%COMPTYPEv=[10,11];

for COMPTYPE=COMPTYPEv
	
	EXTRASTR='';
	if COMPTYPE==1					 % '64', [32,16,8,4], [1,2,3,4]
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='64';  %8,16,64
		Mratio=8; %16,8,4
		RUNTYPEv=[1,2,3,4];
		%RUNTYPEv=[4];
	elseif COMPTYPE==10				% '16', 8, 4 +s2/s4
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='16';  %8,16,64
		Mratio=8; %16,8,4
		RUNTYPEv=4;
		EXTRASTR='s2';
	elseif COMPTYPE==11
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='16';  %8,16,64
		Mratio=8; %16,8,4
		RUNTYPEv=4;
		EXTRASTR='s4'; %(av2)
	end
	
	
	for RUNTYPE = RUNTYPEv
		
		% 1=w353, 2=w533 wavelet
		% 3=f353, 4=f533 fourier
		
		% 	f_cs_res1{k0,k}=fhap;
		% 	f_os{k0,k}=f;
		% 	E_cs_res1{k0,k}=fhap;
		% 	E_os{k0,k}=fspecf;
		% 	k_os{1}=kvec;
		
		if Mratio==16
			FIGQ2='_16';
			M1=num2str(2066);
			M2=num2str(1033);
		elseif Mratio==8
			FIGQ2='_8';  % ''=m8
			M1=num2str(4132);
			M2=num2str(2066);
		elseif Mratio==4
			FIGQ2='_4';
			M1=num2str(8263);
			M2=num2str(4132);
		end
		if RUNTYPE==1
			RUNSTR=['_w353',EXTRASTR];
			%var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_31_w353.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			%var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_1_w353.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);
			n_good=41; % for QOM
			n_bad=57; %24
			
		elseif RUNTYPE==2
			RUNSTR=['_w533',EXTRASTR];
			% 	var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_31_w533.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			% 	var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_1_w533.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);
			n_good=27; % for QOM
			n_bad=34; %36
			
		elseif RUNTYPE==3
			RUNSTR=['_f353',EXTRASTR];
			%var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_31_f353.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			%var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_1_f353.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);
			n_good=8; % for QOM
			n_bad=14; %28
			
		elseif RUNTYPE==4
			RUNSTR=['_f533',EXTRASTR];
			% 	var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_31_f533.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312',RUNSTR,'.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32',RUNSTR,'.mat']);
			% 	var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_8_synth_284_1_1_f533.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82',RUNSTR,'_d2if.mat']);
			n_good=5; % for QOM
			n_bad=50; %81
			if strcmp(RUNSTR,'_f533s2')
				n_good=3; % for QOM
				n_bad=4;
			elseif strcmp(RUNSTR,'_f533s4')
				n_good=2; % for QOM
				n_bad=4;
			end
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
			if COMPTYPE==10		
				C1=60;C2=15*60e8;p1=-5/3;p2=-3;
			elseif COMPTYPE==11		
				C1=60;C2=6.5*60e8;p1=-5/3;p2=-3;
			end
		end
		
		
		%==============================
		
		N=length(var1b.f_os{1});
		J=log2(N);
		kvec=var1b.k_os{1};
		[SIMN0,~]=size(var1b.f_os);
		
		
		% find best/worst cases
% 			E_av2=var1b.E_cs_res1{1,1}*0;
% 			E_avo=E_av2*0;
% 			af=@(x)log(x);afi=@(x)exp(x); %log
% 			%af=@(x)x;afi=af;	% linear
% 			for i=1:SIMN0  % average spectra
% 				E_av2=E_av2+af(var4.E_cs_res1{i,1});
% 				E_avo=E_avo+af(var1b.E_os{i,1});
% 			end
% 		
% 			E_av2=afi(E_av2/SIMN0);
% 			E_avo=afi(E_avo/SIMN0);
% 		
% 		
% 			wj=find(kvec>2^9 & kvec<=2^(J-2));
% 			for i=1:SIMN0
% 				d{i}=log(var4.E_cs_res1{i,1}(wj)./E_av2(wj));  % E_av{4}
% 				Vd(i)=var(d{i});  % a smooth spectrum (E_av*C) has 0 error
% 		
% 				do{i}=log(var4.E_cs_res1{i,1}(wj)./E_avo(wj));
% 				Vdo(i)=norm(do{i});  % compare to org norm
% 			end
% 		
% 			[~,Vdi]=sort(Vd);
% 			[~,Vdoi]=sort(Vdo);
% 		
% 			for i=1:5
% 				figure,set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
% 				loglog(kvec,var1b.E_os{1,1})
% 				hold all
% 				loglog(kvec,var4.E_cs_res1{Vdi(i),1})
% 				loglog(kvec,var4.E_cs_res1{Vdoi(i),1})
% 				ylim(ylimit)
% 				xlim([50 N/2])
% 				title(['best ',num2str(i)])
% 			end
% 			for i=1:5
% 				figure,set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
% 				loglog(kvec,var1b.E_os{1,1})
% 				hold all
% 				loglog(kvec,var4.E_cs_res1{Vdi(end+1-i),1})
% 				loglog(kvec,var4.E_cs_res1{Vdoi(end+1-i),1})
% 				ylim(ylimit)
% 				xlim([50 N/2])
% 				title(['worst ',num2str(i)])
% 			end
% 			Vdoi(1)
% 			Vdoi(end)
% 		
% 		return
		
		figure;
		set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
		
		% original
		hh0=loglog(kvec,var1b.E_os{1,1},'-','Color',[0.7 0 0.7],'LineWidth',1.4);
		hold all
		
		% M-best
		hh1b=loglog(kvec,var1b.E_cs_res1{1,1},'-','Color',[0.2 0.2 1],'LineWidth',0.7);
		
		% Shannon
		hh2=loglog(kvec,var2.E_cs_res1{1,1},':','Color',[1 0.2 0.2],'LineWidth',1.3);
		
		% QOM bad
		hh4b=loglog(kvec,var4.E_cs_res1{n_bad,1},'--','Color',[0.6 0.6 0.6],'LineWidth',0.8);
		
		% QOM good
		hh4=loglog(kvec,var4.E_cs_res1{n_good,1},'-','Color',[0 0 0],'LineWidth',0.8);
		
		
		if COMPTYPE==1	
			loglog([100,1400],C1*3e-3*[100,1400].^p1,'k--','LineWidth',0.8)
			loglog([700,1e4],C2*3e-7*[700,1e4].^p2,'k--','LineWidth',0.8)
		elseif COMPTYPE==10		
			loglog([100,1e3*8*1.2],C1*3e-3*[100,1e3*8*1.2].^p1,'k--','LineWidth',0.8)
			loglog([1e3*8/1.2,1.5e4],C2*3e-7*[1e3*8/1.2,1.5e4].^p2,'k--','LineWidth',0.8)
		elseif COMPTYPE==11		
			loglog([100,1.04e3*4*1.2],C1*3e-3*[100,1.043e3*4*1.2].^p1,'k--','LineWidth',0.8)
			loglog([1.03e3*4/1.2,1.5e4],C2*3e-7*[1.03e3*4/1.2,1.5e4].^p2,'k--','LineWidth',0.8)
		end
		
		
		axis tight
		% \hspace{-15mm} $ \hspace{10mm}
		ylabel('\hspace{-2mm} $ E(k)$','FontSize', 11);
		xlabel('$k$','FontSize', 11);
		%set(gca,'Outerposition',[-0.08,0,1.1,1])  %posx,posy,widx,widy
		%set(gca,'XTick',jvec)
		
		ylim(ylimit)
		xlim([50 N/2])
		
		
		legend([hh0,hh1b,hh2,hh4,hh4b],{'original','$M$-best','Shannon',[CSNAME,' best'],[CSNAME,' worst']}...
			,'interpreter', 'latex','Location','southwest');
		
		
		if SAVEFIG==1
			%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,FIGQ2,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,FIGQ2,RUNSTR], '-eps', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'c1',FIGQUAL,FIGQ2,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'c1',FIGQUAL,FIGQ2,RUNSTR], '-eps', '-transparent', FIGCOLOR)
			
		end
		
		
	end
	
end


% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
