
%========== averaged spectra =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-rgb';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
FIGQUAL='';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;

COMPTYPEv=[1,2,3,4,10,11];
%COMPTYPEv=[4];

for COMPTYPE=COMPTYPEv
	
	EXTRASTR='';
	if COMPTYPE==1					 % '64', [32,16,8,4], [1,2,3,4]
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='64';  %8,16,64
		Mratio=4; %16,8,4
		RUNTYPEv=[1,2,3,4,10];
	elseif COMPTYPE==2
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='64';  %8,16,64
		Mratio=8; %16,8,4
		RUNTYPEv=[1,2,3,4,10];
	elseif COMPTYPE==3
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='64';  %8,16,64
		Mratio=16; %16,8,4
		RUNTYPEv=[1,2,3,4,10];
	elseif COMPTYPE==4
		SAVEDLOC='./spectra_var/if/';
		SIMN0t='64';  %8,16,64
		Mratio=32; %16,8,4
		RUNTYPEv=[1,2,3,4];  %not with 10
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
	
	%SAVEDLOC='./spectra_var/';
	% SAVEDLOC='./spectra_var/if/';
	% SIMN0t='16';  %8,16,64
	% Mratio=8; %16,8,4
	
	
	% 1=w353, 2=w533 wavelet
	% 3=f353, 4=f533 fourier
	% 10=rhh
	%RUNTYPEv=10;		%1,2,3,4,5
	for RUNTYPE = RUNTYPEv
		
		if Mratio==32
			AVM='_32';
			M1=num2str(1033);
			M2=num2str(517);
		elseif Mratio==16
			AVM='_16';
			M1=num2str(2066);
			M2=num2str(1033);
		elseif Mratio==8
			AVM='_8';
			M1=num2str(4132);
			M2=num2str(2066);
		elseif Mratio==4
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
		
		E_av(1:6)={var1.E_cs_res1{1,1}*0};
		E_min=E_av{4}+1e20;								% minimum
		af=@(x)log(x);afi=@(x)exp(x); %log
		%af=@(x)x;afi=af;	% linear
		% average spectra
		for i=1:SIMN0
			E_av{1}=E_av{1}+af(var1.E_cs_res1{i,1});
			E_av{2}=E_av{2}+af(var1b.E_cs_res1{i,1});
			E_av{3}=E_av{3}+af(var3.E_cs_res1{i,1});
			E_av{4}=E_av{4}+af(var4.E_cs_res1{i,1});
			E_min=min([E_min';var4.E_cs_res1{i,1}'])';	% minimum
			
			E_av{5}=E_av{5}+af(var2.E_cs_res1{i,1});
			E_av{6}=E_av{6}+af(var1.E_os{i,1});
		end
		
		for i=1:6
			E_av{i}=afi(E_av{i}/SIMN0);
		end
		
		
		%      logscale difference from average
		E_av2=E_av{4};
		nrem_max=round(0.3*SIMN0);
		nrem_max=0;  % org
		nremk=44; % k-means
		d={};Vd=[];
		for nrem=1:nrem_max   % 24
			wj=find(kvec>2^9 & kvec<=2^(J-2));
			for i=1:SIMN0
				d{i}=log(var4.E_cs_res1{i,1}(wj)./E_av2(wj));  % E_av{4}
				Vd(i)=var(d{i});  % a smooth spectrum (E_av*C) has 0 error
				%Vd(i)=mean(d{i}.^2);
				%d{i}=log(var4.E_cs_res1{i,1}(wj)./E_min(wj));     % min
				%Vd(i)=var(d{i});
				%Vd(i)=mean(d{i}.^2);
				%d{i}=log(var4.E_cs_res1{i,1}(wj)./E_av{6}(wj));
				%Vd(i)=norm(d{i});  % compare to org norm
			end
			%figure,plot(ones(SIMN0,1),Vd,'*')
			
			[~,Vind0]=sort(Vd);
			
			%====== least Var
			Vind=Vind0(1:end-nrem);  % remove nrem worst
			%====== kmeans Var
			% 		data=Vd(:);  %column
			% 		Classes=dcKMeans(data,2,[min(data);median(data)]);
			% 		xvec=1:length(Vd);
			% 		figure
			% 		hold all,
			% 		plot(xvec(Classes==1),data(Classes==1),'b+')
			% 		plot(xvec(Classes==2),data(Classes==2),'r+')
			% 		sum(Classes==1)
			% 		Vind=xvec(Classes==1);
			%====== k-means k wise
			% 		addpath('C:\Users\Lenovo\Documents\MATLAB\MAC\clusteringtoolbox','-end');
			% 		[score,Ek] = kmeans_spectra(var4.E_cs_res1,wj,@log);
			% 		[~,Vind]=sort(score);
			% 		Vind=Vind(1+nremk:end);
			% 		for c=1:SIMN0
			% 			Ek{c}=exp(Ek{c});
			% 		end
			%%%var4.E_cs_res1=Ek;	% bad freq. changed
			%======
			
			E_av2=E_av{4}*0;
			af=@(x)log(x);afi=@(x)exp(x); %log
			%af=@(x)x;afi=af;	% linear
			for i=Vind
				E_av2=E_av2+af(var4.E_cs_res1{i,1});
				%E_av2=E_av2+af(Ek{i,1});
			end
			E_av2=afi(E_av2/length(Vind));
			
		end
		% sorted Var
		%figure,	plot(Vd(Vind0),'*'),hold all,plot(Vd(Vind),'*')
		
		
		
		% plot all cases
% 		figure
% 		for i=setdiff(1:SIMN0,Vind)
% 			loglog(kvec,var4.E_cs_res1{i,1},'+')
% 			hold all
% 		end
% 		for i=Vind
% 			loglog(kvec,var4.E_cs_res1{i,1},'+k')
% 		end
% 		loglog(kvec,var4.E_os{1,1},'+r')
% 		loglog(kvec,E_av2,'g')
% 		
% 		
% 		% different averages types
% 		figure;
% 		set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
% 		
% 		% original
% 		hh0=loglog(kvec,E_av{6},'-','Color',[0.7 0 0.7],'LineWidth',1.4);
% 		hold all
% 		% QOM
% 		hh4=loglog(kvec,E_min,'--','Color',[0.6 0.6 0.6],'LineWidth',0.8);
% 		hh4=loglog(kvec,E_av{4},'-','Color',[0 0 0],'LineWidth',0.8);
% 		
% 		hh4=loglog(kvec,E_av2,'-','Color',[1 0 0.2],'LineWidth',0.8);
		
		
		% 	hh4=loglog(kvec,var4.E_cs_res1{15,1},'-','Color',[1 0 0],'LineWidth',0.8);
		% 	hh4=loglog(kvec,var4.E_cs_res1{14,1},'-','Color',[0 1 0],'LineWidth',0.8);
		
		%============================================================
		
		
		figure;
		if COMPTYPE==1 || COMPTYPE==2 || COMPTYPE==3 %three figures.
			if COMPTYPE==2 || COMPTYPE==3
				x_size=0.58*FIG_SIZE(1)-16;
			else
				x_size=0.58*FIG_SIZE(1);
			end
			if RUNTYPE==1||RUNTYPE==2||RUNTYPE==3
				y_size=0.58*FIG_SIZE(2)-16;
			else
				y_size=0.58*FIG_SIZE(2);
			end
			set(gcf, 'Position', [700 300 x_size y_size]);
		else
			set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
		end
		
		
		% original
		hh0=loglog(kvec,E_av{6},'-','Color',[0.7 0 0.7],'LineWidth',1.4);
		hold all
		
		% M-best
		hh1=loglog(kvec,E_av{2},'--','Color',[0.2 0.2 1],'LineWidth',1.3);
		
		% M/2
		hh1b=loglog(kvec,E_av{1},'--','Color',[0.2 0.2 1],'LineWidth',0.8);
		
		% shannon
		hh2=semilogy(kvec,E_av{5},':','Color',[1 0.2 0.2],'LineWidth',1.3);
		
		% LOMP
		hh3=loglog(kvec,E_av{3},':','Color',[0 0.6 0],'LineWidth',0.8);
		
		% QOM
		%hh4=loglog(kvec,E_av{4},'-','Color',[0 0 0],'LineWidth',0.8);
		hh4=loglog(kvec,E_av2,'-','Color',[0 0 0],'LineWidth',0.8);
		
		
		% loglog([100,1400],C1*3e-3*[100,1400].^p1,'k--','LineWidth',0.8)
		% loglog([700,1e4],C2*3e-7*[700,1e4].^p2,'k--','LineWidth',0.8)
		
		
		axis tight
		% \hspace{-15mm} $ \hspace{10mm}
		if ~(COMPTYPE==2 || COMPTYPE==3)
			ylabel('\hspace{-2mm} $ E(k)$','FontSize', 11);
		end
		if ~((COMPTYPE==1||COMPTYPE==2||COMPTYPE==3) && (RUNTYPE==1||RUNTYPE==2||RUNTYPE==3))
			xlabel('$k$','FontSize', 11);
			set(gca,'XTick',[10^2,10^3,10^4])
		end
		%set(gca,'Outerposition',[-0.08,0,1.1,1])  %posx,posy,widx,widy
		%set(gca,'XTick',jvec)
		
		ylim(ylimit)
		xlim([50 N/2])
		
		
		legend([hh0,hh1,hh1b,hh2,hh3,hh4],{'original','$M$-best','$M/2$-best','Shannon','LOMP',CSNAME}...
			,'interpreter', 'latex','Location','southwest');
		
		
		if SAVEFIG==1
			%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,'av',AVM,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,'av',AVM,RUNSTR], '-eps', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'a1',FIGQUAL,AVM,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'a1',FIGQUAL,AVM,RUNSTR], '-eps', '-transparent', FIGCOLOR)
			disp([FIGLOC,'a1',FIGQUAL,AVM,RUNSTR])
		end
		
		
	end
	
end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
