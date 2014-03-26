
%========== averaged spectra 2d =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-rgb';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
FIGQUAL='';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;

COMPTYPEv=[1,2,10];
%COMPTYPEv=1;

for COMPTYPE=COMPTYPEv

if COMPTYPE==1			% N/M=16		 % '8'*2, [4,8], [1,5] (if)
	SAVEDLOC='./spectra_var/if/';
	SIMN0t='8'; 
	SEEDR=1:2;
	Mratio=4; 
	RUNTYPEv=[1,2,3,5,6,7,15,16];
	%RUNTYPEv=[15];
elseif COMPTYPE==2		% N/M=64
	SAVEDLOC='./spectra_var/if/';
	SIMN0t='8'; 
	SEEDR=1:2;
	Mratio=8; 
	RUNTYPEv=[1,5];
elseif COMPTYPE==10		% N/M=16		% '1'*1, 4, 10
	SAVEDLOC='./spectra_var/';
	SIMN0t='1'; 
	SEEDR=1;
	Mratio=4; 
	RUNTYPEv=10;
elseif COMPTYPE==12			% N/M=16	*testing* less coef
	SAVEDLOC='./spectra_var/if/';
	SIMN0t='8'; %2 actually 
	SEEDR=1;
	Mratio=4; 
	RUNTYPEv=15;
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
	
	ics=[1,2,4,5,6];
	for seedn=SEEDR
		seedstr=['_',num2str(seedn)];
		var1t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M1,'_',SIMN0t,'_',mstr,'_1_2_31',seedstr,RUNSTR,'.mat']);
		var2t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M1,'_',SIMN0t,'_',mstr,'_1_2_312',seedstr,RUNSTR,'.mat']);
		var5t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M1,'_',SIMN0t,'_',mstr,'_1_2_32',seedstr,RUNSTR,'.mat']);
		var4t=load([SAVEDLOC,'spectra_error_2d_saved_1048576_',M2,'_',SIMN0t4,'_',mstr,'_1_2_82',seedstr,RUNSTR4,'.mat']);
		EEt{1}=var1t.E_cs_res1;
		EEt{2}=var2t.E_cs_res1;
		EEt{4}=var4t.E_cs_res1;
		EEt{5}=var5t.E_cs_res1;
		EEt{6}=var1t.E_os;
		for i=ics
			EE{i}=[EE{i};EEt{i}];
		end
		
	end
	if RUNTYPE==15
		RUNSTR='_jh'; %jh
	end
	
	% ===========================
	
	Ns=size(var1t.f_os{1},1);
	J=log2(Ns);
	kvec=var1t.k_os{1};
	[SIMN0,~]=size(EE{4});
	
	E_av=cell(6,1);
	E_av(1:6)={EE{1}{1,1}*0};
	E_min=E_av{4}+1e20;								% minimum
	af=@(x)log(x);afi=@(x)exp(x); %log
	%af=@(x)x;afi=af;	% linear
	% average spectra
	for i=1:SIMN0
		for j=ics
			E_av{j}=E_av{j}+af(EE{j}{i,1});
		end
		E_min=min([E_min';EE{4}{i,1}'])';	% minimum
	end
	
	for i=1:6
		E_av{i}=afi(E_av{i}/SIMN0);
	end
	
	
	%      logscale difference from average
	E_av2=E_av{4};
	nrem_max=round(0.3*SIMN0);
	nrem_max=0;  % org
	nremk=44; % k-means
	Vind=1:SIMN0;
	d={};Vd=[];
	for nrem=1:nrem_max   % 24
		wj=find(kvec>2^5 & kvec<=2^(J-2));
		for i=1:SIMN0
			d{i}=log(EE{4}{i,1}(wj)./E_av2(wj));  % E_av{4}
			Vd(i)=var(d{i});  % a smooth spectrum (E_av*C) has 0 error
			%Vd(i)=mean(d{i}.^2);
			%d{i}=log(EE{4}{i,1}(wj)./E_min(wj));     % min
			%Vd(i)=var(d{i});
			%Vd(i)=mean(d{i}.^2);
			%d{i}=log(EE{4}{i,1}(wj)./E_av{6}(wj));  % E_av{4}
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
			E_av2=E_av2+af(EE{4}{i,1});
			%E_av2=E_av2+af(Ek{i,1});
		end
		E_av2=afi(E_av2/length(Vind));
		
	end
	% sorted Var
	%figure,	plot(Vd(Vind0),'*'),hold all,plot(Vd(Vind),'*')
	
	
	
	% plot all cases
	figure
	for i=setdiff(1:SIMN0,Vind)
		loglog(kvec,EE{4}{i,1},'+')
		hold all
	end
	for i=Vind
		loglog(kvec,EE{4}{i,1},'+k')
		hold all
	end
	loglog(kvec,EE{6}{1,1},'+r')
	loglog(kvec,E_av2,'g')
	
% 	
% 	% different averages types
% 	figure;
% 	set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
% 	
% 	% original
% 	hh0=loglog(kvec,E_av{6},'-','Color',[0.7 0 0.7],'LineWidth',1.4);
% 	hold all
% 	% QOM
% 	hh4=loglog(kvec,E_min,'--','Color',[0.6 0.6 0.6],'LineWidth',0.8);
% 	hh4=loglog(kvec,E_av{4},'-','Color',[0 0 0],'LineWidth',0.8);
% 	
% 	hh4=loglog(kvec,E_av2,'-','Color',[1 0 0.2],'LineWidth',0.8);
	
	
	% 	hh4=loglog(kvec,var4.E_cs_res1{15,1},'-','Color',[1 0 0],'LineWidth',0.8);
	% 	hh4=loglog(kvec,var4.E_cs_res1{14,1},'-','Color',[0 1 0],'LineWidth',0.8);
	
	%============================================================
	
	
	figure;
	set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
	
	% original
	hh0=loglog(kvec,moving_average(E_av{6},1,5),'-','Color',[0.7 0 0.7],'LineWidth',1.4);
	hold all
	
	% M/2
	hh1b=loglog(kvec,moving_average(E_av{1},1,5),'--','Color',[0.2 0.2 1],'LineWidth',0.8);
	
	% M-best
	hh1=loglog(kvec,moving_average(E_av{2},1,5),'--','Color',[0.2 0.2 1],'LineWidth',1.3);
	
	% shannon
	hh2=semilogy(kvec,moving_average(E_av{5},1,5),':','Color',[1 0.2 0.2],'LineWidth',1.3);
	
	% LOMP
	% hh3=loglog(kvec,E_av{3},':','Color',[0 0.6 0],'LineWidth',0.8);
	
	% QOM
	%hh4=loglog(kvec,E_av{4},'-','Color',[0 0 0],'LineWidth',0.8);
	%hh4=loglog(kvec,E_av2,'-','Color',[0 0 0],'LineWidth',0.8);
	hh4=loglog(kvec,moving_average(E_av2,1,5),'-','Color',[0 0 0],'LineWidth',0.8);
	
	loglog([118,Ns/2-20],C1*3e-3*[118,Ns/2-20].^p1,'k--','LineWidth',0.8)
	if RUNTYPE==1 || RUNTYPE==5
		loglog([20,Ns/2-20],C2*3e-7*[20,Ns/2-20].^p2,'k--','LineWidth',0.8)
	elseif RUNTYPE==16
		loglog([2,40],C2*3e-7*[2,40].^p2,'k--','LineWidth',0.8)
	elseif RUNTYPE==10
		loglog([20,Ns/4],C2*3e-7*[20,Ns/4].^p2,'k--','LineWidth',0.8)
	else
		loglog([10,138],C2*3e-7*[10,138].^p2,'k--','LineWidth',0.8)
	end
	
	axis tight
	% \hspace{-15mm} $ \hspace{10mm}
	ylabel('\hspace{-2mm} $ E(k)$','FontSize', 11);
	xlabel('$k$','FontSize', 11);
	%set(gca,'Outerposition',[-0.08,0,1.1,1])  %posx,posy,widx,widy
	%set(gca,'XTick',jvec)
	
	ylim(ylimit)
	if RUNTYPE==16
		xlim([1 Ns/2-6])
	else
		xlim([7 Ns/2-6])
	end
	
	legend([hh0,hh1,hh1b,hh2,hh4],{'original','$M$-best','$M/2$-best','Shannon',CSNAME}...
		,'interpreter', 'latex','Location','southwest');
	
	
	if SAVEFIG==1
		%export_fig([FIGLOC,'spectra_error_2d',FIGQUAL,'av',AVM,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
		%export_fig([FIGLOC,'spectra_error_2d',FIGQUAL,'av',AVM,RUNSTR], '-eps', '-transparent', FIGCOLOR)
		export_fig([FIGLOC,'a2',FIGQUAL,AVM,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
		export_fig([FIGLOC,'a2',FIGQUAL,AVM,RUNSTR], '-eps', '-transparent', FIGCOLOR)
	end
	
	if RUNTYPE==10  % radial image
		figure;
		set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
		imagesc(1:Ns,1:Ns,var1t.f_os{1}(1:2:end,1:2:end))
		axis square;
		colormap(flipud(gray))
		set(gca,'XTickLabel',{' '})
		set(gca,'YTickLabel',{' '})
		if SAVEFIG==1
			%export_fig([FIGLOC,'spectra_error_2d',FIGQUAL,'im',AVM,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			%export_fig([FIGLOC,'spectra_error_2d',FIGQUAL,'im',AVM,RUNSTR], '-eps', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'a2',FIGQUAL,'i',AVM,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
			export_fig([FIGLOC,'a2',FIGQUAL,'i',AVM,RUNSTR], '-eps', '-transparent', FIGCOLOR)
		end
	end
	
end

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
