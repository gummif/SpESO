
%========== error vs M =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-rgb';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
SAVEDLOC='./spectra_var/';
FIGQUAL='';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;


COMPTYPE=1;

if COMPTYPE==1					 % '64', [32,16,8,4], [1,2,3,4]
	SAVEDLOC='./spectra_var/if/';
	SIMN0t='64';  %8,16,64
	Mratio=8; %16,8,4
	RUNTYPEv=[1,2,3,4];
	%RUNTYPEv=[1];
end


for RUNTYPE = RUNTYPEv;
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
		RUNSTR='_w353';
		%do1vec=[8]; (-3,-5/3); qmf=MakeONFilter('Coiflet',3);
		var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_w353.mat']);
		var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_w353.mat']);
		var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_w353.mat']);
		var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_w353.mat']);
		var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_w353_d2if.mat']);
	elseif RUNTYPE==2
		RUNSTR='_w533';
		%do1vec=[8]; (-5/3,-3); qmf=MakeONFilter('Coiflet',3);
		var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_w533.mat']);
		var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_w533.mat']);
		var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_w533.mat']);
		var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_w533.mat']);
		var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_w533_d2if.mat']);
	elseif RUNTYPE==3
		RUNSTR='_f353';
		%do1vec=[8]; (-5/3,-3); qmf=MakeONFilter('Coiflet',3);
		var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_f353.mat']);
		var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_f353.mat']);
		var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_f353.mat']);
		var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_f353.mat']);
		var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_f353_d2if.mat']);
	elseif RUNTYPE==4
		RUNSTR='_f533';
		%do1vec=[8]; (-5/3,-3); qmf=MakeONFilter('Coiflet',3);
		var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_f533.mat']);
		var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_f533.mat']);
		var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_f533.mat']);
		var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_f533.mat']);
		var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_f533_d2if.mat']);
	end
	
	
	kvec=var1.k_os{1};
	%plot 1 spectrum
	ee=var1.E_os{1,1};
	ee1=var1.E_cs_res1{1,1};
	figure,loglog(kvec,ee),hold all,loglog(kvec,ee1)
	ee1b=var1b.E_cs_res1{1,1}; loglog(kvec,ee1b)
	ee2=var2.E_cs_res1{1,1}; loglog(kvec,ee2)
	ee3=var3.E_cs_res1{1,1}; loglog(kvec,ee3)
	ee4=var4.E_cs_res1{1,1}; loglog(kvec,ee4)
	
	
	%==============================
	
	N=length(var1.f_os{1});
	J=log2(N);
	kvec=var1.k_os{1};
	jvec=7:J-1;
	
	[SIMN0,~]=size(var1.f_os);
	err1=zeros(SIMN0,J-1);
	err1b=err1;err2=err1;
	err3=err1;err4=err1;
	for j=jvec
		wj=find(kvec>2^(j-1) & kvec<=2^j);
		for i=1:SIMN0
			err1(i,j) =norm(log(var1.E_os{i,1}(wj))-log(var1.E_cs_res1{i,1}(wj)))/length(wj);
			err1b(i,j)=norm(log(var1b.E_os{i,1}(wj))-log(var1b.E_cs_res1{i,1}(wj)))/length(wj);
			err2(i,j) =norm(log(var2.E_os{i,1}(wj))-log(var2.E_cs_res1{i,1}(wj)))/length(wj);
			err3(i,j) =norm(log(var3.E_os{i,1}(wj))-log(var3.E_cs_res1{i,1}(wj)))/length(wj);
			err4(i,j) =norm(log(var4.E_os{i,1}(wj))-log(var4.E_cs_res1{i,1}(wj)))/length(wj);
		end
	end
	err1=err1(:,jvec);
	err1b=err1b(:,jvec);
	err2=err2(:,jvec);
	err3=err3(:,jvec);
	err4=err4(:,jvec);
	
	%============================
	
	figure;
	set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
	
	
	% M/2 best
	hh1=semilogy(jvec,mean(err1,1),'--','Color',[0.2 0.2 1],'LineWidth',0.8);
	hold all
	errorbar(jvec,mean(err1,1),std(err1,1,1),std(err1,1,1),'x','Color',[0.2 0.2 1],'LineWidth',0.8)
	
	% M best
	hh1b=semilogy(jvec,mean(err1b,1),'--','Color',[0.2 0.2 1],'LineWidth',1.3);
	errorbar(jvec,mean(err1b,1),std(err1b,1,1),std(err1b,1,1),'x','Color',[0.2 0.2 1],'LineWidth',0.8)
	
	% shannon
	hh2=semilogy(jvec,mean(err2,1),':','Color',[1 0.2 0.2],'LineWidth',1.3);
	errorbar(jvec,mean(err2,1),std(err2,1,1),std(err2,1,1),'x','Color',[1 0.2 0.2],'LineWidth',0.8)
	
	% LOMP
	lstd=std(err3,1,1);
	temp=mean(err3,1);
	lstd(lstd>=mean(err3,1))=temp(lstd>=mean(err3,1))-1e-10;
	hh3=semilogy(jvec,mean(err3,1),':','Color',[0 0.6 0],'LineWidth',0.8);
	errorbar(jvec,mean(err3,1),lstd,std(err3,1,1),'x','Color',[0 0.6 0],'LineWidth',0.8)
	
	
	% QOM
	lstd=std(err4,1,1);
	temp=mean(err4,1);
	lstd(lstd>=mean(err4,1))=temp(lstd>=mean(err4,1))-1e-10;
	hh4=semilogy(jvec,mean(err4,1),'-','Color',[0 0 0],'LineWidth',1.3);
	errorbar(jvec,mean(err4,1),lstd,std(err4,1,1),'x','Color',[0 0 0],'LineWidth',0.8)
	
	% lines
	semilogy([11.5,11.5],[1e-5,1e3],'-.','Color',[0 0 0],'LineWidth',0.6);
	%semilogy([10.5,10.5],[1e-5,1e3],':','Color',[0 0 0],'LineWidth',0.8);
	
	grid on
	axis tight
	% \hspace{-15mm} $ \hspace{10mm}
	ylabel('\hspace{-2mm} $  \|  \log E-\log E^{\star} \|_{w_j} $','FontSize', 11);
	xlabel('$j$','FontSize', 11);
	set(gca,'Outerposition',[-0.08,0,1.1,1])
	set(gca,'XTick',jvec)
	
	ylim([9.5e-5 8e0])
	xlim([jvec(1)-0.2 jvec(end)+0.2])
	
	
	legend([hh1,hh1b,hh2,hh3,hh4],{'$M/2$-best','$M$-best','Shannon','LOMP',CSNAME}...
		,'interpreter', 'latex','Location','southeast');
	
	
	if SAVEFIG==1
		%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,FIGQ2,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
		%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,FIGQ2,RUNSTR], '-eps', '-transparent', FIGCOLOR)
		export_fig([FIGLOC,'e1',FIGQUAL,FIGQ2,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
		export_fig([FIGLOC,'e1',FIGQUAL,FIGQ2,RUNSTR], '-eps', '-transparent', FIGCOLOR)

	end
	
end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
