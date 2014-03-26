
%========== slope errors table (slope range 1 (128-1000)) =========
set(0,'defaulttextinterpreter','latex');
FIGCOLOR='-cmyk';  % RGB,CMYK
FIG_SIZE(1)=560; %x (*1.8)<-(1,2)
FIG_SIZE(2)=420; %y
FIGLOC='./spectra_fig/';
FIGQUAL='_load3';
CSNAME='SpESO'; % Spectrum Estimation by Sparse Optimization
SAVEFIG=1;
FIGQ2='_s1';  % ''=m8



COMPTYPEv=1;

for COMPTYPE=COMPTYPEv
	
if COMPTYPE==1					 % '64', [32,16,8,4], [1,2,3,4]
	SAVEDLOC='./spectra_var/if/';
	SIMN0t='64';  %8,16,64
	Mrvec=[16,8,4];
	RUNTYPEv=[1,2,3,4];
	%RUNTYPEv=[1];
end



for RUNTYPE = RUNTYPEv; %1,2,3,4
	% 1=w353, 2=w533 wavelet
	% 3=f353, 4=f533 fourier
	
	%Mrvec=8;
	rr=1;
	for Mratio=Mrvec
		
		if Mratio==16
			M1=num2str(2066);
			M2=num2str(1033);
		elseif Mratio==8
			M1=num2str(4132);
			M2=num2str(2066);
		elseif Mratio==4
			M1=num2str(8263);
			M2=num2str(4132);
		end
		if RUNTYPE==1
			RUNSTR='_w353';
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_w353.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_w353.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_w353.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_w353.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_w353_d2if.mat']);
		elseif RUNTYPE==2
			RUNSTR='_w533';
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_w533.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_w533.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_w533.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_w533.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_w533_d2if.mat']);
		elseif RUNTYPE==3
			RUNSTR='_f353';
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_f353.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_f353.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_f353.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_f353.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_f353_d2if.mat']);
		elseif RUNTYPE==4
			RUNSTR='_f533';
			var1=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_31_f533.mat']);
			var1b=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_312_f533.mat']);
			var2=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_32_f533.mat']);
			var3=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M1,'_',SIMN0t,'_synth_284_1_1_f533.mat']);
			var4=load([SAVEDLOC,'spectra_error_1d_saved_32768_',M2,'_',SIMN0t,'_synth_284_1_82_f533_d2if.mat']);
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
		end
		
		
		N=length(var1.f_os{1});
		J=log2(N);
		kvec=var1.k_os{1};
		% jvec=J-5:J-3;
		% kmax=N/4;
		jvec=J-5;
		kmin=2^7;
		%inguess=[5/3,3e-6];
		
		
		[SIMN0,~]=size(var1.f_os);
		err=zeros(4,J-1);
		orgslope=err;
		% err1b=err1;err2=err1;
		% err3=err1;err4=err1;
		
		
		E_av(1:6)={var1.E_cs_res1{1,1}*0};
		%E_av{4}=E_av{4}+1e20;	% minimum
		af=@(x)log(x);afi=@(x)exp(x); %log
		%af=@(x)x;afi=af;	% linear
		% average spectra
		for i=1:SIMN0
			E_av{1}=E_av{1}+af(var1.E_cs_res1{i,1});
			E_av{2}=E_av{2}+af(var1b.E_cs_res1{i,1});
			E_av{3}=E_av{3}+af(var3.E_cs_res1{i,1});
			E_av{4}=E_av{4}+af(var4.E_cs_res1{i,1});
			%E_av{4}=min([E_av{4}';var4.E_cs_res1{i,1}'])';	% minimum
			E_av{5}=E_av{5}+af(var2.E_cs_res1{i,1});
			E_av{6}=E_av{6}+af(var1.E_os{i,1});	
		end
		%for i=[1,2,3,5]	% minimum
		for i=1:6
			E_av{i}=afi(E_av{i}/SIMN0);
		end
		
		%E_avM{rr}=E_av;
		
		% plot all
		figure
		loglog(kvec,var4.E_os{1,1},'r')
		hold all
		for i=1:SIMN0
			loglog(kvec,var4.E_cs_res1{i,1},'Color',[0.7 0.7 0.7])
		end
		loglog(kvec,var4.E_os{1,1},'r')
		loglog(kvec,E_av{4},'k')
		
		% for original
		wj=find(kvec>kmin & kvec<2^jvec(end));
		ys=[E_av{6}(wj(1)),E_av{6}(wj(end))];  % a line through range
		ks=kvec([wj(1),wj(end)]);
		inguesso(1)=-log(ys(1)/ys(2))/log(ks(1)/ks(2)); % exponent
		inguesso(2)=ys(1)/ks(1)^(-inguesso(1));			% constant
		
		for j=jvec
			%wj=find(kvec>2^j & kvec<=kmax);
			wj=find(kvec>kmin & kvec<2^j);
			for i=1:5
				
				ys=[E_av{i}(wj(1)),E_av{i}(wj(end))];
				ks=kvec([wj(1),wj(end)]);
				inguess(1)=-log(ys(1)/ys(2))/log(ks(1)/ks(2));
				inguess(2)=ys(1)/ks(1)^(-inguess(1));
				%loglog([100,1e4],inguess(2)*[100,1e4].^(-inguess(1)),'k--','LineWidth',0.8)
				
				%fitscale=@(t) sqrt(t);
				fitscale=@(t) log(t);
				s1=exp_fit(E_av{i}(wj),kvec(wj)',fitscale,inguess,'log');
				s1o=exp_fit(E_av{6}(wj),kvec(wj)',fitscale,inguesso,'log');
				err(i,j)=s1o(1)-s1(1);
				orgslope(i,j)=s1o(1);
				if i==4 && j==jvec(end)
					loglog([128,1e3],s1(2)*[128,1e3].^(-s1(1)),'b','LineWidth',2)
				end
				
			end
		end
		iorder=[4,5,1,2,3];
		err=err(iorder,jvec);
		orgslope=orgslope(iorder,jvec);
		
		% Mratio <-> rr
		errM{rr}=err;
		orgM{rr}=orgslope;
		rr=rr+1;
	end
	
	% special [16,8,4] at jvec(end)
	em16=errM{1}(:,end);
	em8=errM{2}(:,end);
	%em8(iorder==5)=errM{2}(iorder==5,1);  %shannon @ kvec=2048
	em4=errM{3}(:,end);
	%em4(iorder==5)=errM{3}(iorder==5,2);  %shannon @ kvec=4096
	err=[em16,em8,em4];
	orgslope=[orgM{1}(:,end),orgM{2}(:,end),orgM{3}(:,end)];
	
	%tablerr=flipud(err')./flipud(orgslope')*100;  % relative percentage
	tablerr=flipud(err');  %difference
	
	disp(tablerr)
	
	
	%============================
	
% 	figure;
% 	set(gcf, 'Position', [700 300 0.66*FIG_SIZE(1) 0.66*FIG_SIZE(2)]);
% 	
% 	
% 	% M/2 best
% 	hh1=plot(Mrvec,err(1,:),'--','Color',[0.2 0.2 1],'LineWidth',0.8);
% 	hold all
% 	% errorbar(jvec,mean(err1,1),std(err1,1,1),std(err1,1,1),'x','Color',[0.2 0.2 1],'LineWidth',0.8)
% 	
% 	% M best
% 	hh1b=plot(Mrvec,err(2,:),'--','Color',[0.2 0.2 1],'LineWidth',1.4);
% 	% errorbar(jvec,mean(err1b,1),std(err1b,1,1),std(err1b,1,1),'x','Color',[0.2 0.2 1],'LineWidth',0.8)
% 	%
% 	% % shannon
% 	% % hh2=semilogy(jvec,mean(err2,1),':','Color',[1 0.2 0.2],'LineWidth',1.3);
% 	% % errorbar(jvec,mean(err2,1),std(err2,1,1),std(err2,1,1),'x','Color',[1 0.2 0.2],'LineWidth',0.8)
% 	%
% 	% % LOMP
% 	% lstd=std(err3,1,1);
% 	% temp=mean(err3,1);
% 	% lstd(lstd>=mean(err3,1))=temp(lstd>=mean(err3,1))-1e-10;
% 	hh3=plot(Mrvec,err(3,:),':','Color',[0 0.6 0],'LineWidth',0.8);
% 	% errorbar(jvec,mean(err3,1),lstd,std(err3,1,1),'x','Color',[0 0.6 0],'LineWidth',0.8)
% 	%
% 	%
% 	% % QOM
% 	% lstd=std(err4,1,1);
% 	% temp=mean(err4,1);
% 	% lstd(lstd>=mean(err4,1))=temp(lstd>=mean(err4,1))-1e-10;
% 	hh4=plot(Mrvec,err(4,:),'-','Color',[0 0 0],'LineWidth',1.3);
% 	% errorbar(jvec,mean(err4,1),lstd,std(err4,1,1),'x','Color',[0 0 0],'LineWidth',0.8)
% 	
% 	% lines
% 	%semilogy([11.5,11.5],[1e-5,1e3],'-.','Color',[0 0 0],'LineWidth',0.6);
% 	%semilogy([10.5,10.5],[1e-5,1e3],':','Color',[0 0 0],'LineWidth',0.8);
% 	
% 	grid on
% 	axis tight
% 	% \hspace{-15mm} $ \hspace{10mm}
% 	ylabel('\hspace{-2mm} $    s-s^{\star}  $','FontSize', 11);
% 	xlabel('$N/M$','FontSize', 11);
% 	%set(gca,'Outerposition',[-0.08,0,1.1,1])
% 	set(gca,'XTick',fliplr(Mrvec))
% 	
% 	ylim([-4 4])
% 	xlim([Mrvec(end)-0.2 Mrvec(1)+0.2])
% 	
% 	legend([hh1,hh1b,hh3,hh4],{'$M/2$-best','$M$-best','LOMP',CSNAME}...
% 		,'interpreter', 'latex','Location','best');
	
	
	if SAVEFIG==1
		% not used
		%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,FIGQ2,RUNSTR], '-pdf', '-transparent', FIGCOLOR)
		%export_fig([FIGLOC,'spectra_error_1d',FIGQUAL,FIGQ2,RUNSTR], '-eps', '-transparent', FIGCOLOR)
		% table
		latextable(tablerr,...
			'format','%5.2f',...  % %4.0f
			'Horiz',{CSNAME,'Shannon','$M/2$-best','$M$-best','LOMP'},...
			'Vert',{'        ','4';RUNSTR,'8';'       ','16'},...
			'name',[FIGLOC,'table1d',FIGQUAL,FIGQ2,RUNSTR,'.tex'])
		
	end
	
	
end

end

% ================================================================








% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
