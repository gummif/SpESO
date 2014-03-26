% Calulations

if SAVEVAR==1
	% variables to save results
	f_cs_res1{SIMN0,SIMN}=[];
	f_os{SIMN0,SIMN}=[];
	E_cs_res1{SIMN0,SIMN}=[];
	E_os{SIMN0,SIMN}=[];
	k_os{1}=[];
	info{SIMN0,SIMN}=[];
end

%global CSPsi CSPsit CSA CSAt

% ======== wavelet transform operators, fwt/iwt : [N,1] -> [N,1]
if DWTTYPE==0
	% CDF 9/7 (not recommended)
	fwt=@(f) cdf97_1d(f,+1);
	iwt=@(f) cdf97_1d(f,-1);
	fwts=@(f) cdf97_1d(f,+1);
	iwts=@(f) cdf97_1d(f,-1);
elseif DWTTYPE==1
	% WAVELAB
	qmf=MakeONFilter('Coiflet',3);   % Coiflet 18
	%qmf=MakeONFilter('Symmlet',10);
	fwt=@(f) FWT_PO(f,0,qmf);
	iwt=@(f) IWT_PO(f,0,qmf);
	
	qmf2=MakeONFilter('Symmlet',6);
	%qmf2=MakeONFilter('Coiflet',2);
	%qmf2=MakeONFilter('Battle',1);
	fwts=@(f) FWT_PO(f,0,qmf2);
	iwts=@(f) IWT_PO(f,0,qmf2);
	
elseif DWTTYPE==2
	% WAVELET TOOLBOX
	f=zeros(N,1);
	dwtmode('per');
	Jmaxlev=log2(N)-4;
	[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('coif3');   % Coiflet 18
	[~,Lfwt] = wavedec(f,Jmaxlev,Lo_D,Hi_D);
	fwt=@(f) wavedec(f',Jmaxlev,Lo_D,Hi_D)';
	iwt=@(f) waverec(f',Lfwt,Lo_R,Hi_R)';

	[Lo_D2,Hi_D2,Lo_R2,Hi_R2] = wfilters('sym6');
	[~,Lfwts] = wavedec(f,Jmaxlev,Lo_D2,Hi_D2);
	fwts=@(f) wavedec(f',Jmaxlev,Lo_D2,Hi_D2)';
	iwts=@(f) waverec(f',Lfwts,Lo_R2,Hi_R2)';
else
	disp('unknown DWTTYPE')
	return
end

% identity operator
Iop=@(f) f;

% Spectrum operator
iN=2;
Kb=N/2048*2;
wintype='hanning';
Espec=@(fx) block_spectrum(fx,Kb,iN,wintype);


timing1(SIMN0,SIMN)=0;
disp(' ');

for k=1:SIMN
	M=Mvec(k);
	
	for k0=1:SIMN0
		
		% % initialize the RNG
		rng('default')
		temp=rand(1e4*k0,1);	% change seed with k0
		
		if JH==1
			f=load('./data/JH_2');
			f=f.jhu3;
			f=f-mean(f)+0.0;
			fhat=fft(f);
		elseif RHHDA==1
			f=load('./data/rhhda');
			f=f.rhhda;
			Ndif=length(f)-N;
			Nkl=floor(Ndif/SIMN0);  %change with k0
			f=f(1+(k0-1)*Nkl:N+(k0-1)*Nkl);
			f=f-mean(f)+0.0;
			fhat=fft(f);
		elseif SYNTH==1
			
			%rng('default')
			rand(2*N+2066,1);  % go some steps in seed forward
			
			get_first_dig=@(xx) floor(xx/10^floor(log10(xx)));
			if synthtype==1
				% FFT based
				syntype=2;
				if get_first_dig(syntslope)==1
					synextra=[1,-3];
				elseif get_first_dig(syntslope)==2
					synextra=[1,-5/3];
				end
				synp = synthetic_param(N,'1d',syntype,synextra);
				[f1,kvec1,E_org1] = synthetic_signal(N,'1d',synp);
				fhat=fft(f1);
				
				if get_first_dig(syntslope)==1
					synextra=[1,-5/3];
				elseif get_first_dig(syntslope)==2
					synextra=[1,-3];
				end
				synp = synthetic_param(N,'1d',syntype,synextra);
				[f2,kvec2,E_org2] = synthetic_signal(N,'1d',synp);
				fhat2=fft(f2);
				
				if syntslope>9 && mod(syntslope,10)==4
					k0s=N/4/2;
				elseif syntslope>9 && mod(syntslope,10)==2
					k0s=N/2/2;
				else
					k0s=N/16/2;
				end
				ik=find(kvec1>=k0s,1,'first');
				sratio=sqrt(mean(E_org1(ik-2:ik+2))/mean(E_org2(ik-2:ik+2)));

				fhat(k0s:N-k0s+2)=fhat2(k0s:N-k0s+2)*sratio;
				f=real(ifft(fhat));
				
				
			elseif synthtype==2
				% wavelet based
				syntype=2;
				
				%synex=[6/6,8/6,7/3];
				if get_first_dig(syntslope)==1
					synex=[8/6,10/6,9/3];  % -3 slope
				elseif get_first_dig(syntslope)==2
					synex=[4/6,6/6,5/3];  % 5/3 law
				end
				
				synp = synthetic_param_w(N,'1d',syntype,synex);
				f = synthetic_signal_w(N,'1d',synp);
				f=iwts(f);
				fhat=fft(f);
				
				iNs=2;
				Kbs=N/2048*8;
				%wintypes='welch';
				wintypes='hanning';
				[kvec,fspecf]=block_spectrum(f,Kbs,iNs,wintypes);
				fhat1=fhat;
				
				% second slope
				if get_first_dig(syntslope)==1
					synex=[4/6,6/6,5/3];  % 5/3 law
				elseif get_first_dig(syntslope)==2
					synex=[8/6,10/6,9/3];  % -3 slope
				end
				
				synp = synthetic_param_w(N,'1d',syntype,synex);
				f = synthetic_signal_w(N,'1d',synp);
				f=iwts(f);
				fhat2=fft(f);
				
				[kvec,fspecf2]=block_spectrum(f,Kbs,iNs,wintypes);
				if syntslope>9 && mod(syntslope,10)==4
					k0s=N/4/2;
				elseif syntslope>9 && mod(syntslope,10)==2
					k0s=N/2/2;
				else
					k0s=N/16/2;
				end
				ik=find(kvec>=k0s,1,'first');
				sratio=sqrt(fspecf(ik)/fspecf2(ik));

				fhat(k0s:N-k0s+2)=fhat2(k0s:N-k0s+2)*sratio;
				f=real(ifft(fhat));
				fhat=fft(f);
				
				% lower energy at high k
				k0s=N/4;
				f=synth_spec_scaledown(f,k0s,1.2); %40
				fhat=fft(f);
				% lower energy at small k
				k0L=100;
				fhat(2:k0L+1)=fhat(2:k0L+1)./linspace(3,1,k0L)'.^3;
				fhat(N:-1:N-k0L+1)=fhat(N:-1:N-k0L+1)./linspace(3,1,k0L)'.^3;

				k0L=40;
				fhat(2:k0L+1)=fhat(2:k0L+1)./linspace(3,1,k0L)'.^2;
				fhat(N:-1:N-k0L+1)=fhat(N:-1:N-k0L+1)./linspace(3,1,k0L)'.^2;
				fhat([2,3,N,N-1])=fhat([2,3,N,N-1])/5;
				f=real(ifft(fhat));
				
				[kvec,fspecf3]=block_spectrum(f,Kbs/8,iNs,wintypes);
			end
			
			% FFTs
			fhat=fft(f);
			
		end
		
		
		if MATRIX==1
			% Filter
			h = sign(2*rand(K,1)-1); %{-1,1}
			
			CSPsi=@(x) CS_rancon_db(x,M,N,h,0,fwt,iwt);  % any M
			CSPsit=@(x) CS_rancon_db(x,M,N,h,1,fwt,iwt);
			CSA=@(x) CS_rancon_db(x,M,N,h,0,Iop,Iop);
			CSAt=@(x) CS_rancon_db(x,M,N,h,1,Iop,Iop);
			
			CSA_fil=CSA;
			CSAt_fil=CSAt;
			% Filter
			%K_d=210;  %210,170
			%M_d=ceil((N+K_d-3)./do1vec(k));%-floor(K_d/do1vec(k));
			K_d=K;  % same
			M_d=M;  % same
			h = sign(2*rand(K_d,1)-1); %{-1,1}
			% dual measurement
			CSPsi_d=@(x) CS_rancon_db(x,M_d,N,h,0,fwt,iwt);  % any M
			CSPsit_d=@(x) CS_rancon_db(x,M_d,N,h,1,fwt,iwt);
			CSA_d=@(x) CS_rancon_db(x,M_d,N,h,0,Iop,Iop);
			CSAt_d=@(x) CS_rancon_db(x,M_d,N,h,1,Iop,Iop);
			
		elseif MATRIX==2
			% Convolution
			% 			h=exp(1i*rand(N,1)*2*pi); %random phase
			% 			h(1)=sign(imag(h(1)));
			% 			h(N/2+1)=sign(imag(h(N/2+1)));
			% 			h(N/2+2:end)=conj(h(N/2:-1:2));
			% 			ind=randperm(N);
			% 			ind=sort(ind(1:M));
			% 			CSPsi=@(x) CS_ranpha(x,M,N,h,0,fwt,iwt,ind);
			% 			CSPsit=@(x) CS_ranpha(x,M,N,h,1,fwt,iwt,ind);
			% 			CSA=@(x) CS_ranpha(x,M,N,h,0,Iop,Iop,ind);
			% 			CSAt=@(x) CS_ranpha(x,M,N,h,1,Iop,Iop,ind);
			disp('conv. not used')
			return
			
		end
		
		fapp=f*0;
		
		% ============= measurements ============
		sigma = 0.0005*0;
		e0 = sigma*randn(N,1);
		e = e0(1:M); % rand indp. of M
		% observations
		y = CSA(f) + e;
		
		% ========== dual measurements ==========
		sigma_d = 0.0005*0;
		e0_d = sigma_d*randn(N,1);
		e_d = e0_d(1:M_d);
		y_d=CSA_d(f) + e_d;
		
		% solve system
		if DECODER==1
			
			L=20;
			tol=1e-4;
			tolf=1e-6;
			tic;
			% LOMP
			fhap = OMP_gfa(y,CSPsi,CSPsit,N,1.3*M/2,L,tol,tolf);
			
			E_cs_res1{k0,k}=Espec(iwt(real(fhap)));
			
		elseif DECODER==31
			tic;
			
			% M/2 Best terms
			fhap = filter_coef(fwt(f),round(M/2));
			
			E_cs_res1{k0,k}=Espec(iwt(real(fhap)));
			
		elseif DECODER==312
			tic;
			
			% M Best terms
			fhap = filter_coef(fwt(f),round(M));
			
			E_cs_res1{k0,k}=Espec(iwt(real(fhap)));
			
		elseif DECODER==32
			tic;
			% Shannon sampling
			temp=f(1:round(N/M):end);
			fapp=interpft(temp,N);
			
			fhap=fwt(fapp);
			
			E_cs_res1{k0,k}=Espec(iwt(real(fhap)));
			
		elseif DECODER==82
			J=log2(N);
			J_M=round(floor(log2(M/2)));
			J0=5;
			
			beta=zeros(J-1,1) + 1;
			beta(J_M:end)=3;%3 %6,2,1
			alph=@(xx) 0.5*std(xx); %0.5
			
			
			tic;
			% Spectrum optimizer with QOMOMP
			
			% Initial L guess
			Ltype=1;
			Lj = Lvec_QOMOMP(N,M,'1d',Ltype,1);
			
			Lj0=Lj;
			if do1vec(k)==64  % N/64
				Lj0(J-1)=0;
				Lj0(J-8)=Lj0(J-8)*0.8;
				Lj0(J-7)=Lj0(J-7)*0.8;
				%Lj0(J-6)=Lj0(J-6)*0.7;
				Lj0(J-6:end)=Lj0(J-6:end)*0.7;
				Lj0(J-5:J-2)=[16,12,8,4];
				paramspec.jvec=J-6:J-2;
			elseif do1vec(k)==32  % N/32
				Lj0(J-1)=0;
				Lj0(J-6)=Lj0(J-6)*0.8;
				Lj0(J-6:end)=Lj0(J-6:end)*0.7;
				Lj0(J-2)=4;  %was 1
				paramspec.jvec=J-6:J-2;
			elseif do1vec(k)==16  % N/16
				Lj0(J-1)=0;
				Lj0(J-5)=Lj0(J-5)*0.5;
				Lj0(J-4)=Lj0(J-4)*1.2;
				Lj0(J-5:end)=Lj0(J-5:end)*0.7;
				paramspec.jvec=J-5:J-2;
			elseif do1vec(k)==8  % N/8
				Lj0(J-4)=Lj0(J-4)*0.6;
				Lj0(J-4:end)=Lj0(J-4:end)*0.7;
				%paramspec.jvec=J-4:J-2; Lj0(J-1)=0;
				paramspec.jvec=J-4:J-1;
				
			else
				disp(do1vec);
				return;
			end
			
			Lj0=round(Lj0);
			% QOMOMP operator
			tol=1e-2*2;
			tolf=1e-3/30; %30 %150 %any >100 ~ same
			Dspec=@(Ljx,tolfx) Dspec_function(@QOMOMPt_gfa,iwt,y,CSPsi,CSPsit,N,J0,Ljx,beta,alph,tol,tolfx);
			Dspec_d=@(Ljx,tolfx) Dspec_function(@QOMOMPt_gfa,iwt,y_d,CSPsi_d,CSPsit_d,N,J0,Ljx,beta,alph,tol,tolfx);
			
			tolf=tolf*10;
			specID=tic;
			% Optimize spectrum
			paramspec.tolf=tolf*1;  %*10
			paramspec.d=0.6;%0.4 %0.8
			paramspec.df=0.1;
			paramspec.iter=3;
			paramspec.res=5;
			paramspec.resf=1; %6 %final average
			[Lopt,Eopt] = spectrum_optim(Lj0,CSA_d,CSAt_d,Dspec,Espec,y_d,paramspec);
			toc(specID)
			
			fhap=Dspec(Lopt,tolf);  %opt solution
			Eopt=Espec(fhap);
			
			specID=tic;
			% dual
			[Lopt_d,Eopt_d] = spectrum_optim(Lj0,CSA,CSAt,Dspec_d,Espec,y,paramspec);
			fhap_d=Dspec_d(Lopt_d,tolf);
			Eopt_d=Espec(fhap_d);
			toc(specID)
			
			% ==== more clever than averaging ====
			[kvec,ff]=Espec(f);
			wj=find(kvec>2^7 & kvec<=2^(J-2));
			[~,E1]=Espec(CSAt_d(CSA_d(fhap)));
			[~,E2]=Espec(CSAt(CSA(fhap_d)));
			[~,Ef1]=Espec(CSAt_d(y_d));
			[~,Ef2]=Espec(CSAt(y));
			ne1=norm(E1(wj)-Ef1(wj))/norm(Ef1(wj));
			ne2=norm(E2(wj)-Ef2(wj))/norm(Ef2(wj));
			ne1L=norm(log(E1(wj)./Ef1(wj)))/norm(log(Ef1(wj)));
			ne2L=norm(log(E2(wj)./Ef2(wj)))/norm(log(Ef2(wj)));
			if ne1<ne2/2
				E_cs_res1{k0,k}=Eopt;
			elseif ne2<ne1/2
				E_cs_res1{k0,k}=Eopt_d;
			else
				E_cs_res1{k0,k}=Eopt/2+Eopt_d/2;
			end
			
			norms.ne1=ne1;norms.ne2=ne2;
			norms.ne1L=ne1L;norms.ne2L=ne2L;
			info{k0,k}=norms;
			
			fprintf('\n');
			
		end
		
		
		timing1(k0,k)=toc;
		[kvec,fspecf]=Espec(f);
		
		if SAVEVAR==1
			f_cs_res1{k0,k}=fhap;
			f_os{k0,k}=f;
			E_os{k0,k}=fspecf;
			k_os{1}=kvec;
		end
		
		fprintf('k0=%d ',k0);
	end
	fprintf('\n');
	fprintf('k=%d done\n',k);
end


if SAVEVAR==1
	% info in 82
	timing={timing1};
	save(SAVEFILE, 'f_cs_res1','f_os','E_cs_res1','E_os','k_os','timing','info');
end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
