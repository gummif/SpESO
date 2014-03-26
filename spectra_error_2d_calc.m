
% Calulations


if LOADVAR==1
	load(SAVEFILE);   % 'f_cs_res1', 'timing'
elseif SAVEVAR==1
	% variables to save results
	f_cs_res1{SIMN0,SIMN}=[];
	f_os{SIMN0,SIMN}=[];
	E_cs_res1{SIMN0,SIMN}=[];
	E_os{SIMN0,SIMN}=[];
	k_os{1}=[];
	info{SIMN0,SIMN}=[];
end

global CSPsi CSPsit CSA CSAt

% ======== wavelet transform operators, fwt/iwt : [Ns,Ns] -> [Ns,Ns]
% % WAVELAB
% % Coiflet 12
% qmf=MakeONFilter('Coiflet',2);
% fwt=@(f) FWT2_PO(f,0,qmf);
% iwt=@(f) IWT2_PO(f,0,qmf);
% % Symmlet 10
% %qmf2=MakeONFilter('Symmlet',5);
% %fwt2=@(f) FWT_PO(f,0,qmf2);
% %iwt2=@(f) IWT_PO(f,0,qmf2);
% qmfs=MakeONFilter('Symmlet',6);
% %qmf2=MakeONFilter('Battle',1);
% fwts=@(f) FWT2_PO(f,0,qmfs);
% iwts=@(f) IWT2_PO(f,0,qmfs);
% ========================
f=zeros(Ns,Ns);
dwtmode('per');
Jmaxlev=log2(Ns)-3;
% WAVELET TOOLBOX
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('coif2');
[~,Lfwt] = wavedec2(f,Jmaxlev,Lo_D,Hi_D);
fwt=@(f) waveletvec2mat(wavedec2(f,Jmaxlev,Lo_D,Hi_D),Lfwt);
iwt=@(f) waverec2(waveletmat2vec(f,Lfwt),Lfwt,Lo_R,Hi_R);

% WAVELET TOOLBOX
[Lo_D2,Hi_D2,Lo_R2,Hi_R2] = wfilters('sym6');
[~,Lfwts] = wavedec2(f,9,Lo_D2,Hi_D2);
fwts=@(f) waveletvec2mat(wavedec2(f,9,Lo_D2,Hi_D2),Lfwts);
iwts=@(f) waverec2(waveletmat2vec(f,Lfwts),Lfwts,Lo_R2,Hi_R2);
% ========================

% identity operator
Iop=@(f) f;
% the (:) operator
col=@(f) f(:);

% Spectrum operator
Espec=@(fx) spectrum_2d(fx);


timing1(SIMN0,SIMN)=0;
disp(' ');

for k=1:SIMN
	M=Mvec(k);
	Ms=Msvec(k);
	
	for k0=1:SIMN0
		
		% % initialize the RNG
		rng('default')
		temp=rand(1.44e5*k0,1);	% change seed with k0
		temp=rand(2.1e6*SEED,1);	% change seed with SEED
		
		
		if JH==1
			data=load('./data/JH_3d_1p1z4u1vor');
			f=data.f;
			vor=data.vor;
			clear data;
			if VORT==1
				f=vor;
			end
			
			% FFTs
			fhat=fft2(f);
			vorhat=fft2(vor);
			
		elseif SYNTH==1
			
			get_first_dig=@(xx) floor(xx/10^floor(log10(xx)));
			if synthtype==1
				% FFT based
				if get_first_dig(syntslope)==1  % one slope
					syntype=3;
					synextra=1;
					synp = synthetic_param(Ns,'2d',syntype,synextra);
					[f,kvec_org,E_org] = synthetic_signal(Ns,'2d',synp);
				elseif get_first_dig(syntslope)==2 || get_first_dig(syntslope)==3
					% (3,5/3) or (5/3,3)
					syntype=4;
					if get_first_dig(syntslope)==2
						synextra=[1,-3];
					elseif get_first_dig(syntslope)==3
						synextra=[1,-5/3];
					end
					synp = synthetic_param(Ns,'2d',syntype,synextra);
					[f1,kvec1,E_org1] = synthetic_signal(Ns,'2d',synp);
					fhat=fft2(f1);
					
					if get_first_dig(syntslope)==2
						synextra=[1,-5/3];
					elseif get_first_dig(syntslope)==3
						synextra=[1,-3];
					end
					synp = synthetic_param(Ns,'2d',syntype,synextra);
					[f2,kvec2,E_org2] = synthetic_signal(Ns,'2d',synp);
					fhat2=fft2(f2);
					
					if syntslope>9 && mod(syntslope,10)==8
						k0s=Ns/8/2;
					elseif syntslope>9 && mod(syntslope,10)==2
						k0s=Ns/2/2;
					else  % default
						k0s=Ns/4/2;
					end
					ik=find(kvec1>=k0s,1,'first');
					sratio=sqrt(mean(E_org1(ik-2:ik+2))/mean(E_org2(ik-2:ik+2)));
					kwave=round(fft2wavenumber(Ns,2));
					fhat(kwave>=k0s)=fhat2(kwave>=k0s)*sratio;
					f=real(ifft2(fhat));
					
				end
			elseif synthtype==2
				% wavelet based
				if get_first_dig(syntslope)==1  % one slope
					syntype=2;		% type of signal
					enslaw=5/3;
					synex=[4/6+1/2,6/6+1/2,enslaw];
					synp = synthetic_param_w(Ns,'2d',syntype,synex);
					f = synthetic_signal_w(Ns,'2d',synp);
					f=iwts(f);
					
					c=100;ks0=200;
					%f=synth_spec_scaledown(f,ks0,c);
				elseif get_first_dig(syntslope)==2 || get_first_dig(syntslope)==3
					% (3,5/3) or (5/3,3)
					syntype=2;
					if get_first_dig(syntslope)==2
						synex=[8/6+1/2,10/6+1/2,9/3];
					elseif get_first_dig(syntslope)==3
						synex=[4/6+1/2,6/6+1/2,5/3];
					end
					synp = synthetic_param_w(Ns,'2d',syntype,synex);
					f = synthetic_signal_w(Ns,'2d',synp);
					f=iwts(f);
					fhat=fft2(f);
					[kvec,fspecf]=spectrum_2d(f);
					
					% second slope
					if get_first_dig(syntslope)==2
						synex=[4/6+1/2,6/6+1/2,5/3];
					elseif get_first_dig(syntslope)==3
						synex=[8/6+1/2,10/6+1/2,9/3];
					end
					synp = synthetic_param_w(Ns,'2d',syntype,synex);
					f = synthetic_signal_w(Ns,'2d',synp);
					f=iwts(f);
					fhat2=fft2(f);
					[kvec,fspecf2]=spectrum_2d(f);
					
					if syntslope>9 && mod(syntslope,10)==8
						k0s=Ns/8/2;
					elseif syntslope>9 && mod(syntslope,10)==2
						k0s=Ns/2/2;
					else  % default
						k0s=Ns/4/2;
					end
					ik=find(kvec>=k0s,1,'first');
					sratio=sqrt(mean(fspecf(ik-2:ik+2))/mean(fspecf2(ik-2:ik+2)));
					kwave=round(fft2wavenumber(Ns,2));
					fhat(kwave>=k0s)=fhat2(kwave>=k0s)*sratio;
					f=real(ifft2(fhat));
					
				end
			elseif synthtype==3
				% radial
				Lr=13; %2.7  %circle
				f=synth_radial_2d(@(r) double(r<1),Ns,Lr*1.2,[300,600],1);
				f=f+synth_radial_2d(@(r) double(r<1),Ns,Lr/1.3,[800,200],1);
				f=f+synth_radial_2d(@(r) double(r<1),Ns,Lr,[600,800],1);
			end
			
			% FFTs
			fhat=fft2(f);
			
		end
		
		
		
		
		
		
		
		
		if MATRIX==1
			% Filter
			h = sign(2*rand(K,K)-1); %{-1,1}
			h=triu(h)+triu(h,1)';
			%h = round(7*rand(K,1)-3.5); %{-3,...,3} spikes in spectrum/more error
			% 				CSPsi=@(x) CS_rancon_2d(x,M,N,h,0,fwt,iwt);
			% 				CSPsit=@(x) CS_rancon_2d(x,M,N,h,1,fwt,iwt);
			% 				CSA=@(x) CS_rancon_2d(x,M,N,h,0,Iop,Iop);
			% 				CSAt=@(x) CS_rancon_2d(x,M,N,h,1,Iop,Iop);
			
			CSPsi=@(x) CS_rancon_db_2d(x,M,N,h,0,fwt,iwt);  % any M
			CSPsit=@(x) CS_rancon_db_2d(x,M,N,h,1,fwt,iwt);
			CSA=@(x) CS_rancon_db_2d(x,M,N,h,0,Iop,Iop);
			CSAt=@(x) CS_rancon_db_2d(x,M,N,h,1,Iop,Iop);
			
			% 				h = sign(2*rand(K,1)-1); %{-1,1}
			% 				CSPsi=@(x) CS_rancon_db_2ds(x,M,N,h,0,fwt,iwt);  % 1D stream
			% 				CSPsit=@(x) CS_rancon_db_2ds(x,M,N,h,1,fwt,iwt);
			% 				CSA=@(x) CS_rancon_db_2ds(x,M,N,h,0,Iop,Iop);
			% 				CSAt=@(x) CS_rancon_db_2ds(x,M,N,h,1,Iop,Iop);
			disp('dual missing')
			
		elseif MATRIX==2
			% Convolution
			h=exp(1i*rand(Ns,Ns)*2*pi); %random phase
			h = mat2conjsym(h);
			
			ind=randperm(N);
			ind=sort(ind(1:M));
			CSPsi=@(x) CS_ranpha_2d(x,M,N,h,0,fwt,iwt,ind);
			CSPsit=@(x) CS_ranpha_2d(x,M,N,h,1,fwt,iwt,ind);
			CSA=@(x) CS_ranpha_2d(x,M,N,h,0,Iop,Iop,ind);
			CSAt=@(x) CS_ranpha_2d(x,M,N,h,1,Iop,Iop,ind);
			
			% dual measurements
			h=exp(1i*rand(Ns,Ns)*2*pi); %random phase
			h = mat2conjsym(h);
			
			ind=randperm(N);
			ind=sort(ind(1:M));
			CSPsi_d=@(x) CS_ranpha_2d(x,M,N,h,0,fwt,iwt,ind);
			CSPsit_d=@(x) CS_ranpha_2d(x,M,N,h,1,fwt,iwt,ind);
			CSA_d=@(x) CS_ranpha_2d(x,M,N,h,0,Iop,Iop,ind);
			CSAt_d=@(x) CS_ranpha_2d(x,M,N,h,1,Iop,Iop,ind);
			
			% check for bias
			%figure,imagesc(abs(CSAt(CSA(f))))
			%figure,imagesc(abs(CSPsit(CSPsi(f))))
		end
		
		fapp=f*0;
		
		% ============= measurements ============
		sigma = 0.0005*0;
		e = sigma*randn(Ms,Ms);
		%e = sigma*randn(M,1);
		% observations
		y = CSA(f) + e;
		% initial guess = min energy
		%x0 = CSPsit(y);
		
		% ========== dual measurements ==========
		sigma_d = 0.0005*0;
		e_d = sigma*randn(Ms,Ms);
		%e = sigma*randn(M,1);
		% observations
		y_d = CSA_d(f) + e_d;
		
		
		J=log2(Ns);
		%Ljall=pow2(2*(1:J-1))'*3;   % ???
		
		% solve system
		if LOADVAR==0
			if DECODER==1
				L=20;
				tol=1e-3;
				tolf=1e-4;
				tic;
				% LOMP
				fhap = OMP_gfa(y,CSPsi,CSPsit,N,1.3*M/2,L);
				
			elseif DECODER==2
				J=log2(Ns);
				J_M=floor(log2(Ms/2));
				J0=J_M-4;
				Lj = zeros(J-1,1);
				Lj(1:J0-1)=pow2(2*(1:J0-1))*3;
				Lj(J0:J-1) = pow2(2*(J0:J-1))*3;
				Lj(J_M) = pow2(2*(J_M)-2)*3;
				Lj(J_M+1) = pow2(2*(J_M)-4)*3;
				Lj(J_M+2:end) = 0;
				if VORT==0
					Lj(J_M) = Lj(J_M)*0.8;
					Lj(J_M+1)=Lj(J_M+1)/8;
				else
					Lj(J_M) = Lj(J_M)*1.5; %1.5
					
					Lj(J_M+1) = Lj(J_M+1)*1.4;
					Lj(J_M+2:end) = 200;
				end
				
				%Lj(J_M+1:J-1)=round(pow2(J_M+(-5:-2.5:(-5-2.5*(J-J_M-2)))));
				
				%Lj=[ 2,4,8,16,28,63,125,235,133,12,0];  % opt4
				%Lj=[ 2,4,8,16,32,58,125,121,199,7,1];  % opt 5
				
				Lj=round(Lj);
				
				tic;
				% QOMOMP
				if ii==1 || (ii==2 && SHOW2==1)
					tol=1e-4; tolf=1e-7;
					tol=1e-2; tolf=1e-3;
					fhap = QOMOMP_2d_gfa(y,CSPsi,CSPsit,N,J0,Lj,tol,tolf);
					
				else % not convolution
					fhap=fwt(f_org);
					fhap(N/2:end)=0;
				end
			elseif DECODER==3
				tic;
				if ii==1
					% M/2 Best terms
					fhap = filter_coef(fwt(f_org),round(M/2));
				elseif ii==2
					% Shannon sampling
					temp=f_org(1:do1vec(k):end,1:do1vec(k):end);
					fapp=interpft(temp,Ns,1);
					fapp=interpft(fapp,Ns,2);
					
					fhap=fwt(fapp);
				end
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
				temp=f(1:round(Ns/Ms):end,1:round(Ns/Ms):end);
				fapp=interpft(temp,Ns,1);
				fapp=interpft(fapp,Ns,2);
				fhap=fwt(fapp);
				
				E_cs_res1{k0,k}=Espec(iwt(real(fhap)));
				
			elseif DECODER==8
				J=log2(Ns);
				J_M=floor(log2(Ms/2));
				J0=J_M-4;
				Lj = zeros(J-1,1);
				Lj(1:J0-1)=pow2(2*(1:J0-1))*3;
				Lj(J0:J-1) = pow2(2*(J0:J-1))*3;
				Lj(J_M) = pow2(2*(J_M)-2)*3;
				Lj(J_M+1) = pow2(2*(J_M)-4)*3;
				Lj(J_M+2:end) = 0;
				if VORT==0
					% AtA adj
					% 						Lj(J_M-1) = Lj(J_M-1)*1.5; %1.0
					% 						Lj(J_M) = Lj(J_M)*3.0; %0.8
					% 						Lj(J_M+1)=Lj(J_M+1)*0.5; %0.15
					
					Lj(J_M-1) = Lj(J_M-1)*1; %1.0
					Lj(J_M) = Lj(J_M)*0.8*4.8; %0.8
					Lj(J_M+1)=Lj(J_M+1)*0.15*200; %0.15
					
					Lj(J_M+2:end) = 5500*1;
					Lj(J-1) = 1000*0.8;
				else
					% 						Lj(J_M-1) = Lj(J_M-1)*0.85; %1
					% 						Lj(J_M) = Lj(J_M)*1.3; %1.5
					% 						Lj(J_M+1) = Lj(J_M+1)*1.7; %1.4
					% 						Lj(J_M+2:end) = 200*1.5; %1.0
					Lj(J_M-1) = Lj(J_M-1)*0.82*1.1; %1
					Lj(J_M) = Lj(J_M)*1.3*2.2; %1.8
					Lj(J_M+1) = Lj(J_M+1)*1.85*3.6; %3.2
					Lj(J_M+2:end) = 200*1.8*1.7; %1.0
					
				end
				
				%Lj(J_M+1:J-1)=round(pow2(J_M+(-5:-2.5:(-5-2.5*(J-J_M-2)))));
				
				%Lj=[ 2,4,8,16,28,63,125,235,133,12,0];  % opt4
				%Lj=[ 2,4,8,16,32,58,125,121,199,7,1];  % opt 5
				
				%Lj(J0:J-1) = pow2(2*(J0:J-1))*3;
				Lj=round(Lj);
				
				q=Lj*0-1;
				q(J-1)=0.014; % approx. like best M
				q(J-2)=0.14;
				q(J-3)=0.27;
				
				%vel
				%0.2590    0.0610         0
				%0.3253    0.2372    0.0384
				
				%vort
				%0.2775    0.1492    0.0171
				%0.3241    0.2577    0.1006
				
				beta=Lj*0+1;
				beta(J_M+1:end)=1; %6,2,1
				alph=@(xx) 0.5*std(xx); %0.5
				
				tic;
				% QOMOMP
				if ii==1 || (ii==2 && SHOW2==1)
					tol=1e-4; tolf=1e-7;
					tol=1e-2; tolf=1e-3;
					fhap = QOMOMPt_2d_gfa(y,CSPsi,CSPsit,N,J0,Lj,q,beta,alph,tol,tolf);
					
				else % not convolution
					fhap=fwt(f_org);
					fhap(N/2:end)=0;
					fhap(1:end)=0;
				end
			elseif DECODER==82
				J=log2(Ns);
				J_M=round(floor(log2(Ms/2)));
				J0=4;
				
				beta=zeros(J-1,1) + 1;
				beta(J_M:end)=2;%3 %6,2,1
				alph=@(xx) 0.5*std(xx(:)); %0.5
				
				
				tic;
				% Spectrum optimizer with QOMOMP
				
				% Initial L guess
				Ltype=1;
				Lj = Lvec_QOMOMP(Ns,Ms,'2d',Ltype,1);
				
				Lj0=Lj;
				if do1vec(k)==5.66  % N/5.66
					Lj0(end)=Lj0(end)*1.5;
					paramspec.jvec=J-4:J-1; %Lj0(J-1)=0;
					% more coef edit:
					Lj0(end-4)=Lj0(end-4)*0.93;
					Lj0(end-3)=Lj0(end-3)*1.0;
					Lj0(end-2)=Lj0(end-2)*1.8;
					Lj0(end-1)=Lj0(end-1)*2.2*1.5;
					Lj0(end)=Lj0(end)*4*2.3;
					
					if JH==1  %try less coefs
						if strcmp(EXTRA,'Le1')
							Lj0(end-2)=Lj0(end-2)*0.7;
							Lj0(end-1)=Lj0(end-1)*0.4;
							Lj0(end)=Lj0(end)*0.2;
						elseif strcmp(EXTRA,'Le2')
							Lj0(end-3)=Lj0(end-3)*0.7;
							Lj0(end-2)=Lj0(end-2)*0.4;
							Lj0(end-1)=Lj0(end-1)*0.2;
							Lj0(end)=Lj0(end)*0.1;		
						elseif strcmp(EXTRA,'Le3')
							Lj0(end-3)=Lj0(end-3)*0.7;
							Lj0(end-2)=Lj0(end-2)*0.3;
							Lj0(end-1)=Lj0(end-1)*0.1;
							Lj0(end)=Lj0(end)*0.0;	
							paramspec.jvec=J-4:J-2;
						end
					end
					
					paramspec.dimsc=2;  % allow more coefs. var.
				elseif do1vec(k)==11.31  % N/11.31
					Lj0(end)=0;
					paramspec.jvec=J-5:J-2; 
					% more coef edit:
					Lj0(end-5)=Lj0(end-5)*0.98;
					Lj0(end-4)=Lj0(end-4)*0.97;
					Lj0(end-3)=Lj0(end-3)*2;
					Lj0(end-2)=Lj0(end-2)*3.5;
					Lj0(end-1)=Lj0(end-1)*7.2;
					
					paramspec.dimsc=2;  % allow more coefs. var.
				else
					disp(do1vec);
					return;
				end
				
				Lj0=round(Lj0);
				q=Lj*0-1;
				q(J-1)=0.014; % approx. like best M
				q(J-2)=0.14;
				q(J-3)=0.27;
				% QOMOMP operator
				tol=1e-2*2;
				tolf=1e-3/30; %30 %150 %any >100 ~ same
				%Dspec=@(Ljx,tolfx) iwt(QOMOMPt_gfa(y,CSPsi,CSPsit,N,J0,Ljx,beta,alph,tol,tolfx));
				%Dspec_d=@(Ljx,tolfx) iwt(QOMOMPt_gfa(y_d,CSPsi_d,CSPsit_d,N,J0,Ljx,beta,alph,tol,tolfx));

				Dspec=@(Ljx,tolfx) Dspec_function(@QOMOMPt_2d_gfa,iwt,y,CSPsi,CSPsit,N,J0,Ljx,beta,alph,tol,tolfx,q);
				Dspec_d=@(Ljx,tolfx) Dspec_function(@QOMOMPt_2d_gfa,iwt,y_d,CSPsi_d,CSPsit_d,N,J0,Ljx,beta,alph,tol,tolfx,q);
				
				%QOMOMPt_2d_gfa(y,CSPsi,CSPsit,N,J0,Lj,q,beta,alph,tol,tolf);
				
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
				%[Lopt,Eopt] = spectrum_optim(Lj0,CSA,CSAt,Dspec,Espec,y,paramspec);
				toc(specID)
				
				%fhap=fwt(Dspec(Lj0,tolf));  %inital guess solution
				fhap=Dspec(Lopt,tolf);  %opt solution
				Eopt=Espec(fhap);  % no average!
				%E_cs_res1{k0,k}=Eopt;

				specID=tic;
				% dual
				[Lopt_d,Eopt_d] = spectrum_optim(Lj0,CSA,CSAt,Dspec_d,Espec,y,paramspec);
				fhap_d=Dspec_d(Lopt_d,tolf);
				Eopt_d=Espec(fhap_d);
				%E_cs_res1{k0,k}=Eopt/2+Eopt_d/2;
				toc(specID)
				
				% ==== more clever than averaging ====
				[kvec,ff]=Espec(f);
				wj=find(kvec>2^3 & kvec<=2^(J-2));
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
				%E_cs_res1{k0,k}=Eopt/2+Eopt_d/2;  % not av2 (worse)
				
% 				figure,loglog(kvec,Eopt,'r'),hold all,loglog(kvec,Eopt_d,'b'),loglog(kvec,ff,'k')
% 				lognormerror(E1,Ef1,kvec,wj)
% 				lognormerror(E2,Ef2,kvec,wj)
% 				figure,plot(kvec(wj),(E1(wj)-Ef1(wj))/norm(Ef1(wj)))
% 				figure,plot(kvec(wj),(E2(wj)-Ef2(wj))/norm(Ef2(wj)))
				% ====
				
				norms.ne1=ne1;norms.ne2=ne2;
				norms.ne1L=ne1L;norms.ne2L=ne2L;
				[~,E1]=Espec(CSAt_d(y_d));
				[~,E2]=Espec(CSAt_d(CSA_d(fhap)));
				[~,E1_d]=Espec(CSAt(y));
				[~,E2_d]=Espec(CSAt(CSA(fhap_d)));
				norms.ATAE1=E1;norms.ATAE2=E2;
				norms.ATAE1_d=E1_d;norms.ATAE2_d=E2_d;
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
			
		else
% 			fhap=f_cs_res1{k0,k};
% 			timing1(k0,k,ii)=timing{1}(k0,k,ii);
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
