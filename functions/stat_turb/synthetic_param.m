function p = synthetic_param(N,dim,type,extra)
% Parameters synthetic signal synthetic_signal(N,dim,param)

if strcmp('1d',dim)==1
	if type==1
		synexpd=extra(1);
		p=struct('k0',50,'p0',1,'k1',round(5500*synexpd),'p1',-5/3,...
				'Ce',2.5e-5/synexpd,'pe',1.2,'C1',1e6,'n',0.1);	
	elseif type==2
		synexpd=extra(1);
		synslope=extra(2);
		p=struct('k0',50,'p0',1,'k1',round(N/3*synexpd),'p1',synslope,...
				'Ce',1/N/2/synexpd,'pe',1.2,'C1',1e6,'n',0.2);
	end
	
elseif strcmp('2d',dim)==1
	if type==1
		synexpd=extra(1);
		p=struct('k0',4,'p0',1,'k1',70*synexpd,'p1',-5/3,...
				'Ce',0.002/synexpd,'pe',1.2,'C1',1e6,'n',0.1,'nd',3);
	elseif type==2  % straight all way
		synexpd=extra(1);
		p=struct('k0',4,'p0',1,'k1',N/2-1,'p1',-5/3,...
				'Ce',0.002/synexpd,'pe',1.2,'C1',1e6,'n',0.1,'nd',5);
	elseif type==3  % straight with steep end
		synexpd=extra(1);
		p=struct('k0',4,'p0',1,'k1',N/2-100,'p1',-5/3,...
				'Ce',0.0003/synexpd,'pe',1.8,'C1',1e6,'n',0.1,'nd',5);
	elseif type==4  % straight with steep end, different slope
		synexpd=extra(1);
		synslope=extra(2);
		p=struct('k0',4,'p0',1,'k1',N/2-100,'p1',synslope,...
				'Ce',0.0003/synexpd,'pe',1.8,'C1',1e6,'n',0.1,'nd',5);
	end
	
end

end



% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
