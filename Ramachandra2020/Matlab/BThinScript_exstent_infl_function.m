function BThinSript_exstent_infl_function(rate_degrade)
%function BThinSript_exstent_function
flag_exstent=1;
exs_thickness_fold=0.25;
exs_modulus_fold=10;
%exs_modulus_fold
exs_sizing=0.05;
flag_exs_degrade=1;

%_ves - vessel
%_exs - external stent
% rio_ves  = 0.647952178626704;						% inner radius at o (s = 0)
% roo_ves = 0.688107119933475;						% M-A radius at o (s = 0)

rio_ves  = 0.6468;						% inner radius at o (s = 0)
roo_ves = 0.6870;						% M-A radius at o (s = 0)

rio_exs=roo_ves+exs_sizing;%0.698107119933475;
roo_exs  = rio_exs+exs_thickness_fold*(roo_ves-rio_ves);	% outer radius at o (s = 0)
h=roo_exs-rio_exs

%rio_exs=0.698107119933475;
%roo_exs  = 0.708107119933475;						% outer radius at o (s = 0)
%
lzoo = 1;                           % axial stretch from o
hMo = roo_ves-rio_ves;						% medial thickness at o (s = 0)
hAo = roo_exs-rio_exs;						% adventitial thickness at o (s = 0)
%
%** PAR 
%

ce  = 89.7100;						% elastin modulus
Get = 1.9000;						% circumferential deposition stretch elastin
Gez = 1.6200;						% axial deposition stretch elastin
% Bt  = 0.0560;						% fraction of circumferential collagen
% within the adventitiaf
% Bz  = 0.0670;						% fraction of axial collagen within the adventitia
alp = 0.5220;						% orientation of diagonal collagen wrt axial direction
%
Ge = [1/Get/Gez Get Gez];			% elastin deposition stretches [r t z] in media and adventitia
%
%betaM = [Bz 1-Bz];					% medial betas [BzM 2*BdM]
%betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [BtA BzA 2*BdA]
%
%** parT0
%
c1m0 = 261.4000;					% c1t muscle at day 0
c2m0 = 0.24;					% c2t muscle
c1c0 = 234.9000;					% c1t collagen
c2c0 = 4.0800;					% c2t collagen
Gm0  = 1.20;					% circumferential deposition stretch (combined medial collagen and smc)
Gc0  = 1.25;					% deposition stretch (collagen)
%
%** parT4
%
% c1m4 = c1m0;					% c1t muscle at day 28 (week 4)
% c2m4 = c2m0;					% c2t muscle
% c1c4 = c1c0;					% c1t collagen
% c2c4 = c2c0;					% c2t collagen
% Gm4  = Gm0;					% circumferential deposition stretch (combined medial collagen and smc)
% Gc4  = Gc0;					% deposition stretch (collagen)
Tqoc = 80;							% collagen half-life (hypertensive)
etaq=1;
Tqom = (1/etaq)*Tqoc;				% smooth muscle half-life
kmco = [1/Tqom,1/Tqoc,1/Tqoc,1/Tqoc];%,...	% ko circ m M, ax c M, diag c M, diag c M
		%1/Tqoc,1/Tqoc,1/Tqoc,1/Tqoc];	% ko circ c A, ax c A, diag c A, diag c A
%
kmc = [kmco; kmco];					% rate parameters
%
KicM= 2;
%KfcM= 1.7404;
KfcM= 1.0; % ABR edited - set to 1, consistent with microstructure-independent infl.
           %Szafron et al., Tissue Eng. C manuscript

KscM = 2.5;
K = [KicM KscM KfcM];				% medial collagen gains [KicM KscM KfcM]
%
etaUm=0.8;
Kim = etaUm*KicM;					% smooth muscle gains [Kim Ksm Kfm]
Ksm = etaUm*KscM;
Kfm=etaUm*KfcM;
%Kfm = 1;                            %ABR  edits % see KfcM above for reference
%
etaUc=1.6667;
KicA =0.0;% etaUc*KicM;					% adventitial collagen gains [KicA KscA KfcA]
KscA =0.0;% etaUc*KscM;
KfcA =0.0;% etaUc*KfcM;
%
Ds = Tqom/40;						% G&R time increment
%
numHl = 20;							% number of Half Lives (Tqo) to consider for integration
domain = numHl*Tqom;				% total integration domain
invals = nearest(domain/Ds);		% Ds-length intervals within integration domain
breaks = invals+1;					% breaks within integration domain
%
t = (0:Ds:domain)';					% integration times for past histories
%
q = exp(-t*kmco);					% survival function
%
SimpsErr = simpsons(q,t(1),t(end),[])./(1./kmco) - ones(size(kmco)); %* ... should go to zero if properly integrated
%
%emco =[0.4714,0.4714,0.0572,0.0333,0.9667];
phioM = [0.25235 0.2627676 0.03430065 0.45058175];%[emco(1:2) emco(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phioA = 1;%[emco(4)   emco(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
rho = 1050;							% arterial mass density
%
rhooM = phioM*rho;					% medial [e mt cz 2*cd] densities
rhooA = phioA*rho;					% advent [e ct cz 2*cd] densities
%
rhoReM = rhooM(1);					% referential medial elastin density (constant)
rhoReA = rhooA(1);					% referential advent elastin density (constant)
rhoReA_orig=rhoReA;
%
%* original homeostatic referential density of fibers (smc and collagen)
%
rhoRf = ones(breaks,1)*[rhooM([2,3,4,4])];%,...	% circ m M, ax c M, diag c M, diag c M
						%rhooA([2,3,4,4])];		% circ c A, ax c A, diag c A, diag c A
%
mN = (ones(breaks,1)*kmco).*rhoRf;	% nominal mass production rates
%
Ups = ones(size(rhoRf));			% stimulus function
mR = mN.*Ups;						% referential mass production rate
%
ltTauInvM = ones(breaks,1);			% history of medial circum stretches
ltTauInvA = ones(breaks,1);			% history of advent circum stretches
lzTauInv  = ones(breaks,1);			% history of axial stretches
%
lTauInv = [ltTauInvM ltTauInvA lzTauInv];	% stretches history

%%
%[Pfo,so,Pfh,sh] = PreStr(PAR,parT0,parT4,riio,emco,riih,emch);	% P, f, stresses at o and h
c1m = c1m0*ones(breaks,1);			% past history of material properties
c2m = c2m0*ones(breaks,1);
c1c = c1c0*ones(breaks,1);
c2c = c2c0*ones(breaks,1);
Gm  =  Gm0*ones(breaks,1);
Gc  =  Gc0*ones(breaks,1);

PAR=[ce Get Gez alp];
parT0=[c1m0 c2m0 c1c0 c2c0 Gm0 Gc0];
riio=[rio_ves, roo_ves, rio_exs,roo_exs];
% if flag_exstent
% Po    = 12.9535;%13.871900523583244;						% inner pressure at o
% sgmIo = 532.4423;%545.551753553257;						% tr(sgm) at o
% fo    = 56.7882;%58.639973313884234;						% vessel axial force at o
% [Pfo,so] = PreStr_exstent_initialstate_only(PAR,parT0,riio);%,emco);	% P, f, stresses at o and h
% Po    = Pfo(1);						% inner pressure at o
% sgmIo = so(7);						% tr(sgm) at o
% fo    = Pfo(2);						% vessel axial force at o
%else
%[Pfo,so] = PreStr_exstent_initialstate_only(PAR,parT0,riio);%,emco);	% P, f, stresses at o and h
[Pfo,so] =  PreStr_exstent_initialstateOnly_unilayer(PAR,parT0,riio);%,emco);	% P, f, stresses at o and h
Po    = Pfo(1);						% inner pressure at o
sgmIo = so(7);						% tr(sgm) at o
fo    = Pfo(2);						% vessel axial force at o
stress_log(1,:)=[so(1),so(2),so(3),so(4),0.0,0.0,0.0,0.0,0.0];

%end
%
%Ph   = Pfh(1);						% inner pressure at h
fctr = 1.5;%Ph/Po;						% relative increase in P
%
Epso = 1;							% relative cardiac output at o
%
%[Pfo,so,Pfh,sh] = PreStr(PAR,parT0,parT4,riio,emco,riih,emch);	% P, f, stresses at o and h
days = 100*20;						% total simulation time (days)
%
steps = round(days/Ds+1);			% total G&R steps
%
ltM = 1; ltA = 1;					% circum stretches at o
JM  = 1; JA  = 1;					% Jacobians at o
%
%fctxi = 1/3;						% factor for inflammation evolution
%
%** initialize solution arrays
%
s    = 0:Ds:days;					% G&R times
ri   = rio_ves*ones(steps,1);			% inner radius
hM   = hMo*ones(steps,1);			% medial thickness
hA   = hAo*ones(steps,1);			% adventitial thickness
rhoR = ones(steps,1)*[rhooM([2,3,4,4])];%,...	% circ m M, ax c M, diag c M, diag c M
					  %rhooA([2,3,4,4])];	% circ c A, ax c A, diag c A, diag c A
%Upsi = ones(steps,8);				% stimulus function
Upsi = ones(steps,4);				% stimulus function
Eps  = Epso*ones(steps,1);			% relative cardiac output Q/Qo
lz   = ones(steps,1);				% axial stretch relative to state o
f    = fo*ones(steps,1);			% axial force
xi   = zeros(steps,1);				% prescribed inflammatory burden
%
%** prescribe inner pressure
%
wks4 = 28;							% 4 weeks = 28 days
SP = wks4/4;						% auxiliary time
%
P(s<=SP) = Po;%(1 + 1/2*(fctr-1)*(1-cos(pi*s(s<=SP)/SP)))*Po;
P(s> SP) = fctr*Po;%fctr*Po;

%
%** prescribe inflammatory burden
%
% SX = wks4/4; SX2 = wks4; SX3 = SX2+2*SX;	% auxiliary times
% fctrXi = 1/2;								% remaining inflammation
delta_duration1=0.21;                       %from Szafron et al.,ABME,2018
beta_skew=7;
t0_sigmoid=2;
rate_sigmoid=0.01;
m_max=2;
gamma_fn=delta_duration1.^beta_skew.*s.^(beta_skew-1).*exp(-delta_duration1*s);
gamma_fn=gamma_fn/max(gamma_fn);
sigmoid_fn=(m_max-1).*(1-exp(-rate_sigmoid.*(s-t0_sigmoid)./(m_max-1)));

xi=gamma_fn;
%xi=sigmoid_fn;
%xi=0.93*gamma_fn+0.25*sigmoid_fn; % scaling factors arbitrarily chose such that the max amplitude is ~ 1.
 
%
% xi(s<=SX) = 0.0;%1/2*(1-cos(pi*s(s<=SX)/SX));
% xi(s>SX & s<=SX2) = 0.0;%1;
% xi(s>SX2 & s<=SX3) = 0.0;%1 - fctrXi/2*(1-cos(pi*(s(s>SX2 & s<=SX3)-SX2)/(2*SX)));
% xi(s>SX3) = 0.0;%1 - fctrXi;
%
%***** COMPUTE GROWTH & REMODELING FOR s > 0
%
options = optimoptions(@fsolve,'Display','none','FunctionTolerance',1e-9,'StepTolerance',1e-9,'OptimalityTolerance',1e-9);
%
j = 0;								% maximum number of local iterations performed
%
flag_contact=0;
%for k = 2:steps %ABR edit
k=2; % switching from 'for' to 'while' loop because counter needs to be 
%reset for contact and for loop keeps internal counter - hard to reset!
%while is easy to reset
while k<=steps 
	%
	s(k) = s(k-1) + Ds;					% current G&R time
	
 if flag_exs_degrade
   t0=1;
   %rate_degrade=0.1;
   off_set=5;
   init = (1.0./(1.0 + exp(rate_degrade.*(t0 - off_set))));
   sigmoidal_degrade = (1.0./(1.0 + exp(rate_degrade.*(s(k) - off_set))))./init;
   rhoReA=rhoReA_orig*sigmoidal_degrade;
 end
  

%     if k<5
%         k
%     end
	%** update inflammation dependent properties for hereditary integrals
	%
	c1m = c1m0;%[c1m0 + (xi(k)^fctxi)*(c1m4-c1m0)
		   %c1m(1:end-1)];
	c2m = c2m0;%[c2m0 + (xi(k)^fctxi)*(c2m4-c2m0)
		   %c2m(1:end-1)];
	c1c = c1c0;%[c1c0 + (xi(k)^fctxi)*(c1c4-c1c0)
		   %c1c(1:end-1)];
	c2c = c2c0;%[c2c0 + (xi(k)^fctxi)*(c2c4-c2c0)
		   %c2c(1:end-1)];
	Gm  = Gm0;%[ Gm0 + (xi(k)^fctxi)*(Gm4-Gm0)
		   % Gm(1:end-1)];
	Gc  = Gc0;%[ Gc0 + (xi(k)^fctxi)*(Gc4-Gc0)
		    %Gc(1:end-1)];
	%
	%** shift past history arrays so that current state is located first
	%
	qaux = q;
	qaux(2:end,:) = qaux(1:end-1,:);
	kmc(2,:) = kmc(1,:);
	%
	rhoRf(2:end,:) = rhoRf(1:end-1,:);
	mN(2:end,:) = mN(1:end-1,:);
	mR(2:end,:) = mR(1:end-1,:);
	lTauInv(2:end,:) = lTauInv(1:end-1,:);
	Ups(2:end,:) = Ups(1:end-1,:);
	%
	rhoRfCor = zeros(1,size(rhoRf,2));	% trial referential densities
	%
	i = 0;								% local iteration
	%
	f(k) = f(k-1);						% trial axial force
	%
	%** perform local iterations until desired tolerance is attained
	%
	while norm(rhoRf(1,:)-rhoRfCor)/norm(rhoRf(1,:)) > eps(norm(rhoRf(1,:)))
		%
		%* compute circumferential stretches and axial force
		%
		%parL = [P(k),lz(k),sgmIo,JM,JA,rio,hMo,hAo,ce,Ge,rhoReM,rhoReA,alp,rho,t(1),t(end),Ds];
        %ABR edit -- changed to have modulus fold change,
        %flag_degradation,current time - s
        parL = [P(k),lz(k),sgmIo,JM,JA,rio_ves,hMo,hAo,ce,Ge,rhoReM,rhoReA,alp,rho,t(1),t(end),Ds,exs_modulus_fold,flag_exs_degrade,s(k)];
		%
		%[ltf,fval,eF] = fsolve(@(ltf) LaplaceCMM_externalStent(ltf,parL,mR,qaux,kmco,kmc,lTauInv,c1m,c2m,c1c,c2c,Gm,Gc),[ltM ltA f(k)],options);
        if flag_contact==0
        [ltf,fval,eF] = fsolve(@(ltf) LaplaceCMM_externalStent_unilayer(ltf,parL,mR,qaux,kmco,kmc,lTauInv,c1m,c2m,c1c,c2c,Gm,Gc),[ltM ltA f(k)],options);
        else
        [ltf,fval,eF] = fsolve(@(ltf) LaplaceCMM_externalStent(ltf,parL,mR,qaux,kmco,kmc,lTauInv,c1m,c2m,c1c,c2c,Gm,Gc),[ltM ltA f(k)],options);
        end
		%
		ltM  = ltf(1);									% medial circum stretch
		ltA  = ltf(2);									% advent circum stretch
		f(k) = ltf(3);									% vessel axial force
        if flag_contact==0
        stresslog=computeStress_externalStent_unilayer(ltf,parL,mR,qaux,kmco,kmc,lTauInv,c1m,c2m,c1c,c2c,Gm,Gc);
        else
        stresslog=computeStress_externalStent(ltf,parL,mR,qaux,kmco,kmc,lTauInv,c1m,c2m,c1c,c2c,Gm,Gc);
        end
		% stress-t-M, stress-t-A,stress-z-M, stress-z-A, theta-eqbm, z-eqbm, compatibility
        %
		%* additional variables from converged solution
		%
		lrM = JM/ltM/lz(k);								% medial radial strech
		lrA = JA/ltA/lz(k);								% advent radial strech
		%
		hM(k) = lrM*hMo;								% medial thickness
		hA(k) = lrA*hAo;								% advent thickness
		ri(k) = (ltM*(2*rio_ves+hMo)-hM(k))/2;				% inner radius
		if (ri(k)+hM(k))>= rio_exs && flag_contact==0 % if there is contact and contact is not detected then break out of the while loop
          break
        end
        % check for contact
		lTauInv(1,:) = [1/ltM	1/ltA	1/lz(k)];		% update current stretches
		%
		%* intramural and shear stress stimuli
		%
		%Dsgm = ( P(k)*ri(k)/(hM(k)+hA(k)) + f(k)/(pi*(hM(k)+hA(k))*(2*ri(k)+hM(k)+hA(k))) )/sgmIo-1;
        if flag_contact==0
        Dsgm = (stresslog(1) + stresslog(3))/sgmIo-1;
        else
        Dsgm = (stresslog(1) + stresslog(3))/(sgmIo)-1;
%        Dsgm = ( P(k)*ri(k)/(hM(k)+hA(k)) + f(k)/(pi*(hM(k)+hA(k))*(2*ri(k)+hM(k)+hA(k))) )/sgmIo-1;
        end
		
        %Dsgm = (stresslog(1) + stresslog(3))/sgmIo-1;
		Dtau = Eps(k)*rio_ves^3/ri(k)^3-1;
		%
		%* rate parameters, survival functions, and stimulus functions
		%
		kmc(1,:)   = kmco;%*(1+Dsgm^2);
		q(2:end,:) = qaux(2:end,:).*( ones(breaks-1,1) * exp(-Ds*mean(kmc)) );
		Ups(1,:)   = [ 1 +  Kim*Dsgm -  Ksm*Dtau +  Kfm*xi(k) ...
					  (1 + KicM*Dsgm - KscM*Dtau + KfcM*xi(k))*ones(1,3)];% ...
					  %(1 + KicA*Dsgm - KscA*Dtau + KfcA*xi(k))*ones(1,4)];    %* a row array
		%
		mR(1,:) = mN(1,:).*Ups(1,:);					% mass production before mass update
		%
		%* update/integrate mass densities and recompute mass production
		%
		rhoRfCor = simpsons(mR.*q,t(1),t(end),[]);		% corrected referential mass densities
		%
		mN(1,:) = kmc(1,:).*rhoRfCor;					% update nominal rates
		mR(1,:) = mN(1,:).*Ups(1,:);					% update mass production
		%
		rhoRf(1,:) = simpsons(mR.*q,t(1),t(end),[]);	% update referential mass densities
		%
		JM = ( rhoReM + sum(rhoRf(1,1:3)) ) / rho;		% update medial [e mt cz 2*cd] volume ratios
		if flag_exs_degrade
            JA=( rhoReA )/rho;  %for degrading polymer
        else
            JA = 1;%( rhoReA )/rho;%+ sum(rhoRf(1,5:7)) ) / rho;		% update advent [e ct cz 2*cd] volume ratios
        end
		%
		i = i+1;										% update iteration counter
		%
        
	end
	%
	%** converged values needed below for plots
	%
    rhoR(k,:) = rhoRf(1,:);				% referential mass densities
	Upsi(k,:) = Ups(1,:);				% stimulus functions
    stress_log(k,:)=[stresslog,Dsgm,Dtau];
    Jacob_stretch_log(k,:)=[JM, JA, lrM, lrA,ltM, ltA];
	%
    if (ri(k)+hM(k))>= rio_exs && flag_contact==0 
    k
    flag_contact=1;
    fprintf('\n Contact occured at step %d\n', k);
    k=k-2;% reset to previous step, so on increment you are at current step for recalc at current step
    end
    j = max(i,j);						% update max number of iterations
	%
    k=k+1;
       
end
%
%% a few plots
%
close all
scrsz = get(0,'ScreenSize');
figure
set(gcf,'position',[0.2*scrsz(3) 0.25*scrsz(4) 0.6*scrsz(3) 0.5*scrsz(4)])
%
%** inner pressure & cytokine dimensionless factor
%
subplot(4,4,1)
hold on
grid on
plot(s,P/P(1),'k','linew',1)
ylabel('$P/P_o$','interpreter','latex')
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days],'XTickLabel',{'','','','',''})
set(gca,'fontsize',13)
%
% subplot(4,4,5)
% hold on
% grid on
% plot(s,xi,'linew',1)
% set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days],'ylim',[-0.2 1.2])
% xlabel('$s$ [days]','interpreter','latex')
% ylabel('$\varrho_f/\varrho_{fm}$','interpreter','latex')
% set(gca,'fontsize',13)
%
%** stimulus function Upsilon for medial smooth muscle
%
subplot(242)
hold on
grid on
plot(s,Upsi(:,1),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\Upsilon^m_{ves}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** referential mass density of medial smc relative to homeostatic
%
subplot(243)
hold on
grid on
plot(s,rhoR(:,1)/rhoR(1,1),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\rho^m_{vesR}/\rho^m_{veso}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** medial wall thickness
%
subplot(244)
hold on
grid on
plot(s,hM/hMo,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$h_{veso}/h_{veso}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** luminal radius
%
subplot(245)
hold on
grid on
plot(s,ri/rio_ves,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$a/a_o$','interpreter','latex')
set(gca,'fontsize',13)
%
%** stimulus function Upsilon for adventitial collagen
% %
% subplot(246)
% hold on
% grid on
% plot(s,Upsi(:,5),'linew',1)
% set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
% xlabel('$s$ [days]','interpreter','latex')
% ylabel('$\Upsilon^c_A$','interpreter','latex')
% set(gca,'fontsize',13)
% %
% %** referential mass density of adventitial collagen relative to homeostatic
% %
% subplot(247)
% hold on
% grid on
% plot(s,rhoR(:,5)/rhoR(1,5),'linew',1)
% set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
% xlabel('$s$ [days]','interpreter','latex')
% ylabel('$\rho^c_{AR}/\rho^c_{Ao}$','interpreter','latex')
% set(gca,'fontsize',13)
%
%** deviation in stress
%
% %
subplot(246)
hold on
grid on
plot(s,stress_log(:,8),s,stress_log(:,9),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('\Delta stress')
legend('\sigma','\tau')
set(gca,'fontsize',13)

subplot(248)
hold on
grid on
plot(s,hA/hAo,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$h_{exs}/h_{exso}$','interpreter','latex')
set(gca,'fontsize',13)


figure
subplot(241)
plot(s,stress_log(:,1))
ylabel('ves - \theta')

subplot(242)
plot(s,stress_log(:,2))
ylabel('exs - \theta')

subplot(243)
plot(s,stress_log(:,3))
ylabel('ves-z')

subplot(244)
plot(s,stress_log(:,4))
ylabel('exs - z')

subplot(245)
plot(s,stress_log(:,5))
ylabel('\theta - eqbm')

subplot(246)
plot(s,stress_log(:,6))
ylabel('z - eqbm')

subplot(247)
plot(s,stress_log(:,7))
ylabel('compatibility')

subplot(248)
plot(s, Jacob_stretch_log(:,1),s,Jacob_stretch_log(:,2),'--',s,Jacob_stretch_log(:,3),s,Jacob_stretch_log(:,4),'--',s,Jacob_stretch_log(:,5),s,Jacob_stretch_log(:,6),'--')
legend('J-ves','J-exs','lr-ves','lr-exs','lt-ves','lt-exs')

%filename=['exstent_summary_',num2str(flag_exstent),'_gamma_',num2str(fctr),'_epso_',num2str(Epso),'_',num2str(exs_modulus_fold)];
%filename=['exstent_summary_',num2str(flag_exstent),'_gamma_',num2str(fctr),'_epso_',num2str(Epso),'_thickness_',num2str(exs_thickness_fold)];
%filename=['exstent_gamma_',num2str(fctr),'_epso_',num2str(Epso),'_thickness_percent_',num2str(exs_thickness_fold*100),'_modulus_percent_',num2str(exs_modulus_fold)];
%  filename=['exstent_gamma_',num2str(fctr),'_epso_',num2str(Epso),'_thickness_percent_',num2str(exs_thickness_fold*100),'_modulus_fold_',...
%      num2str(exs_modulus_fold)];
% filename=['exstent_gamma_',num2str(fctr),'_epso_',num2str(Epso),'_thickness_percent_',num2str(exs_thickness_fold*100),'_modulus_fold_',...
%     num2str(exs_modulus_fold),'_degrade_',num2str(rate_degrade*1e4)];


filename=['_exstent_gamma_',num2str(fctr),'_epso_',num2str(Epso),'_thickness_percent_',num2str(exs_thickness_fold*100),'_modulus_fold_',...
    num2str(exs_modulus_fold),'_degrade_',num2str(rate_degrade*1e4)];
save(filename)


end