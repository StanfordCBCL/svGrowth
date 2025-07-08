function P = LaplaceFull_externalStent(ro,par,mR,q,kmco,kmc,lTauInv,c1mt,c2mt,c1ct,c2ct,Gm,Gc,c1mc,c2mc,c1cc,c2cc)
%
%**	compute inner Pressure (P) for outer radius (ro) from Laplace equation
%   for full constrained mixture model with G&R 'frozen' at time s

%
%** RETRIEVE KNOWN VALUES
%
P      = par(1);							% inner pressure
lz     = par(2);							% axial stretch
sgmIo  = par(3);							% tr(sgm) at o
JM     = par(4);							% medial Jacobian
JA     = par(5);							% advent Jacobian
rio    = par(6);							% innner radius at o
hMo    = par(7);							% medial thickness at o
hAo    = par(8);							% advent thickness at o
ce     = par(9);							% elastin shear modulus
Ge     = par(10:12);						% deposition stretch of elastin
rhoReM = par(13);							% ref mass density of medial elastin
rhoReA = par(14);							% ref mass density of advent elastin
alp    = par(15);							% diagonal collagen orientation
rho    = par(16);							% artery mass density
t      = par(17:18);						% hereditary integral limits
Ds     = par(19);							% time step increment
%ABR edits
exs_modulus_fold=par(20);                                               % fold increase in modulus
ric    = par(23);							% current inner radius
hMc    = par(24);							% current med thickness
hAc    = par(25);							% current adv thickness
fz     = par(26);							% axial force on vessel
%
%** KINEMATICS / GEOMETRY
%
roo  = rio + hMo + hAo;						% outer radius at o
%
rMAc = ric + hMc;							% current M-A interface radius
roc  = ric + hMc + hAc;						% current outer radius
%
ri  = sqrt(ro^2+1/lz*(ric^2-roc^2));		% inner radii at P-d test
rMA = sqrt(ro^2+1/lz*(rMAc^2-roc^2));		% M-A radii at P-d test
%
hM = rMA-ri;								% medial thicknesses at P-d test
hA = ro-rMA;								% adventitial thicknesses at P-d test
%
lrMoc = hMc/hMo;							% incremental radial stretch media
lrAoc = hAc/hAo;							% incremental radial stretch adventitia
%
ltMoc = (2*ric+hMc)/(2*rio+hMo);			% incremental circumferential stretch media
ltAoc = (2*roc-hAc)/(2*roo-hAo);			% incremental circumferential stretch adventitia
%
ltM = (2*ri+hM)/(2*ric+hMc)*ltMoc;			% circumferential stretch media
ltA = (2*ro-hA)/(2*roc-hAc)*ltAoc;			% circumferential stretch adventitia
%
lrM = (1/ltM/lz)*lrMoc;						% radial stretch media
%lrA = (1/ltA/lz)*lrAoc;						% radial stretch adventitia
%
ltTauInvM = lTauInv(:,1);					% lTauInv(1,1) = 1/ltM
ltTauInvA = lTauInv(:,2);					% lTauInv(1,2) = 1/ltA
lzTauInv  = lTauInv(:,3);					% lTauInv(1,3) = 1/lz
%
FM = [lrM,ltM,lz];							% medial deformation grad
%FA = [lrA,ltA,lz];							% advent deformation grad
%
%** ELASTIN
%
FeM = FM.*Ge;								% medial elastin deformation grad
%FeA = FA.*Ge;								% advent elastin deformation grad
%ABR edit
%FeA = FA;
%
SeMod = ce*[1,1,1];							% second P-K stresses at constituent level
sgmWeM = rhoReM/(JM*rho)*FeM.*SeMod.*FeM;	% associated Cauchy stress [r,t,z] in media
%sgmWeA = rhoReA/(JA*rho)*FeA.*SeMod.*FeA*exs_modulus_fold;	% associated Cauchy stress [r,t,z] in advent
pM = sgmWeM(1);								% medial Lagrange multiplier
%pA = sgmWeA(1);								% advent Lagrange multiplier
%
%** SMC and COLLAGEN
%
%* circum, diag, and axial stretches (including G's)
%
lftM = ltM*ltTauInvM.*Gm;
lfdM = sqrt((ltM*ltTauInvM*sin(alp)).^2+(lz*lzTauInv*cos(alp)).^2).*Gc;
%ABR edits
%lftA = ltA*ltTauInvA.*Gc;
%lfdA = sqrt((ltA*ltTauInvA*sin(alp)).^2+(lz*lzTauInv*cos(alp)).^2).*Gc;
%lftA = ltA*ltTauInvA;
%lfdA = sqrt((ltA*ltTauInvA*sin(alp)).^2+(lz*lzTauInv*cos(alp)).^2);
lfz  = lz*lzTauInv.*Gc;
%
%* assign tension or compression constants
%
c1tM = c1mt; c2tM = c2mt;
c1tM(lftM<1) = c1mc(lftM<1);
c2tM(lftM<1) = c2mc(lftM<1);
c1dM = c1ct; c2dM = c2ct;
c1dM(lfdM<1) = c1cc(lfdM<1);
c2dM(lfdM<1) = c2cc(lfdM<1);
% c1tA = c1ct; c2tA = c2ct;
% c1tA(lftA<1) = c1cc(lftA<1);
% c2tA(lftA<1) = c2cc(lftA<1);
% c1dA = c1ct; c2dA = c2ct;
% c1dA(lfdA<1) = c1cc(lfdA<1);
% c2dA(lfdA<1) = c2cc(lfdA<1);
c1z = c1ct; c2z = c2ct;
c1z(lfz<1) = c1cc(lfz<1);
c2z(lfz<1) = c2cc(lfz<1);
%
%* auxiliary second-P-K-like stresses
%
SfModM(:,1) = c1tM.*(lftM.^2-1).*exp(c2tM.*(lftM.^2-1).^2).*Gm.^2;
SfModM(:,2) =  c1z.*(lfz.^2 -1).*exp( c2z.*(lfz.^2 -1).^2).*Gc.^2;
SfModM(:,3) = c1dM.*(lfdM.^2-1).*exp(c2dM.*(lfdM.^2-1).^2).*Gc.^2*sin(alp)^2;
SfModM(:,4) = c1dM.*(lfdM.^2-1).*exp(c2dM.*(lfdM.^2-1).^2).*Gc.^2*cos(alp)^2;
%SfModA(:,1) = c1tA.*(lftA.^2-1).*exp(c2tA.*(lftA.^2-1).^2).*Gc.^2;
%SfModA(:,2) =  c1z.*(lfz.^2 -1).*exp( c2z.*(lfz.^2 -1).^2).*Gc.^2;
%SfModA(:,3) = c1dA.*(lfdA.^2-1).*exp(c2dA.*(lfdA.^2-1).^2).*Gc.^2*sin(alp)^2;
%SfModA(:,4) = c1dA.*(lfdA.^2-1).*exp(c2dA.*(lfdA.^2-1).^2).*Gc.^2*cos(alp)^2;

%ABR edits
%SfModA(:,1) = 0.0;
%SfModA(:,2) = 0.0;
%SfModA(:,3) = 0.0;
%SfModA(:,4) = 0.0;

%
%* Cauchy stresses at constituent level
%
lfM = [ltM*ltTauInvM lz*lzTauInv ltM*ltTauInvM lz*lzTauInv];
%lfA = [ltA*ltTauInvA lz*lzTauInv ltA*ltTauInvA lz*lzTauInv];
%
sgmWfInt = [ 1/(JM*rho)*lfM.*SfModM.*lfM];%	1/(JA*rho)*lfA.*SfModA.*lfA ];
%
%* intramural stress stimulus, rate parameters, and survival functions
%
%Dsgm       = ( P*ric/(hMc+hAc) + fz/(pi*(hMc+hAc)*(2*ric+hMc+hAc)) )/sgmIo - 1;
breaks     = size(q,1);
kmc(1,:)   = kmco;%*(1+Dsgm^2);
q(2:end,:) = q(2:end,:).*( ones(breaks-1,1) * exp(-Ds*mean(kmc)) );
%
%* Cauchy stresses at mixture level
%
sgmWf  = simpsons(mR.*q.*sgmWfInt,t(1),t(end),[]);	% all of them
sgmWfM = sgmWf(1:4);								% media
%sgmWfA = [0 0 0 0];%sgmWf(5:8);								% adventitia
%
%** MIXTURE
%
sgmtM = sgmWeM(2) + sgmWfM(1) + sgmWfM(3) - pM;		% medial circum stress
%sgmtA = sgmWeA(2) + sgmWfA(1) + sgmWfA(3) - pA;		% advent circum stress
%
%** Pressure from corresponding equilibrium equation
%
P = (sgmtM*hM)/ri;						% luminal Pressure
%P = (sgmtM*hM + sgmtA*hA)/ri;						% luminal Pressure
%
end
