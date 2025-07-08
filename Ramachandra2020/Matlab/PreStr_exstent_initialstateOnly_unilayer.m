function [Pfo,so] = PreStr_exstent_initialstateOnly_unilayer(PAR,parT0,riio)%,emco)
%
%*** computes mechanobiologically equilibrated stresses, pressure, and
%    transducer axial force at original (o) states
%

%
%** PAR
%
c   = PAR(1);						% c elastin
Get = PAR(2);						% circumferential deposition stretch elastin
Gez = PAR(3);						% axial deposition stretch elastin
alp = PAR(4);						% orientation of diagonal collagen wrt axial direction
%
% betaM = [Bz 1-Bz];					% medial betas [bzM 2*bdM]
% betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [btA bzA 2*bdA]
%
%-------------------------------  STATE o  --------------------------------
%
%** parT0
%
c1m = parT0(1);						% c1t muscle
c2m = parT0(2);						% c2t muscle
c1c = parT0(3);						% c1t collagen
c2c = parT0(4);						% c2t collagen
Gm  = parT0(5);						% circumferential deposition stretch (combined medial collagen and smc)
Gc  = parT0(6);						% deposition stretch (collagen)
%
%** emco and riio
%
phiMo = [0.25235 0.2627676 0.03430065 0.45058175];%[emco(1:2) emco(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
%phiAo = 1; % [emco(4)   emco(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
% rio  = riio(1);						% inner radius
% rMAo = riio(2);						% M-A radius
% roo  = riio(3);						% outer radius
%
rio_ves  = riio(1);						% inner radius at o (s = 0)
%rMAo = 0.676311018458466;						% M-A radius at o (s = 0)
roo_ves = riio(2);						% M-A radius at o (s = 0)
%roo  = 0.688107119933475;						% outer radius at o (s = 0)
% rio_exs=riio(3);
% roo_exs  = riio(4);	


lzoo = 1;							% axial stretch from o
%
% hMo = rMAo-rio;						% medial thickness
% hAo = roo-rMAo;						% adventitial thickness
hMo = roo_ves-rio_ves;						% medial thickness at o (s = 0)
%hAo = roo_exs-rio_exs;						% adventitial thickness at o (s = 0)
%

%
SMo = pi/lzoo*(roo_ves^2-rio_ves^2);		% medial cross-sectional area
%SAo = pi/lzoo*(roo^2-rMAo^2);		% adventitial cross-sectional area
%
%Ge = [1/Get/Gez 1/Get/Gez Get Get Gez Gez];	% [GerM GerA GetM GetA GezM GezA]
%ABR edit
Ge = [1/Get/Gez 1 Get 1 Gez 1];	% [GerM GerA GetM GetA GezM GezA]

%Ge = [1/Get/Gez 1 Get 1 Gez 1];	% [GerM GerA GetM GetA GezM GezA]
%
stMo = phiMo(1)*c*(Ge(3)^2-Ge(1)^2) + ...
	   phiMo(2)*c1m*(Gm^2-1)*exp(c2m*(Gm^2-1)^2)*Gm^2 + ...
	   phiMo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
%stAo = phiAo(1)*c*(Ge(4)^2-Ge(2)^2);
% + ...
% 	   phiAo(2)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
% 	   phiAo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
szMo = phiMo(1)*c*(Ge(5)^2-Ge(1)^2) + ...
	   phiMo(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiMo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
%szAo = phiAo(1)*c*(Ge(6)^2-Ge(2)^2);
% + ...
% 	   phiAo(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
% 	   phiAo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
sto = stMo;%(stMo*hMo+stAo*hAo)/(hMo+hAo);
szo = szMo;%(szMo*SMo+szAo*SAo)/(SMo+SAo);
%
Pfo(1,1) = (stMo*hMo)/rio_ves;%(stMo*hMo + stAo*hAo)/rio;			% P at o
Pfo(1,2) = szMo*SMo;%szMo*SMo + szAo*SAo;					% f at o
%
so = [stMo 0.0 szMo 0.0 sto szo sto+szo];		% stresses at o

end