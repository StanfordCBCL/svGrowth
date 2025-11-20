K_c1=0.010000; % rate parameter - collage intramural
K_c2=8.006000; % rate parameter - collagen, shear
K_m1=2.009000; % rate parameter - smooth muscle, intramural
K_m2=6.007000; % rate parameter - smooth muscle, shear
C_ratio=20;
G_ch=1.039000;
G_mh=1.310000;
G_et=1.75;
G_ez=1.75;
kq_m=0.901000;
kq_c=0.010000;
T_M=109000.000000;


% Days to computer pressure-diameter curves
format short e

% Transmural pressure
P_h = 5.0*133.32237*10.0;  %(dyn/cm s^2)

% Vessel geometry
a_M = 0.2;           % (cm)
%L_M = 10;              % (cm)

% Viscosity and density of blood
mu = 0.037;             % (g/cm.s)
%rho_f = 1.050;         % (g/cm^3)

%K_f = 0.5;

% Homeostatic wall shear stress
%tau_wh = 50.6;%7.426;         % (g/cm.s^2) or (dynes/cm^2) or (0.1 Pa)
tau_wh = 6;
% Calculate Homeostatic flow rate and blood velocity
Q_Mh = tau_wh*pi*a_M^3/(4.0*mu);    %(mL/s)
%v_Mh = Q_Mh / (pi* a_M^2);              %(cm/s)


%-----------------------------------------------
% Solid material parameters
%-----------------------------------------------
nfiberfly=4;
% Mass density of solid constituents
rho_s = 1.050; % (g/cm^3)

% Fraction of collagen fibers (axial, cir, helical, helical)
c_frac = [0.1, 0.1, 0.4, 0.4];

% Alignment of collagen fibers
alpha_ckh = [0.0, 90.0, 45.0, 135.0]*(pi/180.0);

% Homeostatic circumferential stress
%sigma_h = 100.0*1000.0*10.0;      %N/cm^2
sigma_h = 18.0*1000.0*10.0;      %N/cm^2

% Calculate homeostatic thickness
h_h = P_h*a_M/sigma_h;
a_mid = a_M + h_h/2;

% Mass fractions
phi_c  = 0.42; %was 0.22
phi_ck = phi_c*c_frac;
phi_e  = 0.1; %was 0.02
phi_m  = 1.0 - phi_c - phi_e;

% Homeostatic masses
M_h   = h_h*rho_s;
M_ckh = M_h*phi_ck;
%M_ch = sum(M_ckh);
M_mh  = phi_m*M_h;
M_eh  = phi_e*M_h;

% sigma_ch = 100.0*1000.0*10.0;
% sigma_mh = 100.0*1000.0*10.0;  % pas + act

sigma_ch = 18.0*1000.0*10.0;
sigma_mh = 18.0*1000.0*10.0;  % pas + act
sigma_eh = (1.0/phi_e)*(sigma_h - phi_c*sigma_ch - phi_m*sigma_mh);


% Muscle activation parameters
Lambda_M=1.35;
Lambda_0=0.50;
C_basal=0.68;
T_S0=C_ratio*C_basal;
%      T_M = 58.0*1000.0*10.0;
% Passive response parameters and free parameter calculation

% Calculate remaining material parameters (g/cm.s^2)
c_m3=0.05025;
c_c2(1)=24022.65/phi_ck(1);
c_c2(2)=1000/phi_ck(2);
c_c2(3)=1.0;
c_c2(4)=1.0;
c_c3(1)=0.1;
c_c3(2)=0.05025;
c_c3(3)=1.035;
c_c3(4)=1.035;
     
c_e=sigma_eh/((G_et^2.0)-(((G_et^2.0)*(G_ez^2.0))^(-1)));
%      sigma_m_pas_h=c_m2*((G_mh**2.0)*(G_mh**(2.0)-1.0))*
%     2 EXP(c_m3*((G_mh**(2.0)-1.0)**2.0))  
%      sigma_act_h=sigma_mh-sigma_m_pas_h
%      temp=1.0-((Lambda_M-1.0)/(Lambda_M-Lambda_0))**2.0
%      T_M=sigma_act_h/((1.0-EXP(-(C_basal**2.0)))*temp)
sigma_act_h=T_M*(1.0 - exp(-C_basal^2.0))*(1.0 - ((Lambda_M - 1.0)/(Lambda_M - Lambda_0))^2.0);
disp(sigma_act_h);
if (sigma_act_h < 0.0)
       sigma_act_h = 0.0;  %only tension is possible
end

temp=c_m3*((G_mh^2.0-1.0)^2.0);
c_m2=(sigma_mh-sigma_act_h)/((G_mh^2.0)*(G_mh^2.0-1.0)*exp(temp));
disp(c_m2)
if(c_m2 <0.0)
      disp('c_m2 is negative')
      exit
end

for i=1:1:nfiberfly
      str_ch(i)=c_frac(i)*c_c2(i)*(G_ch^2.0)*((G_ch^2.0)-1.0)*exp(c_c3(i)*((G_ch^2.0)-1.0)^2.0)*(sin(alpha_ckh(i)))^2.0;
end
     
c_c2(4)=(sigma_ch-str_ch(1)-str_ch(2))/(2.0*str_ch(4));
c_c2(3)=c_c2(4);
c_c2_c=c_m2;
c_c3_c=c_m3;
c_m2_c=c_m2;
c_m3_c=c_m3;


%checking whether the equlibrium equation is satisfied with material
%parameters

dwdLt_c = 0.0;
for i = 1:4
    dwdLt_c = dwdLt_c + phi_ck(i)*h_h*dWkdLn(G_ch, c_c2(i), c_c3(i), c_c2_c, c_c3_c)*G_ch*(sin(alpha_ckh(i)))^2;
end

dwdLt_m = h_h*phi_m*dWmdLn(G_mh, c_m2, c_m3, c_m2_c, c_m3_c)*G_mh;
dwdLt_e = h_h*phi_e*dWedLn(G_et, G_ez, 1, c_e)*G_et;
dwdLt   = dwdLt_c + dwdLt_m + dwdLt_e;

T_act = sigma_act_h*phi_m*h_h;
F = dwdLt + T_act - P_h*a_M;  % this should be zero at equilibrium