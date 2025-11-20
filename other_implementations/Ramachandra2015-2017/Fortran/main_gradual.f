c <^>
c-----
c(~ ~)
c /(` 
c   )
c Growth and Remodeling main file for the constrained mixture model
c Adapted from a MATLAB code from Prof.Jay Humphrey's lab
c For bugs, issues and errors email : 
c abbangal@ucsd.edu or amarsden@eng.ucsd.edu
c Nomenclature(mostly): 
c               _h - homeostatic value
c               _c - collagen related value
c               _m - smooth muscle related value
c               _e - elastin related value
c Variables in code are (almost) consistent with equations in the 
c publication: Valentin, A., L. Cardamone, S. Baek,   and 
c J. D. Humphrey."Complementary vasoactivity and matrix remodelling 
c in arterial adaptations to altered flow and pressure." Journal of 
c The Royal Society Interface 6, no. 32 (2009): 293-306.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c Notes on version of the code
c 1-Assumed constant Z stretch
c 2-Tension dependent degradation-SMC and collagen
c 3- EC Damage
c 4- SMC damage
      PROGRAM main_parameterstudy
      USE COM_GnR
      IMPLICIT NONE
!Interface all the functions used in the code
      INTERFACE

       FUNCTION q_i(s,kq)
       REAL*8 q_i
       REAL*8 s,kq
       END FUNCTION q_i
       
       FUNCTION dmg_smc(t,t_smcdamage,k_smcdamage,smc_power_coeff)
       REAL*8 dmg_smc
       REAL*8 t,t_smcdamage,k_smcdamage,smc_power_coeff
       END FUNCTION dmg_smc

       FUNCTION f_beta(beta,kq)
       REAL*8 f_beta
       REAL*8 beta,kq
       END FUNCTION f_beta
       
       FUNCTION rm_m(n_SMC,m_basal,dn_stress,dn_tau,Zeta_m,Km2)
       REAL*8 rm_m
       REAL*8 n_SMC,m_basal,dn_stress,dn_tau,Zeta_m,Km2
       END FUNCTION rm_m

       FUNCTION rm_c(n_c,m_basal,dn_stress,dn_tau,Zeta_c,Kc2)
       REAL*8 rm_c
       REAL*8 n_c,m_basal,dn_stress,dn_tau,Zeta_c,Kc2
       END FUNCTION rm_c

       FUNCTION dWedLn(Ln_t,Ln_z,ax,c_e)
       REAL*8 dWedLn
       REAL*8 Ln_t,Ln_z,ax,c_e
       END FUNCTION dWedLn
    
       FUNCTION dWkdLn(Ln,c_c2,c_c3,c_c2_c,c_c3_c)
       REAL*8 dWkdLn
       REAL*8 Ln,c_c2,c_c3,c_c2_c,c_c3_c
       END FUNCTION dWkdLn
     
       FUNCTION dWmdLn(Ln,c_m2,c_m3,c_m2_c,c_m3_c)
       REAL*8 dWmdLn
       REAL*8 Ln,c_m2,c_m3,c_m2_c,c_m3_c
       END FUNCTION dWmdLn

       FUNCTION ddWeddLn(Ln_t,Ln_z,f,c_e)
       REAL*8 ddWeddLn
       REAL*8 Ln_t,Ln_z,f,c_e
       END FUNCTION ddWeddLn
    
       FUNCTION ddWkddLn(Ln,c_c2,c_c3,c_c2_c,c_c3_c)
       REAL*8 ddWkddLn
       REAL*8 Ln,c_c2,c_c3,c_c2_c,c_c3_c
       END FUNCTION

       FUNCTION ddWmddLn(Ln,c_m2,c_m3,c_m2_c,c_m3_c)
       REAL*8 ddWmddLn
       REAL*8 Ln,c_m2,c_m3,c_m2_c,c_m3_c
       END FUNCTION

      END INTERFACE

c modeling & simulation parameters  
      REAL*8 comp_t,delta_P,delta_Q,K_c1,K_c2,K_m1,K_m2,C_ratio
      REAL*8 G_ch,G_mh,G_et,G_ez,kq_m,kq_c


c Variables used for computation in the code
      REAL*8 P_h,tau_wh,sigma_h,h_h,Q_Mh,rho_s,a_M,mu
      REAL*8 phi_c,phi_e,phi_m,M_h,M_mh,M_eh,t,h
      REAL*8 sigma_mh,sigma_eh,sigma_ch,sigma_m_pas_h,sigma_act_h
      REAL*8 c_m2,c_m3,c_e,beta_m,beta_c,beta,m_basal_m
      REAL*8 Lambda_M,Lambda_0,C_basal,T_S0,T_M,T_act,T_S,y_Lkn,y_Lmn
      REAL*8 c_c2_c,c_c3_c,c_m2_c,c_m3_c,k_act,a_t,a_act_p,da_act_p
      REAL*8 P,Q_M,tol_m,tol_a,tau_w,F,n_tau,dn_tau,wt,a_act
      REAL*8 Lt_tau,Lz_tau,L_t,L_z,dt,L_m_act,J,C_t,dC_tda,dT_actda
      REAL*8 total_M, F_a,dF_ada,nm_stress,dn_C,T_t_c,T_z_c,T_t_m
      REAL*8 mean_age_m,mean_age_c,t_gradual
      REAL*8, ALLOCATABLE :: c_frac(:),alpha_ckh(:),phi_ck(:),
     2 c_c2(:),c_c3(:),str_ch(:),M_ckh(:),sigma_ck(:),alpha_ck(:),
     2 M_ck(:),nc_stress(:)

      INTEGER i,n_t,num_it1,num_it2,tn0,ept
      REAL*8 temp,temp1,temp2,temp3,temp4,temp5(4)
      CHARACTER*50 filename1,filename2,filename3
     

c Variables for calculating stress in collagen
      REAL*8 sQ_ck,dLndLt,dLndLz,ddLnddLt,mc_tau,M_c
      REAL*8, ALLOCATABLE :: dWck_dLn(:),alpha_tau(:),Lc_k_tau(:),
     2 Lc_k(:),Lc_kn(:),m_basal_ck(:),mp_c(:),dWcdLn_h(:)

c Variables for calculating stress in SMC
      REAL*8 Lm_n,dWm_dLn,sQ_m,M_m,mm_tau,dWmdLn_h
c Variables for calculating stress in elastin
      REAL*8 M_e,ax,n_c,n_smc,f_flag,mp_c_file(4),m_m_file,tau_w_file,
     2 tau_w_file2 
c Variables for flags,tension dependent degradation, ec damage, smc damage
      REAL*8 ept_flag,flagfile_exist,flag_tension_degrade,flag_ecdamage,
     2 flag_smcdamage,greek_xi,k_ecdamage,k_smcdamage,t_ecdamage,
     3 t_smcdamage,zeta_smcdegrade_c,zeta_collagendegrade_c,
     4 ec_power_coeff,smc_power_coeff,y_smcdamage,
     5 flag_tidentify_smcdamage
cRead input parameters of the Growth and Remodeling from a file
      OPEN(1,FILE='GnRinput',STATUS='OLD',FORM='FORMATTED')
      READ (1,*) ept ! evlauation count
c      READ (1,*) flagfile_exist
      READ (1,*) comp_t ! time of computation
      READ (1,*) delta_P! pressure perturbation
      READ (1,*) delta_Q ! flow perturbation
      READ (1,*) K_c1 ! gain parameter collagen - sigma
      READ (1,*) K_c2 ! gain parameter collagen - tau
      READ (1,*) K_m1 ! gain parameter SMC - sigma
      READ (1,*) K_m2 ! gain parameter SMC - tau
      READ (1,*) C_ratio ! scaling constant for constrictor to dilator
      READ (1,*) G_ch ! prestretch collagen 
      READ (1,*) G_mh ! prestretch muscle
      READ (1,*) G_et ! prestretch elastin - theta
      READ (1,*) G_ez ! prestretch elastin - z
      READ (1,*) kq_m ! Half life muscle
      READ (1,*) kq_c ! Half life collagen
      CLOSE(1)
! check if the flag for flagged simulations is turned on
c      IF (flagfile_exist .EQ. 1.0) THEN
      OPEN(1,FILE='GnRflaggedinputs',STATUS='OLD',FORM='FORMATTED')
      READ (1,*) ept_flag
      READ (1,*) flag_tension_degrade
      READ (1,*) flag_ecdamage
      READ (1,*) flag_smcdamage
      READ (1,*) zeta_smcdegrade_c      ! smc-critical tension length
      READ (1,*) zeta_collagendegrade_c ! collagen-criticaltensionlength
      READ (1,*) k_ecdamage             ! ecdamage-rate parameter
      READ (1,*) k_smcdamage		! smcdamage- rate parameter
      READ (1,*) t_ecdamage             ! time of insult - ec
c      READ (1,*) t_smcdamage		! time of insult - smc
      READ (1,*) ec_power_coeff		! ec power coefficient
      READ (1,*) smc_power_coeff	! smc power coefficient
      CLOSE(1)
c      END IF

c ~~~~~~~~~~~~~~~~~~~Definitions~~~~~~~~~~~~~~~~~~~~~~~~~~
c number of collagen fiber families
      nfiberfly=4
c Viscosity(g/cms)
      mu=0.037
c Density(g/cm^3)
      rho_s=1.050
c Time stepping parameters      
      n_t_step=20.0
      dt=1.0/n_t_step
      num_DL=n_t_step*age_max+1.0
      num_t=comp_t*n_t_step+1.0
      PRINT *, 'num_t is ', num_t
c Newton Raphson parameter.
      beta=0.3

c~~~~~~~~~~~~~~~~~~ALLOCATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ALLOCATE(c_frac(nfiberfly),alpha_ckh(nfiberfly),
     2 alpha_tau(nfiberfly),c_c2(nfiberfly),c_c3(nfiberfly),
     2 M_ckh(nfiberfly),M_ck(nfiberfly))
      ALLOCATE(str_ch(nfiberfly),sigma_ck(nfiberfly),DQ_m(num_DL),
     2 DQ_c(num_DL),Da(num_t),Dm_c(num_t,nfiberfly),Dm_m(num_t),
     2 Dalpha(num_t),DQ2_c(num_t,nfiberfly),DQ2_m(num_t))
      ALLOCATE(Lc_k_tau(nfiberfly),Lc_k(nfiberfly),Lc_kn(nfiberfly),
     2 phi_ck(nfiberfly),nc_stress(nfiberfly))
      c_frac=(/0.1,0.1,0.4,0.4/)
      alpha_ckh=(/ 0.0,90.0,45.0,135.0/)*pi/180.0


      ALLOCATE (dWck_dLn(nfiberfly),alpha_ck(nfiberfly),
     2 m_basal_ck(nfiberfly),mp_c(nfiberfly),dWcdLn_h(nfiberfly))
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c Transmural Pressure mmHgX1333.322=dyn/cm^2      
      P_h=5.0*1333.22
c Vessel initial radius
      a_M=0.2 ! L: this is in cm
c Homeostatic circumferential stree
      sigma_h=18.0*1000.0*10.0 ! L: kPa*1000 = Pa. Pa*10 = dyn/cm^2
c Homeostatic wall shear stress (dyn/cm^2)
      tau_wh=6
c Calculate homeostatic Flow and thickness
      Q_Mh=tau_wh*pi*(a_M**3.0)/(4.0*mu)
      h_h=P_h*a_M/sigma_h
cFlags!
      n_c=1.0
      n_smc=1.0
      f_flag=1.0
      ax=1.0
      flag_tidentify_smcdamage=0.0
cMass fractions 
      phi_c=0.42
      phi_ck=phi_c*c_frac
      phi_e=0.10
      phi_m=1.00-phi_c-phi_e

c Homeostatic masses
      M_h=h_h*rho_s
      M_ckh=M_h*phi_ck
      M_mh= phi_m*M_h
      M_eh= phi_e*M_h

c Homeostatic stress values - note sigma_mh=passive+active
      sigma_ch=sigma_h
      sigma_mh=sigma_h
      sigma_eh=(sigma_h-phi_c*sigma_ch-phi_m*sigma_mh)/phi_e
c      sigma_eh=sigma_h

c Muscle activation parameters
      Lambda_M=1.35
      Lambda_0=0.50
      C_basal=0.68
      T_S0=C_ratio*C_basal
c Passive response parameters and free parameter calculation
      c_m2=1.0*1000.0*10.0 ! *1.0kPa*1000 = Pa. Pa*10 = dyn/cm^2
      c_m3=0.05025
      c_c2(1)=24022.65/phi_ck(1) ! axial
      c_c2(2)=1000/phi_ck(2) ! circumferential
      c_c2(3)=1.0
      c_c2(4)=1.0
      c_c3(1)=0.1 ! L: c2 axial
      c_c3(2)=0.05025 ! L: c2 circumferential
      c_c3(3)=1.035 ! L: c2 diagonal
      c_c3(4)=1.035 ! L: c2 diagonal
     
      c_e=sigma_eh/((G_et**2.0)-(((G_et**2.0)*(G_ez**2.0))**(-1)))
      sigma_m_pas_h=c_m2*((G_mh**2.0)*(G_mh**(2.0)-1.0))*
     2 EXP(c_m3*((G_mh**(2.0)-1.0)**2.0))  
      sigma_act_h=sigma_mh-sigma_m_pas_h
      temp=1.0-((Lambda_M-1.0)/(Lambda_M-Lambda_0))**2.0
      T_M=sigma_act_h/((1.0-EXP(-(C_basal**2.0)))*temp)
      DO i=1, nfiberfly
      str_ch(i)=c_frac(i)*c_c2(i)*(G_ch**2.0)*((G_ch**2.0)-1.0)*
     2 EXP(c_c3(i)*((G_ch**2.0)-1.0)**2.0)*(SIN(alpha_ckh(i)))**2.0
      ENDDO   
     
      c_c2(4)=(sigma_ch-str_ch(1)-str_ch(2))/(2.0*str_ch(4)) !L: this should be = 71
      c_c2(3)=c_c2(4) !L: this should be = 71
      c_c2_c=c_m2
      c_c3_c=c_m3
      c_m2_c=c_m2
      c_m3_c=c_m3
      
c Check if the mixture is at homeostatsis at t=0
      dwdLt_c=0.0
      DO i=1, nfiberfly
       dwdLt_c=dwdLt_c+phi_ck(i)*h_h*dWkdLn(G_ch,c_c2(i),c_c3(i),
     2 c_c2_c,c_c3_c)*G_ch*(SIN(alpha_ckh(i)))**(2.0)
      ENDDO
      dwdLt_m=h_h*phi_m*dWmdLn(G_mh,c_m2,c_m3,c_m2_c,c_m3_c)*G_mh
      dwdLt_e= h_h*phi_e*dWedLn(G_et,G_ez,ax,c_e)*G_et
      dwdLt=dwdLt_c+dwdLt_m+dwdLt_e
      T_act=sigma_act_h*phi_m*h_h
      F=dwdLt+T_act-P_h*a_M
      PRINT *,'Deviation from homeostasis at t=0 is', F      
c Find a stretch to satisfy:
c  dWkdLn(Ln)/dWkdLn(G_ck) = yielding criteria (about 10)
c  dWmdLn(Ln)/dWndLn(G_mh) = yielding criteria
c  (Newton-Raptson method)
       y_Lkn=3.0
       y_Lmn=3.0
       y_smcdamage=2.0
c Rate parameter of active VSM tone in the ref. config.
       k_act=1.0/20.0
c EC damage coefficients
c       ec_power_coeff=0.83! powercoeff for ec damage function, coming from outside the code
c       t_ecdamage=0.001 coming from outside the code
c       k_ecdamage=0.30 ! coming from outside the code
       greek_xi=1.0
cSurvival Fractions Q(t)
       DQ_m(1)=0.0
       DQ_c(1)=0.0
       
       DO i=2,num_DL
        t=dt*(i-1)
        DQ_m(i)=DQ_m(i-1)+0.5*dt*(q_i(t-dt,kq_m)+q_i(t,kq_m))
        DQ_c(i)=DQ_c(i-1)+0.5*dt*(q_i(t-dt,kq_c)+q_i(t,kq_c))
       END DO
       
       mean_age_m=DQ_m(num_DL)
       mean_age_c=DQ_c(num_DL)
          
       DO i=1,num_DL
        DQ_m(i)=(mean_age_m-DQ_m(i))/mean_age_m
        DQ_c(i)=(mean_age_c-DQ_c(i))/mean_age_c
       END DO
      
       m_basal_ck=M_ckh/mean_age_c
       m_basal_m=M_mh/mean_age_m   

c      PRINT *, 'm_basal_ck is ' , m_basal_ck
c       PRINT*, 'm_basal_m is ', m_basal_m
c Data for nn=1 at t=0
       Da(1)=a_M
c Reference radius for active SMT
       a_act_p=a_M
c da_act/dt
       da_act_p=0.0
c Rate of mass production of collagen(4 families) and muscle
       Dm_c(1,:)=m_basal_ck
       Dm_m(1)=m_basal_m       
cFiber angle of newly produced collagen fibers
       Dalpha(1)=alpha_ckh(3)
       DO i=1,num_t
       DQ2_c(i,:)=(/1.0,1.0,1.0,1.0/)
       DQ2_m(i)=1.0
       ENDDO
       sigma_ck=(/sigma_ch,sigma_ch,sigma_ch,sigma_ch/)

c~~~~~~~~~~ Begin Newton-Raphson and March through time~~~~~~~~~~~~~~
      t_gradual=30.0
      DO n_t=2,num_t
      t=(n_t-1.0)*dt
      
c Predictor Step
      IF (t .GT. 2.0) THEN
         IF(t.LE. t_gradual) THEN
          P= P_h+P_h*(delta_P-1.0)*(t-2.0)/(t_gradual-2.0)
          Q_M = Q_Mh+Q_Mh*(delta_Q-1.0)*(t-2.0)/(t_gradual-2.0)
        ELSE
            P=P_h*delta_P;
            Q_M=Q_Mh*delta_Q;
        END IF
      ELSE
         P=P_h
         Q_M=Q_Mh
      END IF

      IF (n_t .GT. 3.0) THEN
         IF ((ABS(Da(n_t-1)-Da(n_t-2))/Da(n_t-1)).LT. 0.1) THEN
         Da(n_t)=2.0*Da(n_t-1)-Da(n_t-2)
         Dm_c(n_t,:)=2.0*Dm_c(n_t-1,:)-Dm_c(n_t-2,:)
         Dm_m(n_t)=2.0*Dm_m(n_t-1)-Dm_m(n_t-2)
         ELSE 
         Da(n_t)=Da(n_t-1)
         Dm_c(n_t,:)=Dm_c(n_t-1,:)
         Dm_m(n_t)=Dm_m(n_t-1)
         END IF
      ELSE
        Da(n_t)=Da(n_t-1.0)
        Dm_c(n_t,:)=Dm_c(n_t-1.0,:)
        Dm_m(n_t)=Dm_m(n_t-1.0)      
      END IF
      a_t=Da(n_t)
c Predictor for a_act
c use eq 16 Baek 2007 cerebral vasospasm II
c da_act/dt=Kact(a(t)-a_act) & a_act=a_act_p+dt*.5*da_act_p
      a_act=a_act_p+0.5*dt*(da_act_p+k_act*(a_t-(a_act_p+dt*da_act_p)))
       
      tol_m=100.0
      tau_w=tau_wh
      num_it1=0
      DO WHILE((tol_m.GT.Max_error_2).AND.(num_it1.LT.Max_it1))
       num_it1=num_it1+1
       tol_a=100.0
       num_it2=0
       DO WHILE((tol_a.GT.Max_error).AND.(num_it2.LT.Max_it2))
       num_it2=num_it2+1
c Current wall shear stress
       tau_w=4.0*mu*Q_M/(pi*(a_t**3.0))
       T_S=T_S0
       IF (flag_ecdamage .EQ. 1.0) THEN
c Begin EC damage modification
        IF (t .LT. t_ecdamage) THEN
        C_t=C_basal-T_S*(tau_w/tau_wh-1.0)
	ELSE
	greek_xi= 1.0-k_ecdamage*((t-t_ecdamage)**ec_power_coeff)*
     2	EXP(1.0-k_ecdamage*t+k_ecdamage*t_ecdamage)
        C_t=(C_basal-T_S*(tau_w/tau_wh-1.0))*greek_xi
	END IF
c End EC damage modification
       ELSE
        C_t=C_basal-T_S*(tau_w/tau_wh-1.0)
       END IF
       dn_tau=(tau_w/tau_wh-1.0)
       IF (C_t .LT. 0.0) THEN
        C_t=0.0
       END IF
c 2D Stretch - constant axial & changing circumferential stretch
       L_z=1.0
       L_t=a_t/a_M
c Alignment of new collagen- defined by stretches(not stresses!)
       Dalpha(n_t)=ATAN(L_t/L_z*TAN(alpha_ckh(3)))
      
       dwdLt_c     = 0.0
       dwdLz_c     = 0.0
       dwdLt_m     = 0.0
       ddwddLt_c   = 0.0
       ddwddLt_m   = 0.0

       M_ck     = (/0.0, 0.0, 0.0, 0.0/)
       dWck_dLn = (/0.0, 0.0, 0.0, 0.0/)
       M_m      = 0.0
      
       IF(n_t.LE.num_DL) THEN

c############BEGIN: Calc initial stress in collagen###############
      M_ck=M_ckh*DQ_c(n_t)*DQ2_c(1,:)           
      Lc_k=SQRT((L_t**2.0)*((SIN(alpha_ckh))**2.0)+(L_z**2.0)*
     2 ((COS(alpha_ckh))**2.0))
   
      DO i=1,nfiberfly
        sQ_ck=DQ_c(n_t)*DQ2_c(1,i)
        dLndLt=G_ch*L_t*((SIN(alpha_ckh(i)))**2.0)/Lc_k(i)
        dLndLz=G_ch*L_z*((COS(alpha_ckh(i)))**2.0)/Lc_k(i)
        ddLnddLt=G_ch*((L_z*SIN(alpha_ckh(i))*
     2 COS(alpha_ckh(i)))**2.0)/(Lc_k(i)**3.0)
        temp=Lc_k(i)*G_ch
       dWck_dLn(i)=dWkdLn(temp,c_c2(i),c_c3(i),c_c2_c,c_c3_c)
       dwdLt_c=dwdLt_c+(M_ckh(i)/rho_s)*sQ_ck*dWck_dLn(i)*dLndLt
       dwdLz_c=dwdLz_c+(M_ckh(i)/rho_s)*sQ_ck*dWck_dLn(i)*dLndLz
       ddwddLt_c=ddwddLt_c+(M_ckh(i)/rho_s)*sQ_ck*
     2 (ddWkddLn(temp,c_c2(i),c_c3(i),c_c2_c,c_c3_c)*
     2 (dLndLt**2.0)+dWck_dLn(i)*ddLnddLt)
       END DO
c##############END: Calc initial stress in collagen#################

c##############BEGIN: Calc initial stress in SMC####################
      M_m =M_mh*DQ_m(n_t)*DQ2_m(1)
      sQ_m=DQ_m(n_t)*DQ2_m(1)
      Lm_n=L_t*G_mh
      dWm_dLn=dWmdLn(Lm_n,c_m2,c_m3,c_m2_c,c_m3_c)
      dwdLt_m=dwdLt_m+(M_mh/rho_s)*sQ_m*dWm_dLn*G_mh
      ddwddLt_m=ddwddLt_m+(M_mh/rho_s)*sQ_m*ddWmddLn(Lm_n,c_m2,c_m3,
     2 c_m2_c,c_m3_c)*(G_mh**2.0)
c##############END: Calc initial stress in SMC#################
       END IF
c CHECK IF THE STRETCH ON SMC IS > ULTIMATE STRETCH AND SET THE FLAG FOR DAMAGE
c       PRINT*,Lm_n
       IF (Lm_n.GE.y_smcdamage.AND.flag_tidentify_smcdamage.EQ.0.0) THEN
       flag_tidentify_smcdamage=1.0
       t_smcdamage=t
       PRINT*, '~~~ALERT~~~: MUSCLE DAMAGE'
       PRINT*, flag_tidentify_smcdamage,t_smcdamage
       END IF

       IF(n_t .LE. num_DL) THEN
         tn0=1
       ELSE
         tn0=n_t-num_DL+1
       END IF
c********************BEGIN: DO n_tau=tn0,n_t*****************
      DO n_tau=tn0,n_t
       IF((n_tau .EQ. tn0).OR.(n_tau .EQ. n_t)) THEN
        wt=0.5*dt
       ELSE
        wt=dt
       END IF
cBegin calculate stress in collagen
      alpha_tau=(/alpha_ckh(1),alpha_ckh(2),Dalpha(n_tau),
     2 2.0*pi-Dalpha(n_tau)/)
      Lt_tau=Da(n_tau)/a_M
      Lz_tau=1.0
      Lc_k_tau=SQRT((Lt_tau**2.0)*((SIN(alpha_tau))**2.0)+
     2 (Lz_tau**2.0)*((COS(alpha_tau))**2.0)) 
       Lc_k=SQRT((L_t**2.0)*((SIN(alpha_tau))**2.0)+
     2 (L_z**2.0)*((COS(alpha_tau))**2.0)) 
       Lc_kn=(G_ch*Lc_k)/Lc_k_tau
      
      DO i=1,nfiberfly
       mc_tau=Dm_c(n_tau,i)
       IF(Lc_kn(i) .LE. y_Lkn) THEN
       
       sq_ck=q_i((n_t-n_tau)*dt,kq_c)*DQ2_c(n_tau,i)
       M_ck(i)=M_ck(i)+mc_tau*sq_ck*wt
       dLndLt=(G_ch/Lc_k_tau(i))*L_t*((SIN(alpha_tau(i)))**2.0)/
     2 Lc_k(i)
       dLndLz=(G_ch/Lc_k_tau(i))*L_z*((COS(alpha_tau(i)))**2.0)/
     2 Lc_k(i)
      ddLnddLt=(G_ch/Lc_k_tau(i))*((L_z*SIN(alpha_tau(i))*
     2 COS(alpha_tau(i)))**2.0)/(Lc_k(i)**(3.0))
      dWck_dLn(i)=dWkdLn(Lc_kn(i),c_c2(i),c_c3(i),c_c2_c,c_c3_c)
      dwdLt_c=dwdLt_c+(mc_tau/rho_s)*sq_ck*dWck_dLn(i)*dLndLt*wt
      dwdLz_c=dwdLz_c+(mc_tau/rho_s)*sq_ck*dWck_dLn(i)*dLndLz*wt
      ddwddLt_c=ddwddLt_c+(mc_tau/rho_s)*sq_ck*
     2 (ddWkddLn(Lc_kn(i),c_c2(i),c_c3(i),c_c2_c,c_c3_c)*
     2 ((dLndLt)**2.0)+dWck_dLn(i)*ddLnddLt)*wt
     
       END IF
      END DO
cEND calculate stress in collagen
cBEGIN calculate stress in SMC
      mm_tau=Dm_m(n_tau)
      Lm_n=G_mh*L_t/Lt_tau
      
      IF(Lm_n.LE.y_Lmn) THEN

c	IF((Lm_n.GE.y_smcdamage).AND.(flag_smcdamage.EQ.1.0).AND.
c     2    (flag_tidentify_smcdamage.EQ.1.0)) THEN
c       sq_m=q_i((n_t-n_tau)*dt,kq_m)*DQ2_m(n_tau)!*
c       M_m=M_m+mm_tau*sq_m*wt
c	ELSE
      sq_m=q_i((n_t-n_tau)*dt,kq_m)*DQ2_m(n_tau)
      M_m=M_m+mm_tau*sq_m*wt
c         END IF

      dWm_dLn=dWmdLn(Lm_n,c_m2,c_m3,c_m2_c,c_m3_c)
      dwdLt_m=dwdLt_m+(mm_tau/rho_s)*sq_m*dWm_dLn*(G_mh/Lt_tau)*wt
      ddwddLt_m=ddwddLt_m+(mm_tau/rho_s)*sq_m*
     2 ddWmddLn(Lm_n,c_m2,c_m3,c_m2_c,c_m3_c)*
     2 ((G_mh/Lt_tau)**2.0)*wt
      END IF
  
cEND calculate stress in SMC

      END DO
c********************END: DO n_tau=tn0,n_t*****************
      M_c=M_ck(1)+M_ck(2)+M_ck(3)+M_ck(4)


cBEGIN calculate stress in Elastin
       M_e=M_eh
       dwdLt_e=(M_e/rho_s)*dWedLn(L_t*G_et,L_z*G_ez,ax,c_e)*G_et
       ddwddLt_e=(M_e/rho_s)*ddWeddLn(L_t*G_et,L_z*G_ez,ax,c_e)*
     2 (G_et**2.0)
cEND calculate stress in Elastin

cBEGIN: Active Smooth Muscle contribution
      L_m_act=a_t/a_act
      J=L_t*L_z
      temp=((Lambda_M-L_m_act)/(Lambda_M-Lambda_0))**2.0
      T_act=(T_M*M_m/(rho_s*J))*(1.0-EXP(-(C_t**2.0)))*L_m_act*
     2 (1.0-temp) 
c      PRINT*, T_act
c Active muscle can generate only tension
      IF (T_act .LT. 0.0) THEN
         T_act=0.0
      END IF
c recall C_t=C_basal-T_S*(tau_w/tau_wh-1)   
      dc_tda=T_S*12.0*mu*Q_M/(tau_wh*pi*(a_t**4.0))
      temp1=-T_M*M_m*(L_z/(rho_s*(J**2.0)*a_t))*
     2 (1.0-EXP(-(C_t**2.0)))*L_m_act*(1.0-temp)
      temp2=(T_M*M_m/(rho_s*J))*((1.0-EXP(-(C_t**2.0)))/a_act)*
     2 (1.0-temp)
      temp3=(T_M*M_m/(rho_s*J))*(1.0-EXP(-(C_t**2.0)))*L_m_act*
     2 ((Lambda_M - L_m_act)/(((Lambda_M - Lambda_0)**2.0)*a_act))
      temp4=(T_M*M_m/(rho_s*J))*2.0*C_t*dC_tda*EXP(-(C_t**2.0))*
     2 L_m_act*((Lambda_M - L_m_act)/(((Lambda_M - Lambda_0)**2.0)*
     2 a_act))
c       PRINT*, temp1,temp2,temp3,temp4
      dT_actda=temp1+temp2+temp3+temp4
c      PRINT*, dT_actda
c END: Active Smooth Muscle Contribution

      total_M=M_c+M_e+M_m
      dwdLt=dwdLt_c+dwdLt_m+dwdLt_e
      ddwddLt=ddwddLt_c+ddwddLt_m+ddwddLt_e
c      PRINT*,'dwdlt values: ', dwdLt,dwdLt_c,dwdLt_m,dwdLt_e
c      PRINT*,'ddwddlt values: ',ddwddLt,ddwddLt_c,ddwddLt_m,ddwddLt_e
c       PRINT*, dwdLt/l_z,T_act, P, a_t,P*a_t
       F_a=(dwdLt/L_z)+T_act-(P*a_t)
c       PRINT*, 'Fa 1',F_a
cNote: -dPda*a_t has been dropped from dF_ada!
      dF_ada=(1.0/L_z)*ddwddLt*(1.0/a_M)+dT_actda-P


      a_t=a_t-beta*(F_a/dF_ada)
      tol_a=sqrt(((a_t-Da(n_t))**2.0)/(Da(n_t)**2.0))
      Da(n_t)=a_t
      a_act=1.0/(1.0+0.5*dt*k_act)*(a_act_p+0.5*dt*(da_act_p+k_act*a_t))
      L_t=a_t/a_M
      Dalpha(n_t)=ATAN((L_t/L_z)*TAN(alpha_ckh(3)))

      END DO
c~~End of DO WHILE((tol_a.GT.Max_error).AND.(num_it2.LT.Max_it))~~~


c Adjust beta if Newton Raphson does not converge
      IF(num_it2 .EQ. Max_it2) THEN
       beta=0.1
      END IF
      
      T_t_c=(rho_s*L_t*L_z/M_c)*(dwdLt_c)/L_z
      T_z_c=(rho_s*L_t*L_z/M_c)*(dwdLz_c)/L_t
      T_t_m=(rho_s*L_t*L_z/M_m)*((dwdLt_m/L_z)+T_act)
c       PRINT*, T_t_c,T_z_c,T_t_m
      temp1=0.0
      temp2=pi/2.0
      alpha_ck=(/temp1,temp2,Dalpha(n_t),2.0*pi-Dalpha(n_t)/)
c      PRINT*,'alpha_ck',alpha_ck,T_z_c,T_t_c
      
      DO i=1,nfiberfly
       sigma_ck(i)=sqrt(((T_z_c*COS(alpha_ck(i)))**2.0)+
     2 ((T_t_c*SIN(alpha_ck(i)))**2.0))
      nc_stress(i)=(sigma_ck(i)/sigma_ch-1.0)
      END DO
c      PRINT*, 'sigma_ck is ',sigma_ck
      nm_stress=(T_t_m/sigma_mh-1.0)
      dn_C=(C_t/C_basal-1.0)
      temp1=rm_c(n_c,m_basal_ck(1),nc_stress(1),dn_tau,K_c1,K_c2)
      temp2=rm_c(n_c,m_basal_ck(2),nc_stress(2),dn_tau,K_c1,K_c2)
      temp3=rm_c(n_c,m_basal_ck(3),nc_stress(3),dn_tau,K_c1,K_c2)
      temp4=rm_c(n_c,m_basal_ck(4),nc_stress(4),dn_tau,K_c1,K_c2)
      mp_c=(/temp1,temp2,temp3,temp4/)
      m_m=rm_m(n_smc,m_basal_m,nm_stress,dn_tau,K_m1,K_m2)
       
c make sure the production rates are not negative
      DO i=1,nfiberfly
       IF(mp_c(i) .LT. 0.0) THEN
       mp_c(i)=0.0
       END IF
      END DO

      IF(m_m.LT. 0.0) THEN
         m_m=0.0
      END IF
      temp1=((mp_c(1)-Dm_c(n_t,1))**2.0)+((mp_c(2)-Dm_c(n_t,2))**2.0)+
     2 ((mp_c(3)-Dm_c(n_t,3))**2.0)+((mp_c(4)-Dm_c(n_t,4))**2.0)
      temp2=(m_m-Dm_m(n_t))**2.0
      temp3=(mp_c(1)**2.0)+(mp_c(2)**2.0)+(mp_c(3)**2.0)+(mp_c(4)**2.0)
      temp4=m_m**2.0 
      tol_m=SQRT((temp1+temp2)/(temp3+m_m**2.0))
      Dm_c(n_t,:)=mp_c
      Dm_m(n_t)=m_m
      
      END DO
c~~End of  DO WHILE((tol_m.GT.Max_error_2).AND.(num_it1.LT.Max_it))~~

      a_act_p=a_act
      da_act_p=k_act*(a_t-a_act)
      h=total_M/(rho_s*L_t*L_z) 
      mp_c_file=mp_c/m_basal_ck
      m_m_file=m_m/m_basal_m
      tau_w_file=tau_w/tau_wh
      tau_w_file2=tau_w/tau_wh-1.0
c      OPEN(9,FILE='dmglog',FORM='FORMATTED',ACCESS='APPEND')
c      WRITE (9,111) t,Lm_n,y_smcdamage,
c     2 flag_smcdamage,flag_tidentify_smcdamage,flag_smcdamage.EQ.1.0,
c     3 flag_tidentify_smcdamage.EQ.1.0,Lm_n.GE.y_smcdamage
c      CLOSE(9)
c111   format (8f15.9,3L5)      
      WRITE (filename1,'(A15,I0)') "GnRoutput1grad_", ept    
      OPEN(10,FILE=filename1,FORM='FORMATTED',ACCESS='APPEND')
      WRITE (10,'(15f20.9)') t,P,a_t,h,a_act,Q_M,C_t,T_act,
     2 Dalpha(n_t)
      CLOSE(10)
      
      WRITE (filename2,'(A15,I0)') "GnRoutput2grad_", ept
      OPEN(11,FILE=filename2,FORM='FORMATTED',ACCESS='APPEND')
      WRITE (11,'(15f20.9)') t,mp_c_file,m_m_file,M_m,M_e,M_ck
      CLOSE(11)

c      PRINT*, 'm_basal_ck is ', m_basal_ck
c      PRINT*, 'm_basal_m is ', m_basal_m
      WRITE (filename3,'(A15,I0)') "GnRoutput3grad_", ept
      OPEN(12,FILE=filename3,FORM='FORMATTED',ACCESS='APPEND')
      WRITE (12,'(15f20.9)') t,tau_w_file,nc_stress,nm_stress,
     2 m_basal_ck,m_basal_m,greek_xi,Lm_n
      CLOSE(12)


c Update DQ2_c and DQ2_m - tension dependent degradation
      L_t=a_t/a_M
      IF(n_t .LE. num_DL) THEN
      tn0=1
      ELSE
      tn0=n_t-num_DL+1
      END IF

      DO i=1,nfiberfly
      dWcdLn_h(i)=dwkdLn(G_ch,c_c2(i),c_c3(i),c_c2_c,c_c3_c)
      END DO
      dWmdLn_h=dWmdLn(G_mh,c_m2,c_m3,c_m2_c,c_m3_c)
      
      DO n_tau=tn0,n_t
      alpha_tau=(/alpha_ckh(1),alpha_ckh(2),Dalpha(n_tau),2.0*pi-
     2 Dalpha(n_tau)/)
      Lt_tau=Da(n_tau)/a_M
      Lz_tau=1.0
      Lc_k_tau=SQRT((Lt_tau**2.0)*((SIN(alpha_tau))**2.0)+
     2 (Lz_tau**2.0)*((COS(alpha_tau))**2.0)) 
      Lc_k=SQRT((L_t**2.0)*((SIN(alpha_tau))**2.0)+
     2 (L_z**2.0)*((COS(alpha_tau))**2.0)) 
      Lc_kn=(G_ch*Lc_k)/Lc_k_tau
      DO i=1,nfiberfly
       beta_c=dWkdLn(Lc_kn(i),c_c2(i),c_c3(i),c_c2_c,c_c3_c)/
     2 dWcdLn_h(i)
       DQ2_c(n_tau,i)=EXP(-1.0*f_beta(beta_c,kq_c)*dt)
c note :*DQ2_c(n_tau,i) is commented out in matlab code
      END DO
      Lm_n=G_mh*(L_t/Lt_tau)
      beta_m=dWmdLn(Lm_n,c_m2,c_m3,c_m2_c,c_m3_c)/dWmdLn_h

       IF((flag_smcdamage.EQ.1.0).AND.
     2 (flag_tidentify_smcdamage.EQ.1.0)) THEN
      DQ2_m(n_tau)=EXP(-1.0*f_beta(beta_m,kq_m)*dt)*
     2 dmg_smc(t,t_smcdamage,k_smcdamage,smc_power_coeff)
        ELSE
         DQ2_m(n_tau)=EXP(-1.0*f_beta(beta_m,kq_m)*dt)  
	END IF
      END DO
      IF(ISNAN(h)) THEN
      GOTO 21
      END IF
     
      IF(sigma_act_h .LT. 0.0) THEN
       tau_wh=1000.0
       nm_stress=10000.0
       GOTO 21
      END IF
c      PRINT *, 'n_t is ' , n_t

      END DO
c~~~~~~~~End of DO n_t=2,num_t ~~~~~~~~~~~~~~~~~~~~~~~~~

      OPEN(1,FILE='C_fn',FORM='FORMATTED',STATUS='REPLACE')
      WRITE (1,'(15f20.9)') a_t,h,tau_w_file2,nc_stress,nm_stress    
      CLOSE(1)

 21     STOP
      END PROGRAM main_parameterstudy
!########      EXTERNAL FUNCTIONS ########################
      REAL*8 FUNCTION q_i(s,kq)
      USE COM_GnR
      IMPLICIT NONE

      REAL*8 s,kq
      q_i=EXP(-kq*s)
      RETURN
      END FUNCTION q_i


      REAL*8 FUNCTION dmg_smc(t,t_smcdamage,k_smcdamage,
     2 smc_power_coeff)
      USE COM_GnR
      IMPLICIT NONE

      REAL*8 t,t_smcdamage,k_smcdamage,smc_power_coeff,ctemp
      ctemp=k_smcdamage*(t-t_smcdamage)**(smc_power_coeff)
      dmg_smc= 1.0-ctemp*EXP(1.0-k_smcdamage*t+k_smcdamage*t_smcdamage)
c      PRINT *, t,t_smcdamage,ctemp,dmg_smc
      RETURN
      END FUNCTION dmg_smc


      REAL*8 FUNCTION f_beta(beta,kq)
      USE COM_GnR
      IMPLICIT NONE

      REAL*8 beta,kq
      f_beta=ABS(beta-1)*kq
      RETURN
      END FUNCTION f_beta

      REAL*8 FUNCTION rm_m(n_SMC,m_basal,dn_stress,dn_tau,Zeta_m,Km2)
      IMPLICIT NONE

      REAL*8 n_SMC,m_basal,dn_stress,dn_tau,Zeta_m,Km2
      rm_m=n_SMC*m_basal*(Zeta_m*dn_stress-1.0*Km2*dn_tau+1.0)
      RETURN
      END FUNCTION rm_m

      REAL*8 FUNCTION rm_c(n_c,m_basal,dn_stress,dn_tau,Zeta_c,Kc2)
      IMPLICIT NONE
      REAL*8 n_c,m_basal,dn_stress,dn_tau,Zeta_c,Kc2
      rm_c=n_c*m_basal*(Zeta_c*dn_stress-1.0*Kc2*dn_tau+1.0)
      RETURN
      END FUNCTION rm_c

      REAL*8 FUNCTION dWedLn(Ln_t,Ln_z,ax,c_e)
      USE COM_GnR
      IMPLICIT NONE
      REAL*8 Ln_t,Ln_z,ax,c_e
c axis=1 is circumferential, axis =2 is axial
      IF (ax .EQ. 1.0) THEN
      dWedLn=c_e*(Ln_t-1.0/(Ln_t**(3.0)*Ln_z**(2.0)))
      ELSE IF (ax .EQ. 2.0) THEN
      dWedLn=c_e*(Ln_z-1.0/(Ln_z**(3.0)*Ln_t**(2.0)))
      ELSE
      dWedLn=0.0
      PRINT *,'WRONG PARAMETERS IN dWedLnfunction'
      END IF
      RETURN
      END FUNCTION dWedLn
      
      REAL*8 FUNCTION dWkdLn(Ln,c_c2,c_c3,c_c2_c,c_c3_c)
      IMPLICIT NONE
      REAL*8 Ln,c_c2,c_c3,c_c2_c,c_c3_c
      IF(Ln.GE.1) THEN
      dWkdLn=c_c2*Ln*(Ln**(2.0)-1.0)*EXP(c_c3*(Ln**(2.0)-1.0)**(2.0))
      ELSE
      dWkdLn=c_c2_c*Ln*(Ln**(2.0)-1.0)*EXP(c_c3_c*(Ln**(2.0)-1)**(2.0))
      END IF
      RETURN
      END FUNCTION dWkdLn

      REAL*8 FUNCTION dWmdLn(Ln,c_m2,c_m3,c_m2_c,c_m3_c)
      IMPLICIT NONE
      REAL*8 Ln,c_m2,c_m3,c_m2_c,c_m3_c
      IF(Ln.GE.1) THEN
       dWmdLn=c_m2*Ln*(Ln**(2.0)-1.0)*EXP(c_m3*(Ln**(2.0)-1.0)**(2.0))
      ELSE
       dWmdLn=c_m2_c*Ln*(Ln**(2.0)-1.0)*EXP(c_m3_c*(Ln**(2.0)-1.0)
     2 **(2.0))
      END IF
      RETURN
      END FUNCTION dWmdLn

      REAL*8 FUNCTION ddWeddLn(Ln_t,Ln_z,f_flag,c_e)
      USE COM_GnR
      IMPLICIT NONE
      REAL*8 Ln_t,Ln_z,f_flag,c_e
      IF (f_flag.EQ.1.0) THEN
        ddWeddLn=c_e*(1+3.0/(Ln_t**(4.0)*Ln_z**(2.0)))
      ELSE IF (f_flag.EQ.2.0) THEN
        ddWeddLn=c_e*(1+3.0/(Ln_z**(4.0)*Ln_t**(2.0)))
      ELSE IF (f_flag.EQ.3.0) then
         ddWeddLn=c_e*2.0/((Ln_t*Ln_z)**(3.0))
      ELSE
         ddWeddLn=0.0
          PRINT *,'WRONG PARAMETER FOR ddWddLn'
      END IF
      RETURN
      END FUNCTION ddWeddLn

      REAL*8 FUNCTION ddWkddLn(Ln,c_c2,c_c3,c_c2_c,c_c3_c)
      IMPLICIT NONE
      REAL*8 Ln,c_c2,c_c3,c_c2_c,c_c3_c,exp_Q
      
      IF (Ln.GE.1) THEN
        exp_Q=EXP(c_c3*(Ln**(2.0)-1.0)**(2.0))
        ddWkddLn=c_c2*(3.0*Ln**(2.0)-1.0+4.0*c_c3*
     2 ((Ln*(Ln**2.0-1.0))**2.0))*exp_Q
      ELSE 
        exp_Q=EXP(c_c3_c*(Ln**(2.0)-1.0)**(2.0))
        ddWkddLn=c_c2_c*(3.0*Ln**(2.0)-1.0+4.0*c_c3_c*
     2 ((Ln*(Ln**2.0-1.0))**2.0))*exp_Q
      END IF
  
      RETURN
      END FUNCTION ddWkddLn

      REAL*8 FUNCTION ddWmddLn(Ln,c_m2,c_m3,c_m2_c,c_m3_c)
      IMPLICIT NONE
      REAL*8 Ln,c_m2,c_m3,c_m2_c,c_m3_c,exp_Q
      IF (Ln.GE.1) THEN
       exp_Q=EXP(c_m3*(Ln**2.0-1.0)**(2.0))
       ddWmddLn=c_m2*(3.0*Ln**(2.0)-1.0+4.0*c_m3*((Ln*(Ln**2.0-1.0)
     2 )**2.0))*exp_Q
      ELSE
      exp_Q=EXP(c_m3_c*(Ln**(2.0)-1.0)**(2.0))
      ddWmddLn=c_m2_c*(3.0*Ln**(2.0)-1.0+4.0*c_m3_c*(
     2 (Ln*(Ln**(2.0)-1.0))**2.0))*exp_Q
       END IF
      RETURN
      END FUNCTION ddWmddLn

      
