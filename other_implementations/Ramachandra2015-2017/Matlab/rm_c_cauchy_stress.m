function y = rm_c_dtauversion(n_c, m_basal,dn_stress, dn_tau, Zeta_c, Kc2)

y = n_c*m_basal*(Zeta_c*dn_stress - Kc2*dn_tau + 1);

%y = n_c * m_basal * (1 + 3*(1 - exp(-Zeta_c*dn_stress^2) ) + Kc2*dn_C);