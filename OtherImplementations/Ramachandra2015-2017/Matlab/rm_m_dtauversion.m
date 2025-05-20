function y = rm_m_dtauversion(n_SMC, m_basal, dn_stress, dn_tau, Zeta_m, Km2)

y = n_SMC*m_basal*(Zeta_m*dn_stress - Km2*dn_tau + 1);

%y = n_SMC * m_basal * (1 + 3*(1 - exp(-Zeta_m*dn_stress^2) ) + Km2*dn_C);
