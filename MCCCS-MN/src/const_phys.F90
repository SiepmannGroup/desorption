module const_phys
  use var_type, only: double_precision
  use const_math,only:twopi,fourpi
  implicit none

  real,parameter::h_planck=6.62606957E-34_double_precision

  real,parameter::m_unit=9.10938291E-31_double_precision

  real,parameter::R_gas=8.3144621_double_precision&
   ,N_Avogadro=6.02214129E23_double_precision&
   ,k_B=1.3806488E-23_double_precision&
   ,cal2joule=4.184_double_precision,joule2cal=1.0_double_precision/cal2joule

  real,parameter::MPa2SimUnits=1E-24_double_precision/k_B !< conversion factor for Mpa to simulation unit (including MPa-->Pa, Ang^3-->m^3, J-->K)

  real,parameter::debroglie_factor=1.74582182191543E1_double_precision !<h_planck/sqrt(twopi*1E-3_double_precision/N_Avogadro*k_B)*1E10_double_precision !< in Ang

  real,parameter::e_unit=1.602176565E-19_double_precision&
   ,eps_0=8.854187817E-12_double_precision&
   ,eXV_to_K=e_unit/k_B& !< e times V to Kelvin
   ,qqfact=(e_unit**2)*1E10_double_precision/fourpi/eps_0/k_B !< conversion factor for coulomb interactions (including e-->C, m-->Ang, fourpi*epsilon_0, J-->K)
  !< qqfact = 0.99865377 for bohr/hartree - from pot_KAng.f code
end module const_phys
