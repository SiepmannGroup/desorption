module const_math
  use var_type, only: double_precision
  implicit none

  real,parameter::onepi=3.1415926535897932384626433832795028841971693993751058209749_double_precision&
   ,twopi=6.2831853071795864769252867665590057683943387987502116419498_double_precision&
   ,fourpi=12.566370614359172953850573533118011536788677597500423283899_double_precision&
   ,sqrtpi=1.7724538509055160272981674833411451827975494561223871282138_double_precision

  real,parameter::ln10=2.3025850929940456840179914546843642076011014886287729760333_double_precision

  real,parameter::raddeg=180.0_double_precision/onepi,degrad=onepi/180.0_double_precision

  real,parameter::eps=1.0E-9_double_precision !eps=epsilon(1.0_single_precision)
end module const_math
