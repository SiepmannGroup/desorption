      function mbessel(z,nu)

      use var_type,only:double_precision
      use const_math,only:onepi
      implicit none

        real(KIND=double_precision)::z,nu
        real(KIND=double_precision)::mbessel
        
!       mbessel = sqrt(onepi/(2.0d0*z))*exp(-z)*
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))

! -- simple form     
        mbessel = sqrt(onepi/(2.0d0*z))*exp(-z)
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))
        end function mbessel

