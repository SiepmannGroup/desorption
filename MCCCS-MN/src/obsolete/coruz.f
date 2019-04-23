      function coruz(imolty,jmolty,rho,ibox)
!      function coruz(iunit,rho,ibox)

! ***************************************************
! *** tail-corrections in energy with the zeolite ***
! ***************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc'
!$$$      include 'system.inc'
!$$$      include 'expsix.inc'
!$$$      include 'merck.inc'
!$$$      include 'nsix.inc'
!$$$      include 'external.inc'

      real(KIND=double_precision)::coruz,eps,rci3,rho
      integer(KIND=normal_int)::iunit,ibox
      real(KIND=double_precision)::epsilon2,sigma2
      real(KIND=double_precision)::rci1
      integer(KIND=normal_int)::imolty,jmolty,ii,jj, ntii, ntjj, ntij

! --- note works only for alkanes!!!
!      if (iunit.ne.1) then
!        eps=zeps(1,4)+(iunit-2)*zeps(2,4)+zeps(3,4)
!      else
!        eps=zeps(1,4)
!      end if
!      rci3=zsig2(1,4)**(3.d0/2.d0)/rcut(ibox)**3 
!      coruz=8.*onepi*eps*rho*(rci3*rci3*rci3/9.-rci3/3.)

      coruz=0.
      do ii = 1, nunit(imolty) 
         ntii = ntype(imolty,ii)
!         do jj = 1, zntype
            ntjj = ztype(jmolty)
            ntij = (ntii-1)*nntype + ntjj
            rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
            if ( lexpand(imolty) ) then
               sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
               epsilon2 = sqrt(epsilon(imolty,ii)*epsi(ntjj))
            else
               sigma2 = sig2ij(ntij)
               epsilon2 = epsij(ntij)
            end if
            coruz = coruz +  8.0d0 * onepi * epsilon2 *  sigma2**(1.5d0) *rho *  (rci3 * rci3 * rci3 / 9.0d0 - rci3 / 3.0d0)
!            end do
      end do
      return
      end
