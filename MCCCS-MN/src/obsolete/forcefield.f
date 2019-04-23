      subroutine forcefield(rczeo)

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
   
!$$$      include 'zeolite.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'control.inc'

      real(KIND=double_precision)::epsilon,sigma,rczeo
      integer(KIND=normal_int)::         itype,jtype,idz

!     In atomtype.f atoms O of the zeolite have been given id=1
!     Hopefully in some other place the guestmolecules CH4 have been
!     assigned id=2

!      if (zntype.ne.2) call cleanup('** forcefield: ntype ne 2 not allowed **')

!     Force field data:

!     The interaction parameters are in this case calculated from 
!     atomic values

!     Fill the matrices

! Needed are epsilon and sigma such that in Angstroms and kcal/mol:
! Energy = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

      read(25,*)  
      read(25,*)
      idz=1
      do itype = 1,zntype
	 read(25,*) epsilon,sigma
         zeps(itype,idz) = epsilon
         zsig2(itype,idz) = sigma*sigma
      end do

! Write epsilons and sigma's:

      write(16,"(/,' FORCE FIELD DEFINITION:',/, ' ------------------------------------------------',/ ,' Lennard-Jones Parameters (K, Angstrom):')")
      jtype=idz
      do itype = 1,zntype
         write(16,"(  ' epsilon(',i2,',',i2,') = ',f6.1, ' K &  sigma (',i2,',',i2,') = ',f6.1,' Angstrom')") itype,jtype,zeps(itype ,jtype), itype,jtype,sqrt(zsig2(itype,jtype))
      end do
 
! Force field cutoff radius rczeo for interactions between guests
! and zeolite atoms:

      write(16,"(/,' Force field cutoffs: rczeo = ',f6.1,/, '                      rcads = ',f6.1,' Angstrom',/, ' ------------------------------------------------',//)") rczeo
      do itype = 1,zntype
        zrc2(itype,idz)=rczeo**2
      end do
!
! === calculate cutt-off of the potential
!
      if (lshift) then
        write(io_output,*) ' Value at cut-off distance '
        do itype = 1,zntype
          zencut(itype,idz)=4.*zeps(itype,idz)* ( (zsig2(itype,idz)/zrc2(itype,idz))**6 -(zsig2(itype,idz)/zrc2(itype,idz))**3)
          write(io_output,"('    interaction ',i3, '   : ',f8.2,'[K]')") itype ,zencut(itype,idz)
        end do
      end if     

      return
      end

