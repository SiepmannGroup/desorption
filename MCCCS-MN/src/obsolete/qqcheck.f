      subroutine qqcheck(i,ibox,rxuu1,ryuu1,rzuu1)

!     **************************************************************
!     ***  sets up the coulom() array for use with the group based *
!     ***  cutoff for charged molecules     Marcus Martin          *
!     **************************************************************

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

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'poten.inc'
!$$$      include 'system.inc'
!$$$      include 'qqlist.inc'

      integer(KIND=normal_int)::i,ibox,jmolty,j

      real(KIND=double_precision)::rcutsq,rxuu1,ryuu1,rzuu1,rxuij,ryuij,rzuij ,rijsq

      rcutsq = rcut(ibox)*rcut(ibox)

      if ( lpbc ) call setpbc(ibox)

      do j = 1,nchain
         lcoulom(j) = .false.
         jmolty = moltyp(j)

         if ( (nboxi(j) .eq. ibox) .and. (i .ne. j) ) then

! ---  check for the group based qq cutoff
            if ( lelect(jmolty) ) then
               rxuij = rxuu1 - rxu(j,1)
               ryuij = ryuu1 - ryu(j,1)
               rzuij = rzuu1 - rzu(j,1)

! --- minimum image the coulombic bead pair separations ***
               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

               if ( (rijsq .lt. rcutsq) .or. lchgall) then
                  lcoulom(j) = .true.
               end if
            end if

         end if

      end do

      return
      end
