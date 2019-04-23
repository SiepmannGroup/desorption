      subroutine update(nblock,ipos,ibox,value,acmove)
!
! *** this subroutine updates the block averages
!
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
!$$$      include 'system.inc'
!$$$      include 'neigh.inc'
!$$$      include 'neigh2.inc'
!$$$      include 'blkavg.inc'
      integer(KIND=normal_int)::nblock,ipos,ibox
      real(KIND=double_precision)::acmove,dp,dn,value
      
      if (nblock.eq.1) then

! ---     first block

!         write(io,*) 'acmove', acmove       
  
         dn = acmove
         dp = value
         if ( dn .lt. 0.5d0 ) then
            baver(ipos,ibox,nblock) = 0.0d0
         else
            baver(ipos,ibox,nblock) = dp / dn
         end if
         nccold(ipos,ibox) = dn
         bccold(ipos,ibox) = dp
         naccu(ipos,ibox) = naccu(ipos,ibox) + dn
         accum(ipos,ibox) = accum(ipos,ibox) + dp
      else      

! ---   other blocks
!         write(io,*) 'acmove', acmove
         dn = acmove - nccold(ipos,ibox)
         dp = value - bccold(ipos,ibox)
         if ( dn .lt. 0.5d0 ) then
            baver(ipos,ibox,nblock) = 0.0d0
         else
            baver(ipos,ibox,nblock) = dp / dn
         end if
         nccold(ipos,ibox) = acmove
         bccold(ipos,ibox) = value
         naccu(ipos,ibox) = naccu(ipos,ibox) + dn
         accum(ipos,ibox) = accum(ipos,ibox) + dp
      end if
      return
      end



















