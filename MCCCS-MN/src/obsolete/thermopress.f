      subroutine thermopress
 
!    *******************************************************************
!    ** calculates the pressure for a configuration according to the  **
!    ** thermodynamic definition. See J. Chem. Phys. Vol.105 P8469.   **
!    ** written in 1998 by Bin Chen.                                  **
!    *******************************************************************
 
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
 
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'control.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'conver.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'blkavg.inc'
!$$$      include 'inputdata.inc'

      logical::lvol,ovrlap
      integer(KIND=normal_int)::i,j,itest,imolty

      integer(KIND=normal_int)::ibox,ic,ncount
      real(KIND=double_precision)::dv,apress
      real(KIND=double_precision)::bxo(2),byo(2),bzo(2)
      real(KIND=double_precision)::volo(2),vboxo(2),dfac(2),voln(2), vintero(2),vtailo(2),vexto(2),velecto(2)
      real(KIND=double_precision)::kxo(vectormax,2),kyo(vectormax,2), kzo(vectormax,2),prefacto(vectormax,2)
      real(KIND=double_precision)::vboxn(2), vintern(2),vtailn(2),vextn(2),velectn(2)
      real(KIND=double_precision)::rxuo(nmax,numax),ryuo(nmax,numax) ,rzuo(nmax,numax)
      real(KIND=double_precision)::volt,expdv,random,df,dx,dy,dz,v,dele, vinter,vtail,vext,vminim,vdum,velect
      real(KIND=double_precision)::xcmo,ycmo,zcmo
      dimension xcmo(nmax),ycmo(nmax),zcmo(nmax)

! --------------------------------------------------------------------


! --- ghost volume move to calculate the pressure
! --- store old box lengths, energy, configuration etc
      do i = 1, 2
         bxo(i)    = boxlx(i)
         byo(i)    = boxly(i)
         if ( lpbcz ) bzo(i)    = boxlz(i)

         if ( lpbcz ) then
            volo(i)   = bxo(i)*byo(i)*bzo(i)
         else
            volo(i)   = bxo(i)*byo(i)
         end if
	 vboxo(i)    = vbox(i)
	 vintero(i)  = vinterb(i)
         vtailo(i)   = vtailb(i)
	 vexto(i)    = vextb(i)  
         velecto(i)  = velectb(i)

! --- store old k vectors and reciprocal sum
         if ( lewald ) then
            call recip(i,vdum,vdum,3)
            ncount = numvect(i)
            do ic = 1,ncount
               kxo(ic,i) = kx(ic,i)
               kyo(ic,i) = ky(ic,i)
               kzo(ic,i) = kz(ic,i)
               prefacto(ic,i) = prefact(ic,i)
            end do
         end if
      end do

      do i = 1, nchain
         imolty = moltyp(i)
         xcmo(i) = xcm(i)
         ycmo(i) = ycm(i)
         if (lpbcz) zcmo(i) = zcm(i)
         do j = 1, nunit(imolty)
            rxuo(i,j) = rxu(i,j)
            ryuo(i,j) = ryu(i,j)
            if ( lpbcz ) rzuo(i,j) = rzu(i,j)
         end do
      end do
      bnpress = bnpress + 1.0d0
      do itest = 1,ndvtry
         dv = dvtry(itest)
         voln(1) = volo(1) + dv
         voln(2) = volo(2) + dv
         if ( lpbcz ) then
            dfac(1) = (voln(1)/volo(1))**(1.0d0/3.0d0)
            dfac(2) = (voln(2)/volo(2))**(1.0d0/3.0d0)
         else
            dfac(1) = sqrt(voln(1)/volo(1))
            dfac(2) = sqrt(voln(2)/volo(2))
         end if            
         boxlx(1) = boxlx(1) * dfac(1)
         boxly(1) = boxly(1) * dfac(1)
         boxlx(2) = boxlx(2) * dfac(2)
         boxly(2) = boxly(2) * dfac(2)
         if ( lpbcz ) then
            boxlz(1) = boxlz(1) * dfac(1)
            boxlz(2) = boxlz(2) * dfac(2)
         end if

         do i = 1, nchain
            
            imolty = moltyp(i)
            ibox = nboxi(i)
            
            df = dfac(ibox) - 1.0d0
            
            dx = xcm(i) * df
            dy = ycm(i) * df
            if ( lpbcz ) dz = zcm(i) * df
            
            xcm(i) = xcm(i) + dx
            ycm(i) = ycm(i) + dy
            if ( lpbcz ) zcm(i) = zcm(i) + dz
            do j = 1, nunit(imolty)
               rxu(i,j) = rxu(i,j) + dx
               ryu(i,j) = ryu(i,j) + dy
               if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
            end do
         end do
         
         lvol = .true.
         do 400 ibox = 1,2
            call sumup( ovrlap, v, vinter,vtail, vdum,vdum, vdum,vdum,vext,velect,vdum, ibox,lvol)
            if ( ovrlap ) goto 500
            vintern(ibox) = vinter
            vtailn(ibox)  = vtail
            vextn(ibox)   = vext  
            velectn(ibox) = velect
            vboxn(ibox)   = vboxo(ibox) +  (vintern(ibox)-vintero(ibox)) + (vextn(ibox)-vexto(ibox)) +  (velectn(ibox)-velecto(ibox))
            apress = exp(-(vboxn(ibox)-vboxo(ibox))*beta)* ((voln(ibox)/volo(ibox))**nchbox(ibox))
            apresscum(ibox,itest) = apresscum(ibox,itest)  + apress

 400     continue
! --- restore old box lengths
 500     do i = 1, 2
            boxlx(i)   = bxo(i)
            boxly(i)   = byo(i)
            
            if ( lpbcz ) boxlz(i)   = bzo(i)
            
            if ( lewald ) then
!     --- restore old k vectors and reciprocal sum and calp
               ncount = numvect(i)
               call recip(i,vdum,vdum,4)
               do ic = 1,ncount
                  kx(ic,i) = kxo(ic,i)
                  ky(ic,i) = kyo(ic,i)
                  kz(ic,i) = kzo(ic,i)
                  prefact(ic,i) = prefacto(ic,i)
               end do
            end if
         end do

         do i = 1, nchain
            imolty = moltyp(i)
            xcm(i) = xcmo(i)
            ycm(i) = ycmo(i)
            if ( lpbcz ) zcm(i) = zcmo(i)
            do j = 1, nunit(imolty)
               rxu(i,j) = rxuo(i,j)
               ryu(i,j) = ryuo(i,j)
               if ( lpbcz ) rzu(i,j) = rzuo(i,j)
            end do
         end do
      end do

!      write(io_output,*) ' press tail' ,  press

      return
      end











