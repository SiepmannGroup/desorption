!> \brief Select one (or two in charge transfer) molecule and displace
!> the charge magnitude of charge sites on selected molecule(s)
!> according to charge nutrality and preferential strategy.
!> \author rewritten by Bin Chen at 6-25-99.
subroutine flucq (ichoice,boxi)
  use const_phys,only:qqfact
  use util_math,only:erfunc
  use util_random,only:random
  use sim_system
  use sim_cell,only:mimage
  use energy_kspace,only:recip,calp
  use energy_pairwise,only:energy
  implicit none

  logical::linterqt,ovrlap
  integer::i,ibox,iunit,j,imolty,icbu,mainunit,ichoice,qunit,boxi,flagon
  real::dchain,vn(nEnergy),vo(nEnergy),deltv,deltvb,dispbig,displit
  real::qion(numax)
  real::vrecipn,vrecipo
  integer::maini,mainj,jchain
  real::qionj(numax),vflucqjo,vflucqjn,corr,rij,rxuij,ryuij,rzuij,vdum
  real::qoldj2,vnewi,velectni,vinterni,voldi,velectoi,vinteroi

! --------------------------------------------------------------------
#ifdef __DEBUG__
  write(io_output,*) 'start FLUCQ in ',myid
#endif

! select a chain at random ***
      dchain  = random(-1)
      do icbu = 1,nmolty
         if ( dchain .lt. pmfqmt(icbu) ) then
            imolty = icbu
            dchain = 2.0E0_dp
         end if
      end do

      if ( lanes ) then

         if ( ichoice .eq. -1 ) then
! equal probability charge moves in two boxes (for swap)
! preferential charge moves according to favor
            dchain = dble(temtyp(imolty))
 77         i = int( dchain*random(-1) ) + 1
            i = parall(imolty,i)
            if ( random(-1) .gt. favor(i)) goto 77
            ibox = nboxi(i)
         else if ( ichoice .eq. -2 ) then
! equal probability charge moves in two boxes (for swap)
! preferential charge moves according to favor2
            dchain = dble(temtyp(imolty))
 66         i = int( dchain*random(-1) ) + 1
            i = parall(imolty,i)
            if ( random(-1) .gt. favor2(i)) goto 66
            ibox = nboxi(i)
         else if ( ichoice .eq. 0 ) then
! equal probability charge moves in separate box (for trans,rot)
 88         dchain = dble(temtyp(imolty))
            i = int( dchain*random(-1) ) + 1
            i = parall(imolty,i)
            ibox = nboxi(i)
            if ( ibox .ne. boxi ) goto 88
         else if (ichoice .eq. 2) then
            dchain = dble(temtyp(imolty))
            i = int( dchain*random(-1) ) + 1
            i = parall(imolty,i)
            ibox = nboxi(i)
         end if
      else
         dchain = dble(temtyp(imolty))
         i = int( dchain*random(-1) ) + 1
         i = parall(imolty,i)
         ibox = nboxi(i)

      end if

! For charge transfer case, select a second molecule to perform the
! fluctuating-charge move. The second molecule should be in the same
! box as chain i.

      If ( lqtrans(imolty) ) then
         if ( lanes ) then

            if ( ichoice .eq. -1 ) then
! equal probability charge moves in two boxes (for swap)
! preferential charge moves according to favor
               dchain = dble(temtyp(imolty))
 70            jchain = int( dchain*random(-1) ) + 1
               jchain = parall(imolty,jchain)
               if ( random(-1) .gt. favor(jchain)) goto 70
               if ( nboxi(jchain) .ne. ibox ) goto 70
            else if ( ichoice .eq. -2 ) then
! equal probability charge moves in two boxes (for swap)
! preferential charge moves according to favor
               dchain = dble(temtyp(imolty))
 60            jchain = int( dchain*random(-1) ) + 1
               jchain = parall(imolty,jchain)
               if ( random(-1) .gt. favor2(jchain)) goto 60
               if ( nboxi(jchain) .ne. ibox ) goto 60
            else if ( ichoice .eq. 0 ) then
! equal probability charge moves in separate box (for trans,rot)
 80            dchain = dble(temtyp(imolty))
               jchain = int( dchain*random(-1) ) + 1
               jchain = parall(imolty,jchain)
               if ( nboxi(jchain) .ne. boxi ) goto 80
            else if (ichoice .eq. 2) then
 90            dchain = dble(temtyp(imolty))
               jchain = int( dchain*random(-1) ) + 1
               jchain = parall(imolty,jchain)
               if ( nboxi(jchain) .ne. ibox ) goto 90
            end if
         else
 100        dchain = dble(temtyp(imolty))
            jchain = int( dchain*random(-1) ) + 1
            jchain = parall(imolty,jchain)
            if ( nboxi(jchain) .ne. ibox ) goto 100
         end if
      end if

      if ( lqtrans(imolty) .and. (jchain .ne. i ) ) then
         linterqt = .true.
      else
         linterqt = .false.
      end if

! store number of units of i in iunit ***

      iunit = nunit(imolty)
! count the number of charge sites
      qunit = 0
      do j = 1, iunit
         if ( lqchg(ntype(imolty,j)) ) qunit = qunit + 1
      end do

! store the charges in qion
      do  j = 1, iunit
         qion(j)   = qqu(i,j)
      end do

! calculate the polariztion energy of i in the old configuration ***

      call charge(i, qion, vo(ivFlucq))

      If ( linterqt ) then
         do j = 1,iunit
            qionj(j) = qqu(jchain,j)
         end do

! calculate the polarization energy of jchain in the old configuration ***

         call charge(jchain, qionj, vflucqjo)

         vo(ivFlucq) = vo(ivFlucq) + vflucqjo
      end if

! Choose one of the units as the main charge transfer site

 30   mainunit = int( dble(iunit)*random(-1) ) + 1
! for unit which is not a charge site
      if ( .not. lqchg(ntype(imolty,mainunit)) ) goto 30
      bnflcq(imolty,ibox) = bnflcq(imolty,ibox) + 1.0E0_dp
      dispbig = ( 2.0E0_dp*random(-1) - 1.0E0_dp )*rmflcq(imolty,ibox)

      if ( linterqt ) then
! For charge transfer case, i molecule increases by dispbig and
! jchain molecule decreases by dispbig.
! correction for the reptition of the calculation of the
! coulombic real space term between maini and mainj
         maini = mainunit
 32      mainunit = int( dble(iunit)*random(-1) ) + 1
! for unit which is not a charge site
         if ( .not. lqchg(ntype(imolty,mainunit)) ) goto 32
         mainj = mainunit
         corr = (qion(maini)-qionj(mainj))*dispbig  - qion(maini)*qionj(mainj)
         qion(maini) = qion(maini) + dispbig
         qionj(mainj) = qionj(mainj) - dispbig
         corr = corr + qion(maini)*qionj(mainj)
         rxuij = rxu(i,maini)-rxu(jchain,mainj)
         ryuij = ryu(i,maini)-ryu(jchain,mainj)
         rzuij = rzu(i,maini)-rzu(jchain,mainj)
         if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
         rij = sqrt( rxuij*rxuij + ryuij*ryuij + rzuij*rzuij )
         corr = qqfact*erfunc(calp(ibox)*rij)*corr/rij

!> \bug problems here for charge transfer

      else

         displit = (-dispbig) / (dble(qunit)-1.0E0_dp)

! displace the charges on the molecule

         do j = 1, iunit
            if ( j .eq. mainunit ) then
               qion(j) = qion(j) + dispbig
            else if ( lqchg(ntype(imolty,j)) ) then
               qion(j) = qion(j) + displit
            end if
            rxuion(j,1) = rxu(i,j)
            ryuion(j,1) = ryu(i,j)
            rzuion(j,1) = rzu(i,j)
            qquion(j,1) = qqu(i,j)
            rxuion(j,2) = rxuion(j,1)
            ryuion(j,2) = ryuion(j,1)
            rzuion(j,2) = rzuion(j,1)
            qquion(j,2) = qion(j)
         end do
      end if
      moltion(1) = imolty
      moltion(2) = imolty

! calculate the polarization energy of i in the new configuration ***

      call charge(i, qion, vn(ivFlucq))

      if ( linterqt ) then

! calculate the polarization energy of jchain in the new configuration ***

         call charge(jchain, qionj, vflucqjn)
         vn(ivFlucq) = vn(ivFlucq) + vflucqjn

! calculate the energy of i in the new configuration ***
         flagon = 2
         rxuion(maini,flagon) = rxu(i,maini)
         ryuion(maini,flagon) = ryu(i,maini)
         rzuion(maini,flagon) = rzu(i,maini)
         qquion(maini,flagon) = qion(maini)
         qoldj2 = qqu(jchain,mainj)
         qqu(jchain,mainj) = qionj(mainj)

         call energy(i,imolty,vn,flagon,ibox,maini,maini,.false.,ovrlap,.false.,.false.,.false.,.false.)
         vnewi = vn(ivTot) + vn(ivFlucq)
         velectni = vn(ivElect) + vn(ivEwald)
         vinterni = vn(ivInterLJ)

! calculate the energy of i in the old configuration ***
         flagon = 1
         rxuion(maini,flagon) = rxu(i,maini)
         ryuion(maini,flagon) = ryu(i,maini)
         rzuion(maini,flagon) = rzu(i,maini)
         qquion(maini,flagon) = qqu(i,maini)

         call energy(i,imolty,vo,flagon,ibox,maini,maini,.false.,ovrlap,.false.,.false.,.false.,.false.)
         voldi = vo(ivTot) + vo(ivFlucq)
         velectoi = vo(ivElect) + vo(ivEwald)
         vinteroi = vo(ivInterLJ)

! restore the old charges

         qqu(jchain,mainj) = qoldj2

! calculate the energy of jchain in the new configuration ***

         flagon = 2
         rxuion(mainj,flagon) = rxu(jchain,mainj)
         ryuion(mainj,flagon) = ryu(jchain,mainj)
         rzuion(mainj,flagon) = rzu(jchain,mainj)
         qquion(mainj,flagon) = qionj(mainj)
         call energy(jchain,imolty,vn,flagon,ibox,mainj,mainj,.false.,ovrlap,.false.,.false.,.false.,.false.)
         vn(ivTot) = vn(ivTot) + vnewi
         vn(ivElect) = vn(ivElect) + vn(ivEwald) + velectni
         vn(ivInterLJ) = vn(ivInterLJ) + vinterni

! calculate the energy of jchain in the old configuration ***

         flagon = 1
         rxuion(mainj,flagon) = rxu(jchain,mainj)
         ryuion(mainj,flagon) = ryu(jchain,mainj)
         rzuion(mainj,flagon) = rzu(jchain,mainj)
         qquion(mainj,flagon) = qqu(jchain,mainj)

         call energy(jchain,imolty,vo,flagon,ibox,mainj,mainj,.false.,ovrlap,.false.,.false.,.false.,.false.)

         vo(ivTot) = vo(ivTot) + voldi
         vo(ivElect) = vo(ivElect) + vo(ivEwald) + velectoi
         vo(ivInterLJ) = vo(ivInterLJ) + vinteroi

! prepare the rxuion etc for ewald sum and dielectric constant
! and for energy calculation

         if ( linterqt ) then
            rxuion(1,1) = rxu(i,maini)
            ryuion(1,1) = ryu(i,maini)
            rzuion(1,1) = rzu(i,maini)
            qquion(1,1) = qqu(i,maini)
            rxuion(2,1) = rxu(jchain,mainj)
            ryuion(2,1) = ryu(jchain,mainj)
            rzuion(2,1) = rzu(jchain,mainj)
            qquion(2,1) = qqu(jchain,mainj)
            rxuion(1,2) = rxu(i,maini)
            ryuion(1,2) = ryu(i,maini)
            rzuion(1,2) = rzu(i,maini)
            qquion(1,2) = qion(maini)
            rxuion(2,2) = rxu(jchain,mainj)
            ryuion(2,2) = ryu(jchain,mainj)
            rzuion(2,2) = rzu(jchain,mainj)
            qquion(2,2) = qionj(mainj)

            do j = 3, iunit
               rxuion(j,1) = rxuion(j,2)
               ryuion(j,1) = ryuion(j,2)
               rzuion(j,1) = rzuion(j,2)
               qquion(j,1) = qquion(j,2)
            end do
         end if

      else

! calculate the energy of i in the new configuration ***
         flagon = 2
         call energy(i,imolty,vn,flagon,ibox,1,iunit,.false.,ovrlap,.false.,.false.,.false.,.false.)
         if (ovrlap) return
         vn(ivTot) = vn(ivTot) + vn(ivFlucq)
         vn(ivElect) = vn(ivElect) + vn(ivEwald)

! calculate the energy of i in the old configuration ***
         flagon = 1
         call energy(i,imolty,vo,flagon,ibox,1,iunit,.false.,ovrlap,.false.,.false.,.false.,.false.)
         vo(ivTot) = vo(ivTot) + vo(ivFlucq)
         vo(ivElect) = vo(ivElect) + vo(ivEwald)

      end if

! Begin Ewald-sum correction
      if ( lewald ) then
         call recip(ibox,vrecipn,vrecipo,1)
         vn(ivTot) = vn(ivTot) + vrecipn
         vo(ivTot) = vo(ivTot) + vrecipo
         vn(ivElect) = vn(ivElect) + vrecipn
         vo(ivElect) = vo(ivElect) + vrecipo
      end if


! check for acceptance ***

      deltv  = vn(ivTot) - vo(ivTot)
! use the thermostat temperature instead of real temp
      deltvb = fqbeta * deltv

! if ( deltv .lt. -100.0E0_dp) then
! write(io_output,*) i,favor(i),deltv
! end if
      if ( deltvb .gt. (2.3E0_dp*softcut) ) return

      if ( deltv .le. 0.0E0_dp ) then
! accept move
      else if ( exp(-deltvb) .gt. random(-1) ) then
! accept move
      else
         return
      end if

      vbox(ivTot,ibox)     = vbox(ivTot,ibox) + deltv
      vbox(ivElect,ibox)  = vbox(ivElect,ibox)  + (vn(ivElect) - vo(ivElect))
      vbox(ivFlucq,ibox)  = vbox(ivFlucq,ibox) + (vn(ivFlucq) - vo(ivFlucq))
      vbox(ivInterLJ,ibox) = vbox(ivInterLJ,ibox) + (vn(ivInterLJ) - vo(ivInterLJ))
! write(io_output,*) 'this move has been accepted!!!'
      do j = 1,iunit
         qqu(i,j) = qion(j)
         if ( linterqt )  qqu(jchain,j) = qionj(j)
      end do
! update the reciprocal-space sum
      if ( lewald ) then
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
          call dipole(ibox,1)
      end if

      bsflcq(imolty,ibox) = bsflcq(imolty,ibox) + 1.0E0_dp

#ifdef __DEBUG__
      write(io_output,*) 'end FLUCQ in ',myid
#endif
      return
    end subroutine flucq


