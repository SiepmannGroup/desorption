!> \brief Optimize the electronic configuration for trans, rot, and swap
!> (config, swatch in the future) moves and accept/reject the
!> combined move.
!> \author written on June 25/99 by Bin Chen.
subroutine anes(i,ibox,boxrem,mtype,laccept,deltv,vn,vo,vinsta,vremta,vnewflucq,voldflucq,lswapinter)
  use util_random,only:random
  use sim_system
  use sim_particle,only:ctrmas
  use energy_kspace,only:recip
  implicit none

  logical::laccept,lswapinter
  integer::i,ibox,boxrem,mtype,imolty,iunit,ichoiq,ip,ibox2,j
  real::deltv,vn(nEnergy),vo(nEnergy),vinsta,vremta,vnewflucq,voldflucq,vboxo(nbxmax)&
   ,vinterbo(nbxmax),vintrabo(nbxmax),vextbo(nbxmax),velectbo(nbxmax),vflucqbo(nbxmax),vtailbo(nbxmax),vvibbo(nbxmax)&
   ,vtgbo(nbxmax),vbendbo(nbxmax),rxuo(numax),ryuo(numax),rzuo(numax),xcmo,ycmo,zcmo,vdum,wratio,volins,volrem,deltvb,vnewt2,voldt2
  real::qquo(nmax,numax)

#ifdef __DEBUG__
      write(io_output,*) 'START the optimization of the charge configuration in ',myid
#endif

      imolty = moltyp(i)
      iunit = nunit(imolty)

! store the old energy, old coordinates and ewald sum
      do ibox2 = 1,nbox
         vboxo(ibox2) = vbox(ivTot,ibox2)
         vinterbo(ibox2) = vbox(ivInterLJ,ibox2)
         vintrabo(ibox2) = vbox(ivIntraLJ,ibox2)
         vextbo(ibox2) = vbox(ivExt,ibox2)
         vflucqbo(ibox2) = vbox(ivFlucq,ibox2)
         velectbo(ibox2) = vbox(ivElect,ibox2)

         if ( mtype .eq. 3 ) then
            vtailbo(ibox2) = vbox(ivTail,ibox2)
            vvibbo(ibox2) = vbox(ivStretching,ibox2)
            vtgbo(ibox2) = vbox(ivTorsion,ibox2)
            vbendbo(ibox2) = vbox(ivBending,ibox2)
         end if

      end do
      do j = 1,iunit
         rxuo(j) = rxu(i,j)
         ryuo(j) = ryu(i,j)
         rzuo(j) = rzu(i,j)
      end do
      xcmo = xcm(i)
      ycmo = ycm(i)
      zcmo = zcm(i)
! store the old charges
      do ip = 1,nchain
         do j = 1,nunit(moltyp(ip))
            qquo(ip,j) = qqu(ip,j)
         end do
      end do
      if (lewald) then
! store the reciprocal-space sum
         do ibox2 = 1,nbox
            call recip(ibox2,vdum,vdum,3)
         end do
         if ( ldielect ) then
            call dipole(ibox,2)
         end if
      end if

! on the new coordinates, continue to use the fluctuating charge
! algorithm to optimize the charge configurations, update the
! energy, coordinates and the ewald sum

      if ( mtype .eq. 3 ) then
! for swap move
         vbox(ivTot,ibox)     = vbox(ivTot,ibox) + vnew(ivTot)
         vbox(ivInterLJ,ibox)  = vbox(ivInterLJ,ibox) + vnew(ivInterLJ)
         vbox(ivIntraLJ,ibox)  = vbox(ivIntraLJ,ibox) + vnew(ivIntraLJ)
         vbox(ivExt,ibox)    = vbox(ivExt,ibox)   + vnew(ivExt)
         vbox(ivElect,ibox)   = vbox(ivElect,ibox)  + vnew(ivElect)+vnew(ivEwald)
         vbox(ivTail,ibox)   = vbox(ivTail,ibox)   + vinsta
         vbox(ivStretching,ibox)    =  vbox(ivStretching,ibox)   + vnew(ivStretching)
         vbox(ivTorsion,ibox)     = vbox(ivTorsion,ibox)     + vnew(ivTorsion)
         vbox(ivBending,ibox)   = vbox(ivBending,ibox)   + vnew(ivBending)
         vbox(ivFlucq,ibox)  = vbox(ivFlucq,ibox)  + vnewflucq

         vbox(ivTot,boxrem)     = vbox(ivTot,boxrem)     - vold(ivTot)
         vbox(ivInterLJ,boxrem)  = vbox(ivInterLJ,boxrem)  - vold(ivInterLJ)
         vbox(ivTail,boxrem)   = vbox(ivTail,boxrem)   - vremta
         vbox(ivIntraLJ,boxrem)  = vbox(ivIntraLJ,boxrem)  - vold(ivIntraLJ)
         vbox(ivStretching,boxrem)    = vbox(ivStretching,boxrem)    - vold(ivStretching)
         vbox(ivTorsion,boxrem)     = vbox(ivTorsion,boxrem)     - vold(ivTorsion)
         vbox(ivExt,boxrem)    = vbox(ivExt,boxrem)    - vold(ivExt)
         vbox(ivBending,boxrem)   = vbox(ivBending,boxrem)   - vold(ivBending)
         vbox(ivElect,boxrem)  = vbox(ivElect,boxrem)  -  (vold(ivElect)+vold(ivEwald))
         vbox(ivFlucq,boxrem)  = vbox(ivFlucq,boxrem)  - voldflucq
      else

         vbox(ivTot,ibox)     = vbox(ivTot,ibox) + deltv
         vbox(ivInterLJ,ibox)  = vbox(ivInterLJ,ibox) + (vn(ivInterLJ) - vo(ivInterLJ))
         vbox(ivIntraLJ,ibox)  = vbox(ivIntraLJ,ibox) + (vn(ivIntraLJ) - vo(ivIntraLJ))
         vbox(ivExt,ibox)    = vbox(ivExt,ibox)   + (vn(ivExt)   - vo(ivExt))
         vbox(ivElect,ibox)   = vbox(ivElect,ibox)  + (vn(ivElect) - vo(ivElect))
      end if

      do j = 1,iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      end do
! update chain center of mass
      call ctrmas(.false.,ibox,i,mtype)

      if (lewald) then
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
         if (mtype .eq. 3) call recip(boxrem,vdum,vdum,2)
         if ( ldielect ) then
! update the dipole term
            call dipole(ibox,1)
         end if
      end if

! begin to optimize the charge configuration

      if ( mtype .eq. 3 ) then
! swap move
! Bin's recommendations are 1000 total swap moves, 500 biases to be near the swapped molecule,
! 500 without bias
         do ichoiq = 1,nswapq
            call flucq(-1,0) ! call with -1 selects molecule according to favor(i) biasing from either box
         end do

         do ichoiq = 1,500
            call flucq(2,0) ! call with 2 selects molecule without bias from either box
         end do
         deltv = vbox(ivTot,ibox) - vboxo(ibox)
         weight = exp(-deltv*beta)
         deltv = vboxo(boxrem) - vbox(ivTot,boxrem)
         weiold = exp(-deltv*beta)
         volins=boxlx(ibox)*boxly(ibox)*boxlz(ibox)
         volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)

         if ( lswapinter ) then
            if (lgibbs) then
! Note: acceptance based on only molecules of type imolty
               wratio = ( weight / weiold ) * ( volins * dble( ncmt(boxrem,imolty)+1 ) /  ( volrem * dble( ncmt(ibox,imolty) ) ) )
            else if (lgrand) then
               if (ibox.eq.1) then
! molecule added to box 1
                  wratio = (weight /  weiold ) *  volins * B(imolty) / (ncmt(ibox,imolty))
               else
! molecule removed from box 1
                  wratio = (weight /  weiold ) *  (ncmt(boxrem,imolty)+1)/ (B(imolty)*volrem)
               end if
            end if
         else
               wratio = weight / weiold
         end if
         if ( wratio .gt. random(-1) ) then
            laccept = .true.
         else
            laccept = .false.
         end if

      else
! trans/rot moves
         do ichoiq = 1,nchoiq(ibox)
            call flucq(0,ibox)  !  call with 0, ibox selects random molecule (no biases) in the specified box
         end do
         vnewt2 = 0.0E0_dp
         voldt2 = 0.0E0_dp
         do ibox2 = 1, nbox
            vnewt2 = vnewt2 + vbox(ivTot,ibox2)
            voldt2 = voldt2 + vboxo(ibox2)
         end do
         deltv = vnewt2 - voldt2

         deltvb = beta * deltv
         if ( deltvb .gt. (2.3E0_dp*softcut) ) then
            laccept = .false.
         else if ( deltv .le. 0.0E0_dp ) then
            laccept = .true.
         else if ( exp(-deltvb) .gt. random(-1) ) then
            laccept = .true.
         else
            laccept = .false.
         end if

      end if

      if ( laccept ) then

! combined move can be accepted now !!!

      else
! restore the old energy and old coordinates and ewald sum
         do ibox2 = 1,nbox
            vbox(ivTot,ibox2) = vboxo(ibox2)
            vbox(ivInterLJ,ibox2) = vinterbo(ibox2)
            vbox(ivIntraLJ,ibox2) = vintrabo(ibox2)
            vbox(ivExt,ibox2) = vextbo(ibox2)
            vbox(ivElect,ibox2) = velectbo(ibox2)
            vbox(ivFlucq,ibox2) = vflucqbo(ibox2)
            if ( mtype .eq. 3 ) then
               vbox(ivTail,ibox2) = vtailbo(ibox2)
               vbox(ivStretching,ibox2) = vvibbo(ibox2)
               vbox(ivTorsion,ibox2) = vtgbo(ibox2)
               vbox(ivBending,ibox2) = vbendbo(ibox2)
            end if
         end do
         do j = 1,iunit
            rxu(i,j) = rxuo(j)
            ryu(i,j) = ryuo(j)
            rzu(i,j) = rzuo(j)
         end do
         do ip = 1,nchain
            do j = 1,nunit(moltyp(ip))
               qqu(ip,j) = qquo(ip,j)
            end do
         end do
         xcm(i) = xcmo
         ycm(i) = ycmo
         zcm(i) = zcmo
         if (lewald) then
! restore the reciprocal-space sum
            do ibox2 = 1,nbox
               call recip(ibox2,vdum,vdum,4)
            end do
            if ( ldielect ) then
! restore old dipole moment
               call dipole(ibox,3)
            end if
         end if
      end if

#ifdef __DEBUG__
      write(io_output,*) 'END the optimization of the charge configuration in ',myid
#endif

      return
    end subroutine anes





