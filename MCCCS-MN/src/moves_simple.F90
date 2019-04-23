MODULE moves_simple
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use energy_kspace,only:recip
  use energy_pairwise,only:energy
  implicit none
  private
  save
  public::translation,rotation,Atom_translation,init_moves_simple,update_translation_rotation_max_displacement&
   ,averageMaximumDisplacement,output_translation_rotation_stats,read_checkpoint_simple,write_checkpoint_simple

  real,allocatable,public::acntrax(:,:),acntray(:,:),acntraz(:,:),acnrotx(:,:),acnroty(:,:),acnrotz(:,:),acstrax(:,:)&
   ,acstray(:,:),acstraz(:,:),acsrotx(:,:),acsroty(:,:),acsrotz(:,:),bntrax(:,:),bntray(:,:),bntraz(:,:),bstrax(:,:)&
   ,bstray(:,:),bstraz(:,:),bnrotx(:,:),bnroty(:,:),bnrotz(:,:),bsrotx(:,:),bsroty(:,:),bsrotz(:,:)
  real::Abntrax,Abntray,Abntraz,Abstrax,Abstray,Abstraz

contains
!> \brief Makes a translational movement in x,y,or z-direction.
!>
!> The maximum displacement is controlled by \b rmtrax(yz) and the
!> number of successful trial moves is stored in \b bstrax(yz).
!> The attempts are stored in \b bntrax(yz)
  subroutine translation()
    use sim_particle,only:update_neighbor_list_molecule,ctrmas,update_coord_in_tree
    use sim_cell,only:update_linked_cell

      logical::lx,ly,lz,ovrlap

      integer::i,ibox,flagon,iunit,j,imolty,icbu
      real::rx,ry,rz,dchain,ddx,ddy,ddz,vn(nEnergy),vo(nEnergy),deltv,deltvb,rchain,vdum,vrecipo,vrecipn

      logical::laccept
! --------------------------------------------------------------------
#ifdef __DEBUG__
      write(io_output,*) 'start TRANSLATION in ',myid
#endif
      lx=.false.
      ly=.false.
      lz=.false.
      rchain = 3.0E0_dp * random(-1)
      if ( rchain .lt. 1.0E0_dp ) then
         lx=.true.
      else if ( rchain .lt. 2.0E0_dp ) then
         ly=.true.
      else
         lz=.true.
      end if

      ovrlap = .false.
! select a chain at random ***
      rchain  = random(-1)
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            exit
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

      if (temtyp(imolty).eq.0) return

      dchain = real(temtyp(imolty),dp)
      i = int( dchain*random(-1) + 1 )
      i = parall(imolty,i)
      ibox = nboxi(i)

! store number of units of i in iunit ***

      iunit = nunit(imolty)

      do  j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      end do
      moltion(1) = imolty

      !write(io_output,*) 'particle ',i,' is being translated, ',lx,ly,lz,'imolty is:',imolty,' and the box is:',ibox
! move i ***
      if (lx) then
         rx =  ( 2.0*random(-1) - 1.0E0_dp ) * rmtrax(imolty,ibox)
         bntrax(imolty,ibox) = bntrax(imolty,ibox) + 1.0E0_dp
      else
         rx=0
      end if
      if (ly) then
         ry =  ( 2.0*random(-1) - 1.0E0_dp ) * rmtray(imolty,ibox)
         bntray(imolty,ibox) = bntray(imolty,ibox) + 1.0E0_dp
      else
         ry=0
      end if
      if (lz) then
         rz =  ( 2.0*random(-1) - 1.0E0_dp ) * rmtraz(imolty,ibox)
         bntraz(imolty,ibox) = bntraz(imolty,ibox) + 1.0E0_dp
      else
         rz=0
      end if
      do j = 1, iunit
         ddx=0
         ddy=0
         ddz=0
         if ( .not. (lsolid(ibox) .and. .not. lrect(ibox)) ) then
            if (lx) then
               if ( lfold .and. lpbcx ) then
                  if ((xcm(i) + rx) .lt. 0.0E0_dp) ddx = boxlx(ibox)
                  if ((xcm(i) + rx) .gt. boxlx(ibox)) ddx = -boxlx(ibox)
               end if
            end if
            if (ly) then
               if ( lfold .and. lpbcy ) then
                  if ((ycm(i) + ry).lt.0.0E0_dp)  ddy=boxly(ibox)
                  if ((ycm(i) + ry).gt.boxly(ibox))  ddy=-boxly(ibox)
               end if
            end if
            if (lz) then
               if ( lfold .and. lpbcz ) then
                  if ((zcm(i) + rz).lt.0.0E0_dp)  ddz=boxlz(ibox)
                  if ((zcm(i) + rz).gt.boxlz(ibox))  ddz=-boxlz(ibox)
               end if
            end if
         end if
         rxuion(j,2) = rxuion(j,1) + rx + ddx
         ryuion(j,2) = ryuion(j,1) + ry + ddy
         rzuion(j,2) = rzuion(j,1) + rz + ddz
         qquion(j,2) = qquion(j,1)
      end do
      moltion(2) = imolty

! calculate the energy of i in the new configuration ***
      flagon = 2
      call energy(i,imolty,vn,flagon,ibox,1,iunit,.false.,ovrlap,.false. ,.false.,.false.,.false.)
      if (ovrlap) return

! calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vo,flagon,ibox,1,iunit,.false.,ovrlap,.false. ,.false.,.false.,.false.)

      if (ovrlap) then
         call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf of TRANSLATION',myid+1)
      end if

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
         call recip(ibox,vrecipn,vrecipo,1)
         vn(ivElect) = vn(ivElect) + vrecipn
         vo(ivElect) = vo(ivElect) + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
            vrecipn    =   (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipn
            vrecipo    =   (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipo
         else if (lstageb) then
            vrecipn =  etais*vrecipn
            vrecipo =  etais*vrecipo
         else if (lstagec) then
            vrecipn =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipn
            vrecipo =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipo
         end if
         vn(ivTot) = vn(ivTot) + vrecipn
         vo(ivTot) = vo(ivTot) + vrecipo
      end if
! check for acceptance ***

      deltv  = vn(ivTot) - vo(ivTot)
      deltvb = beta * deltv


! For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vn,vo,vdum,vdum,vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) bstrax(imolty,ibox) = bstrax(imolty,ibox) + 1.0E0_dp
            if (ly) bstray(imolty,ibox) = bstray(imolty,ibox) + 1.0E0_dp
            if (lz) bstraz(imolty,ibox) = bstraz(imolty,ibox) + 1.0E0_dp
         end if
         return
      end if

      if ( deltvb .gt. (2.3E0_dp*softcut) ) return

      if ( deltv .le. 0.0E0_dp ) then
! accept move
      else if ( exp(-deltvb) .gt. random(-1) ) then
! accept move
      else
! move rejected
         return
      end if

! write(io_output,*) 'TRANSLATION accepted i',i

      vbox(ivTot,ibox)     = vbox(ivTot,ibox) + deltv
      vbox(ivInterLJ,ibox)  = vbox(ivInterLJ,ibox) + (vn(ivInterLJ) - vo(ivInterLJ))
      vbox(ivIntraLJ,ibox)  = vbox(ivIntraLJ,ibox) + (vn(ivIntraLJ) - vo(ivIntraLJ))
      vbox(ivExt,ibox)    = vbox(ivExt,ibox)   + (vn(ivExt) - vo(ivExt))
      vbox(ivElect,ibox)   = vbox(ivElect,ibox)  + (vn(ivElect) - vo(ivElect))
      vbox(ivIpswb,ibox) = vbox(ivIpswb,ibox) + (vipswn-vipswo)
      vbox(ivWellIpswb,ibox) = vbox(ivWellIpswb,ibox) + (vwellipswn-vwellipswo)
      vipsw = vbox(ivIpswb,ibox)
      vwellipsw = vbox(ivWellIpswb,ibox)
      vbox(iv3body,ibox) = vbox(iv3body,ibox) + (vn(iv3body)-vo(iv3body))

      ! Update beads coordinates in kdtree
      if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
          do j = 1, iunit
              rxu_update(j) = rxuion(j,2)
              ryu_update(j) = ryuion(j,2)
              rzu_update(j) = rzuion(j,2)
          end do

          call update_coord_in_tree(i, iunit, ibox, ibox, .true., .false.)
      end if


      ! Update coordinates in r*u arrays
      do j = 1,iunit
         if (lx) rxu(i,j) = rxuion(j,2)
         if (ly) ryu(i,j) = ryuion(j,2)
         if (lz) rzu(i,j) = rzuion(j,2)
      end do

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
! update the dipole term
         call dipole(ibox,1)
      end if

! update chain center of mass

      call ctrmas(.false.,ibox,i,1)

      if (lx) bstrax(imolty,ibox) = bstrax(imolty,ibox) + 1.0E0_dp
      if (ly) bstray(imolty,ibox) = bstray(imolty,ibox) + 1.0E0_dp
      if (lz) bstraz(imolty,ibox) = bstraz(imolty,ibox) + 1.0E0_dp

      if ( licell .and. (ibox.eq.boxlink)) then
! update linkcell list
         call update_linked_cell(i)
      end if

      if (lneigh.or.lneighbor.or.lgaro) then
         call update_neighbor_list_molecule(i)
      end if

#ifdef __DEBUG__
      write(io_output,*) 'end TRANSLATION in ',myid,i
#endif

      return
    end subroutine translation

!> \brief Makes a rotational movement around "x" space-fixed axis.
!>
!> The maximum displacement is controlled by \b rmrotx and the
!> number of successful rotation is given by \b bsrotx. \n
!> Rotation chooses one of the three space-fixed axes at random
!> and rotates the molecule around this axis by \b dgamma radians.
!> The maximum angular displacement is \b dgamax.
    subroutine rotation()
      use sim_particle,only:update_neighbor_list_molecule,ctrmas,update_coord_in_tree

      logical::lx,ly,lz,ovrlap
      integer::i,ibox,flagon,iunit,j,imolty,iuroty,icbu
      real::rx,ry,rz,dchain,rchain,vn(nEnergy),vo(nEnergy),dgamma,rxorig,ryorig,rzorig,rxnew2,rynew2,rznew2,deltv,deltvb,vdum
! further variable definitions
      real::cosdg,sindg,rmrot
      real::vrecipn,vrecipo
      logical::laccept

! --------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start ROTATION in ',myid
#endif

      lx=.false.
      ly=.false.
      lz=.false.
      rchain = 3.0E0_dp * random(-1)
      if ( rchain .lt. 1.0E0_dp ) then
         lx=.true.
      else if ( rchain .lt. 2.0E0_dp ) then
         ly=.true.
      else
         lz=.true.
      end if

      ovrlap = .false.
! select a chain type at random ***
      rchain  = random(-1)
      do icbu = 1,nmolty
         if ( rchain .lt. pmromt(icbu) ) then
            imolty = icbu
            exit
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

      if (temtyp(imolty).eq.0) return

      dchain = real(temtyp(imolty),dp)
      i = int( dchain*random(-1) + 1 )
      i = parall(imolty,i)
      ibox = nboxi(i)

      iuroty = iurot(imolty)
!kea 6/4/09 -- for multiple rotation centers
      if(iuroty.lt.0) then
         if(nrotbd(imolty).gt.1) then
            rchain = random(-1)
            do icbu = 1,nrotbd(imolty)
               if( rchain .lt. pmrotbd(icbu,imolty)) then
                  iuroty = irotbd(icbu,imolty)
               end if
            end do
         else
            iuroty = irotbd(1,imolty)
         end if
      end if

! store number of units of i in iunit ***
      iunit = nunit(imolty)

! store current positions in old-new array #1# ***
      do j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      end do
      moltion(1) = imolty

      !write(io_output,*) 'particle ',i,' is being rotated, ',lx,ly,lz,'imolty is:',imolty,' and the box is:',ibox
! choose a random angular displacement ***

      if (lx) then
         rmrot = rmrotx(imolty,ibox)
         bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0E0_dp
      else if (ly) then
         rmrot = rmroty(imolty,ibox)
         bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0E0_dp
      else if (lz) then
         rmrot = rmrotz(imolty,ibox)
         bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0E0_dp
      end if
      dgamma = ( 2.0E0_dp*random(-1) - 1.0E0_dp ) * rmrot

! set up the rotation marix ***

      cosdg = cos( dgamma )
      sindg = sin( dgamma )

! Determine the rotation coordinates ***
      if (iuroty .eq. 0) then
! Use the center of mass for rotation
         rxorig = xcm(i)
         ryorig = ycm(i)
         rzorig = zcm(i)
      else
! Use iurot for rotation
         rxorig = rxuion(iuroty,1)
         ryorig = ryuion(iuroty,1)
         rzorig = rzuion(iuroty,1)
      end if

! write(io_output,*) 'before rotating'
! write(io_output,*) xcm(i),ycm(i),zcm(i)
! write(io_output,*) rxu(i,1),ryu(i,1),rzu(i,1)


      if (lx) then

! ROTATE UNITS OF I AROUND X-AXIS ***

         do  j = 1, iunit
            ry = ryuion(j,1) - ryorig
            rz = rzuion(j,1) - rzorig
            rynew2 = cosdg * ry + sindg * rz
            rznew2 = cosdg * rz - sindg * ry

            rxuion(j,2) = rxuion(j,1)
            ryuion(j,2) = ryorig + rynew2
            rzuion(j,2) = rzorig + rznew2
            qquion(j,2) = qquion(j,1)
         end do
      end if
      if (ly) then

! ROTATE UNITS OF I AROUND y-AXIS ***

         do  j = 1, iunit
            rx = rxuion(j,1) - rxorig
            rz = rzuion(j,1) - rzorig
            rxnew2 = cosdg * rx - sindg * rz
            rznew2 = cosdg * rz + sindg * rx

            rxuion(j,2) = rxorig + rxnew2
            ryuion(j,2) = ryuion(j,1)
            rzuion(j,2) = rzorig + rznew2
            qquion(j,2) = qquion(j,1)
         end do
      end if
      if (lz) then

! ROTATE UNITS OF I AROUND z-AXIS ***

         do  j = 1, iunit
            rx = rxuion(j,1) - rxorig
            ry = ryuion(j,1) - ryorig
            rxnew2 = cosdg * rx + sindg * ry
            rynew2 = cosdg * ry - sindg * rx

            rxuion(j,2) = rxorig + rxnew2
            ryuion(j,2) = ryorig + rynew2
            rzuion(j,2) = rzuion(j,1)
            qquion(j,2) = qquion(j,1)
         end do
      end if
      moltion(2) = imolty

! calculate the energy of i in the new configuration ***

      flagon = 2
      call energy(i,imolty,vn,flagon,ibox,1,iunit,.false.,ovrlap,.false.,.false.,.false.,.false.)
      if (ovrlap) return

! calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vo,flagon,ibox,1,iunit,.false.,ovrlap,.false.,.false.,.false.,.false.)

      if (ovrlap) call err_exit(__FILE__,__LINE__,'disaster- overlap for old conf in ROTATION',myid+1)
      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
         call recip(ibox,vrecipn,vrecipo,1)
         vn(ivElect) = vn(ivElect) + vrecipn
         vo(ivElect) = vo(ivElect) + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
           vrecipn  =  (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipn
           vrecipo  =  (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipo
         else if (lstageb) then
           vrecipn  =  etais*vrecipn
           vrecipo  =  etais*vrecipo
         else if (lstagec) then
           vrecipn  =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipn
           vrecipo  =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipo
         end if
         vn(ivTot) = vn(ivTot) + vrecipn
         vo(ivTot) = vo(ivTot) + vrecipo
      end if

! check for acceptance ***

      deltv  = vn(ivTot) - vo(ivTot)
      deltvb = beta * deltv

! For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,2,laccept,deltv,vn,vo,vdum,vdum,vdum,vdum,.false.)
         if (laccept) then
            if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0E0_dp
            if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0E0_dp
            if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0E0_dp
         end if
         return
      end if

      if ( deltvb .gt. (2.3E0_dp*softcut) ) return

      if ( deltv .le. 0.0E0_dp ) then
! move accepted
      else if ( exp(-deltvb) .gt. random(-1) ) then
! move accepted
      else
! move rejected
         return
      end if

! write(io_output,*) 'ROTATION accepted',i
      vbox(ivTot,ibox) = vbox(ivTot,ibox) + deltv
      vbox(ivInterLJ,ibox)  = vbox(ivInterLJ,ibox) + (vn(ivInterLJ) -  vo(ivInterLJ))
      vbox(ivIntraLJ,ibox)  = vbox(ivIntraLJ,ibox) + (vn(ivIntraLJ) - vo(ivIntraLJ))
      vbox(ivExt,ibox)    = vbox(ivExt,ibox)   + (vn(ivExt)-vo(ivExt))
      vbox(ivElect,ibox)  = vbox(ivElect,ibox) + (vn(ivElect)-vo(ivElect))
      vbox(ivIpswb,ibox) = vbox(ivIpswb,ibox) + (vipswn-vipswo)
      vbox(ivWellIpswb,ibox) = vbox(ivWellIpswb,ibox) + (vwellipswn-vwellipswo)
      vipsw = vbox(ivIpswb,ibox)
      vwellipsw = vbox(ivWellIpswb,ibox)

      ! Update coordinates in bead-based kdtree
      if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
          do j = 1, iunit
              rxu_update(j) = rxuion(j,2)
              ryu_update(j) = ryuion(j,2)
              rzu_update(j) = rzuion(j,2)
          end do

          call update_coord_in_tree(i, iunit, ibox, ibox, .true., .false.)
      end if

      ! Update coordinates in r*u arrays
      do j = 1, iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      end do

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
         call dipole(ibox,1)
      end if


! update center of mass
      call ctrmas(.false.,ibox,i,2)

      if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0E0_dp
      if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0E0_dp
      if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0E0_dp

      if (lneigh.or.lneighbor) then
         call update_neighbor_list_molecule(i)
      end if

#ifdef __DEBUG__
      write(io_output,*) 'end ROTATION in ',myid,i
#endif

      return
    end subroutine rotation

!> \brief Makes a translational movement in x,y,or z-direction.
!>
!> The maximum displacement is controlled by \b rAtrax(yz) and the
!> number of successful trial moves is stored in \b Abstrax(yz).
!> The attempts are stored in \b Abntrax(yz)

! -- RKL 4/14
!
!> \brief  allow atom translations only on specified atoms
!
!> Previously the user could only choose to do atom translations on every atom in the system or none at all. Now, users
!> can specify the specific atoms they want to do translations on.

    subroutine Atom_translation()
      use sim_particle,only:update_neighbor_list_molecule,ctrmas
      use sim_cell,only:update_linked_cell
      use energy_kspace,only:recip_atom
      use energy_intramolecular,only:U_bonded

      logical::lx,ly,lz,ovrlap
      integer::i,ibox,flagon,iunit,j,imolty,icbu
      integer::pick_unit, pick_chain
      real::rx,ry,rz,dchain,vn(nEnergy),vo(nEnergy),deltv,deltvb,rchain,vdum,vrecipo,vrecipn
      logical::laccept
      integer::atom_sel, molt_sel, mole_sel, rand_idx

! --------------------------------------------------------------------
#ifdef __DEBUG__
      write(io_output,*) 'start ATOM_TRANSLATION in ',myid
#endif

      lx=.false.
      ly=.false.
      lz=.false.
      rchain = 3.0E0_dp * random(-1)
      if ( rchain .lt. 1.0E0_dp ) then
         lx=.true.
      else if ( rchain .lt. 2.0E0_dp ) then
         ly=.true.
      else
         lz=.true.
      end if

      ovrlap = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Select atom type and get the corresponding molecule type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       rand_idx = random(-1)*(natomtrans_atoms)+1
       atom_sel = atomtrans_atomlst (rand_idx)
       imolty = atomtrans_moleclst(rand_idx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now pick the specific molecule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         dchain = real(temtyp(imolty),dp)
         pick_chain = int( dchain*random(-1) + 1 )
         pick_chain = parall(imolty,pick_chain)
         ibox = nboxi(pick_chain)

! Store the atom of interest, the molecule of interest, and the
! number of atoms in the molecule

      pick_unit = atom_sel
      i         = pick_chain
      iunit     = nunit(imolty)

! Used for old enegy evaluation (non-bonded)

      do j = 1,iunit
        rxuion(j,1) = rxu(i,j)
        ryuion(j,1) = ryu(i,j)
        rzuion(j,1) = rzu(i,j)
        qquion(j,1) = qqu(i,j)
      end do

      moltion(1) = imolty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Select random displacement, increment attempts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (lx) then
         rx =  ( 2.0*random(-1) - 1.0E0_dp ) * Armtrax
         Abntrax = Abntrax + 1.0E0_dp
      else
         rx=0
      end if
      if (ly) then
         ry =  ( 2.0*random(-1) - 1.0E0_dp ) * Armtray
         Abntray = Abntray + 1.0E0_dp
      else
         ry=0
      end if
      if (lz) then
         rz =  ( 2.0*random(-1) - 1.0E0_dp ) * Armtraz
         Abntraz = Abntraz + 1.0E0_dp
      else
         rz=0
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the energy of i in the old configuration
! Nonbonded is calculated in energy and is based on r*uion coords
! Bonded is calculated in U_bonded and is based on r*u coords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      flagon = 1     ! 1 refers to the coordinates for the OLD configuration
      call energy(i,imolty,vo,flagon,ibox,pick_unit,pick_unit,.true.,ovrlap,.false.,.false.,.false.,.true.)
      if (ovrlap) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf of ATOM_TRANSLATION',myid+1)

      call U_bonded(i,imolty,vo(ivStretching),vo(ivBending),vo(ivTorsion))


! Update the CURRENT coordinates with our attempted translation

      do j = 1,iunit
        if (j .eq. pick_unit) then
           rxuion(j,2) = rxuion(j,1) + rx
           ryuion(j,2) = ryuion(j,1) + ry
           rzuion(j,2) = rzuion(j,1) + rz
           qquion(j,2) = qquion(j,1)
        else
           rxuion(j,2) = rxuion(j,1)
           ryuion(j,2) = ryuion(j,1)
           rzuion(j,2) = rzuion(j,1)
           qquion(j,2) = qquion(j,1)
        end if
      end do

      rxu(pick_chain,pick_unit) = rxu(pick_chain,pick_unit) + rx
      ryu(pick_chain,pick_unit) = ryu(pick_chain,pick_unit) + ry
      rzu(pick_chain,pick_unit) = rzu(pick_chain,pick_unit) + rz

      moltion(2) = imolty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the energy of i in the new configuration
! Nonbonded is calculated in energy and is based on r*uion coords
! Bonded is calculated in U_bonded and is based on r*u coords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      flagon = 2     ! 2 refers to the coordinates for the NEW configuration
      call energy(i,imolty,vn,flagon,ibox,pick_unit,pick_unit,.true.,ovrlap,.false.,.false.,.false.,.true.)
      if (ovrlap)  then
            rxu(pick_chain,pick_unit) = rxu(pick_chain,pick_unit) - rx
            ryu(pick_chain,pick_unit) = ryu(pick_chain,pick_unit) - ry
            rzu(pick_chain,pick_unit) = rzu(pick_chain,pick_unit) - rz
            return
      end if

      call U_bonded(i,imolty,vn(ivStretching),vn(ivBending),vn(ivTorsion))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the the recip stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( lewald .and. lelect(imolty) ) then
         call recip_atom(ibox,vrecipn,vrecipo,1,pick_unit)
         vn(ivElect) = vn(ivElect) + vrecipn + vn(ivEwald)
         vo(ivElect) = vo(ivElect) + vrecipo + vo(ivEwald)
         vn(ivTot) = vn(ivTot) + vrecipn
         vo(ivTot) = vo(ivTot) + vrecipo
      end if

! check for acceptance ***

      vn(ivTot) = vn(ivTot) + vn(ivStretching) + vn(ivBending) + vn(ivTorsion)  ! We were missing this part in the last version! Meaning the bonded
      vo(ivTot) = vo(ivTot) + vo(ivStretching) + vo(ivBending) + vo(ivTorsion)  ! potential wasn't being factored into the acceptance criteria!!!

      deltv  = vn(ivTot) - vo(ivTot)
      deltvb = beta * deltv

! For ANES algorithm, do the Fluctuating charge moves.
! For time being it will not work for atom disp [Neeraj]****
      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vn,vo,vdum,vdum,vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) Abstrax = Abstrax + 1.0E0_dp
            if (ly) Abstray = Abstray + 1.0E0_dp
            if (lz) Abstraz = Abstraz + 1.0E0_dp
         end if
         return
      end if

      if ( deltvb .gt. (2.3E0_dp*softcut) ) then
                ! move rejected
            rxu(pick_chain,pick_unit) = rxu(pick_chain,pick_unit) - rx
            ryu(pick_chain,pick_unit) = ryu(pick_chain,pick_unit) - ry
            rzu(pick_chain,pick_unit) = rzu(pick_chain,pick_unit) - rz
            return
      end if

      if ( deltv .le. 0.0E0_dp ) then
               ! accept move
      else if ( exp(-deltvb) .gt. random(-1) ) then
               ! accept move
      else
               ! move rejected
            rxu(pick_chain,pick_unit) = rxu(pick_chain,pick_unit) - rx
            ryu(pick_chain,pick_unit) = ryu(pick_chain,pick_unit) - ry
            rzu(pick_chain,pick_unit) = rzu(pick_chain,pick_unit) - rz
            return
      end if

! write(io_output,*) 'TRANSLATION accepted i',i
      vbox(ivTot,ibox)     = vbox(ivTot,ibox) + deltv
      vbox(ivInterLJ,ibox)  = vbox(ivInterLJ,ibox) + (vn(ivInterLJ) - vo(ivInterLJ))
      vbox(ivIntraLJ,ibox)  = vbox(ivIntraLJ,ibox) + (vn(ivIntraLJ) - vo(ivIntraLJ))
      vbox(ivExt,ibox)    = vbox(ivExt,ibox)   + (vn(ivExt)   - vo(ivExt))
      vbox(ivElect,ibox)   = vbox(ivElect,ibox)  + (vn(ivElect) - vo(ivElect))

      vbox(ivTorsion,ibox) = vbox(ivTorsion,ibox) + (vn(ivTorsion)-vo(ivTorsion))
      vbox(ivBending,ibox) = vbox(ivBending,ibox) + (vn(ivBending)-vo(ivBending))
      vbox(ivStretching,ibox) = vbox(ivStretching,ibox) + (vn(ivStretching)-vo(ivStretching))


      if (lx) rxu(i,pick_unit) = rxuion(pick_unit,2)
      if (ly) ryu(i,pick_unit) = ryuion(pick_unit,2)
      if (lz) rzu(i,pick_unit) = rzuion(pick_unit,2)

      if (lewald .and. lelect(imolty)) then
! update reciprocal-space sum
         call recip_atom(ibox,vdum,vdum,2,pick_unit)
      end if

      if ( ldielect ) then
! update the dipole term
         call dipole(ibox,1)
      end if


! update chain center of mass

      call ctrmas(.false.,ibox,i,10)

      if (lx) Abstrax = Abstrax + 1.0E0_dp
      if (ly) Abstray = Abstray + 1.0E0_dp
      if (lz) Abstraz = Abstraz + 1.0E0_dp

      if ( licell .and. (ibox.eq.boxlink)) then
! update linkcell list
         call update_linked_cell(i)
      end if

      if (lneigh.or.lneighbor) then
         call update_neighbor_list_molecule(i)
      end if

#ifdef __DEBUG__
      write(io_output,*) 'end ATOM_TRANSLATION in ',myid,i
#endif
      return
    end subroutine Atom_translation

  subroutine init_moves_simple(io_input,lprint)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::jerr,i
    real::armtra,rmtra,rmrot
    namelist /mc_simple/ armtra,pm_atom_tra,tatra,pmtra,pmtrmt,rmtra,tarot,pmromt,rmrot

    if (allocated(acntrax)) deallocate(acntrax,acntray,acntraz,acnrotx,acnroty,acnrotz,acstrax,acstray,acstraz,acsrotx,acsroty,acsrotz,bntrax,bntray,bntraz,bstrax,bstray,bstraz,bnrotx,bnroty,bnrotz,bsrotx,bsroty,bsrotz,stat=jerr)
    allocate(acntrax(ntmax,nbxmax),acntray(ntmax,nbxmax),acntraz(ntmax,nbxmax),acnrotx(ntmax,nbxmax),acnroty(ntmax,nbxmax)&
     ,acnrotz(ntmax,nbxmax),acstrax(ntmax,nbxmax),acstray(ntmax,nbxmax),acstraz(ntmax,nbxmax),acsrotx(ntmax,nbxmax)&
     ,acsroty(ntmax,nbxmax),acsrotz(ntmax,nbxmax),bntrax(ntmax,nbxmax),bntray(ntmax,nbxmax),bntraz(ntmax,nbxmax)&
     ,bstrax(ntmax,nbxmax),bstray(ntmax,nbxmax),bstraz(ntmax,nbxmax),bnrotx(ntmax,nbxmax),bnroty(ntmax,nbxmax)&
     ,bnrotz(ntmax,nbxmax),bsrotx(ntmax,nbxmax),bsroty(ntmax,nbxmax),bsrotz(ntmax,nbxmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_moves_simple: allocation failed',jerr)

    Abstrax = 0.0E0_dp
    Abstray = 0.0E0_dp
    Abstraz = 0.0E0_dp
    Abntrax = 0.0E0_dp
    Abntray = 0.0E0_dp
    Abntraz = 0.0E0_dp
    acntrax = 0.E0_dp
    acntray = 0.E0_dp
    acntraz = 0.E0_dp
    acnrotx = 0.E0_dp
    acnroty = 0.E0_dp
    acnrotz = 0.E0_dp
    acstrax = 0.E0_dp
    acstray = 0.E0_dp
    acstraz = 0.E0_dp
    acsrotx = 0.E0_dp
    acsroty = 0.E0_dp
    acsrotz = 0.E0_dp
    bstrax = 0.0E0_dp
    bstray = 0.0E0_dp
    bstraz = 0.0E0_dp
    bsrotx = 0.0E0_dp
    bsroty = 0.0E0_dp
    bsrotz = 0.0E0_dp
    bntrax = 0.0E0_dp
    bntray = 0.0E0_dp
    bntraz = 0.0E0_dp
    bnrotx = 0.0E0_dp
    bnroty = 0.0E0_dp
    bnrotz = 0.0E0_dp

    !> read namelist mc_simple
    do i=1,nmolty
       pmtrmt(i)=real(i,dp)/nmolty
       pmromt(i)=real(i,dp)/nmolty
    end do
    armtra=0.01_dp
    rmtra=0.3_dp
    rmrot=0.4_dp

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=mc_simple,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_simple',jerr)
    end if

    call mp_bcast(armtra,1,rootid,groupid)
    call mp_bcast(pm_atom_tra,1,rootid,groupid)
    call mp_bcast(tatra,1,rootid,groupid)
    call mp_bcast(pmtra,1,rootid,groupid)
    call mp_bcast(pmtrmt,nmolty,rootid,groupid)
    call mp_bcast(rmtra,1,rootid,groupid)
    call mp_bcast(tarot,1,rootid,groupid)
    call mp_bcast(pmromt,nmolty,rootid,groupid)
    call mp_bcast(rmrot,1,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_SIMPLE','------------------------------------------'
       write(io_output,'(A,G16.9)') 'pm_atom_tra: ',pm_atom_tra
       write(io_output,'(A,G16.9)') 'initial maximum displacement for atom translation: ',armtra
       write(io_output,'(A,F4.2)') 'target translation acceptance ratio: ',tatra
       write(io_output,'(A,G16.9)') 'initial maximum displacement for molecule translation: ',rmtra
       write(io_output,'(A,G16.9)') 'pmtra: ',pmtra
       do i=1,nmolty
          write(io_output,'(A,I0,A,F8.4)') '   translation probability for molecule type ',i,' (pmtrmt): ',pmtrmt(i)
       end do
       write(io_output,'(A,F4.2)') 'target rotation acceptance ratio: ',tarot
       write(io_output,'(A,G16.9)') 'initial maximum displacement for molecule rotation: ',rmrot
       write(io_output,'(A,G16.9)') 'pmrot: ',1.0E0_dp
       do i=1,nmolty
          write(io_output,'(A,I0,A,f8.4)') '   rotational probability for molecule type ',i,' (pmromt): ',pmromt(i)
       end do
    end if

    Armtrax=armtra
    Armtray=armtra
    Armtraz=armtra
    rmtrax=rmtra
    rmtray=rmtra
    rmtraz=rmtra
    rmrotx=rmrot
    rmroty=rmrot
    rmrotz=rmrot
  end subroutine init_moves_simple

  subroutine update_translation_rotation_max_displacement(io_output)
    use sim_particle,only:check_neighbor_list
    integer,intent(in)::io_output
    real::rcutmin,ratrax,ratray,ratraz,rttrax,rttray,rttraz,rarotx,raroty,rarotz,rtrotx,rtroty,rtrotz
    integer::im,imolty

    ! adjust the atomic displacements
    rcutmin = rcut(1)
    do im = 2,nbox
       if (rcutmin.gt.rcut(im)) then
          rcutmin = rcut(im)
       end if
    end do

    if ( Abntrax .gt. 0.5E0_dp ) then
       ratrax = Abstrax / Abntrax
       rttrax = Armtrax * ratrax / tatra
       if ( rttrax .gt. 2.0E0_dp * rcutmin ) then
          ! maximum translational displacement
          Armtrax = 2.0E0_dp*rcutmin
       else if (rttrax .lt. 1.0E-10_dp ) then
          ! ratio must have been zero, so divide by 10
          Armtrax = Armtrax/10.0E0_dp
       else
          ! accept new ratio
          Armtrax = rttrax
       end if
    end if

    if ( Abntray .gt. 0.5E0_dp ) then
       ratray = Abstray / Abntray
       rttray = Armtray * ratray / tatra
       if ( rttray .gt. 2.0E0_dp * rcutmin ) then
          ! maximum translational displacement
          Armtray = 2.0E0_dp*rcutmin
       else if (rttray .lt. 1.0E-10_dp ) then
          ! ratio must have been zero, so divide by 10
          Armtray = Armtray/10.0E0_dp
       else
          ! accept new ratio
          Armtray = rttray
       end if
    end if

    if ( Abntraz .gt. 0.5E0_dp ) then
       ratraz = Abstraz / Abntraz
       rttraz = Armtraz * ratraz / tatra
       if ( rttraz .gt. 2.0E0_dp * rcutmin ) then
          ! maximum translational displacement
          Armtraz = 2.0E0_dp*rcutin
       else if (rttraz .lt. 1.0E-10_dp ) then
          ! ratio must have been zero, so divide by 10
          Armtraz = Armtraz/10.0E0_dp
       else
          ! accept new ratio
          Armtraz = rttraz
       end if
    end if

    ! adjust maximum translational displacement for COM***
    do im=1,nbox
       if (myid.eq.rootid) then
          write(io_output,*) 'Box ',im
       end if
       do imolty = 1,nmolty
          ! rmtrax
          ! check whether any x translations have been done for
          ! molecule type imolty in box im
          if ( bntrax(imolty,im) .gt. 0.5E0_dp ) then
             ! compute possible new ratio
             ratrax = bstrax(imolty,im) / bntrax(imolty,im) ! percent acceptance
             rttrax = rmtrax(imolty,im) * ratrax / tatra

             if ( rttrax .gt. 2.0E0_dp * rcut(im)) then
                ! maximum translational displacement
                rmtrax(imolty,im) = 2.0E0_dp*rcut(im)
             else if (rttrax .lt. 1.0E-10_dp ) then
                ! ratio must have been zero, so divide by 10
                rmtrax(imolty,im) = rmtrax(imolty,im)/10.0E0_dp
             else
                ! accept new ratio
                rmtrax(imolty,im) = rttrax
             end if
          end if

          ! rmtray
          ! check whether any y translations have been done for
          ! molecule type imolty in box im
          if ( bntray(imolty,im) .gt. 0.5E0_dp ) then
             ! compute possiblt new ratio
             ratray = bstray(imolty,im) / bntray(imolty,im)
             rttray = rmtray(imolty,im) * ratray / tatra

             if ( rttray .gt. 2.0E0_dp*rcut(im) ) then
                ! maximum translational displacement
                rmtray(imolty,im) = 2.0E0_dp*rcut(im)
             else if (rttray .eq. 0.0E0_dp) then
                ! ratio must have been zero, divide old by 10
                rmtray(imolty,im) = rmtray(imolty,im)/10.0E0_dp
             else
                ! accept new ratio
                rmtray(imolty,im) = rttray
             end if
          end if

          ! rmtraz
          ! check whether any z translations have been done for
          ! molecule type imolty in box im
          if ( bntraz(imolty,im) .gt. 0.5E0_dp ) then
             ! compute possible new ratio
             ratraz = bstraz(imolty,im) / bntraz(imolty,im)
             rttraz = rmtraz(imolty,im) * ratraz / tatra

             if ( rttraz .gt. 2.0E0_dp*rcut(im) ) then
                ! maximum translational displacement
                rmtraz(imolty,im) = 2.0E0_dp*rcut(im)
             else if ( rttraz .lt. 1.0E-10_dp ) then
                ! ratio must have been zero, divide old by 10
                rmtraz(imolty,im) = rmtraz(imolty,im)/10.0E0_dp
             else
                ! accept new ratio
                rmtraz(imolty,im) = rttraz
             end if
          end if

          ! rmrotx
          ! check whether any x-axis rotations have been done for
          ! molecule type imolty in box im
          if ( bnrotx(imolty,im) .gt. 0.5E0_dp ) then
             ! compute possible new ratio
             rarotx = bsrotx(imolty,im) / bnrotx(imolty,im)
             rtrotx = rmrotx(imolty,im) * rarotx / tarot

             if ( rtrotx .lt. 1.0E-10_dp )  then
                ! ratio was zero, divide old by 10
                rmrotx(imolty,im) = rmrotx(imolty,im)/10.E0_dp
             else if (rtrotx .gt. 3.1415E0_dp) then
                ! maximum rotational displacement
                rmrotx(imolty,im) = 3.1415E0_dp
             else
                ! accept trial ratio
                rmrotx(imolty,im) = rtrotx
             end if
          end if

          ! rmroty
          ! check whether any y-axis rotations have been done for
          ! molecule type imolty in box im
          if ( bnroty(imolty,im) .gt. 0.5E0_dp ) then
             ! compute possible new ratio
             raroty = bsroty(imolty,im) / bnroty(imolty,im)
             rtroty = rmroty(imolty,im) * raroty / tarot

             if ( rtroty .lt. 1.0E-10_dp)  then
                ! ratio was zero, divide old by 10
                rmroty(imolty,im) = rmroty(imolty,im)/10.0E0_dp
             else if (rtroty .gt. 3.1415E0_dp) then
                ! maximum rotational displacement
                rmroty(imolty,im) = 3.1415E0_dp
             else
                ! accept trial ratio
                rmroty(imolty,im) = rtroty
             end if
          end if

          ! rmrotz
          ! check whether any y-axis rotations have been done for
          ! molecule type imolty in box im
          if ( bnrotz(imolty,im) .gt. 0.5E0_dp ) then
             ! compute possible new ratio
             rarotz = bsrotz(imolty,im) / bnrotz(imolty,im)
             rtrotz = rmrotz(imolty,im) * rarotz / tarot

             if (rtrotz .eq. 0.0E0_dp)  then
                ! ratio was zero, divide old by 10
                rmrotz(imolty,im) = rmrotz(imolty,im)/10.0E0_dp
             else if (rtrotz .gt. 3.1415E0_dp) then
                ! maximum rotational displacement
                rmrotz(imolty,im) = 3.1415E0_dp
             else
                ! accept trial ratio
                rmrotz(imolty,im) = rtrotz
             end if
          end if

          if (lneigh) call check_neighbor_list(im,imolty)
          call averageMaximumDisplacement(im,imolty)

          ! KM for MPI
          ! only processor 0 writes to output files (except for error messages)
          if (myid.eq.rootid) then
             ! write some ratio update information ***
             write(io_output,"(' Type ',i2,' bn',6(1x,f7.0),/,9x,'bs',6(1x,f7.0))") imolty,bntrax(imolty,im),bntray(imolty,im)&
              ,bntraz(imolty,im),bnrotx(imolty,im),bnroty(imolty,im),bnrotz(imolty,im),bstrax(imolty,im),bstray(imolty,im)&
              ,bstraz(imolty,im),bsrotx(imolty,im),bsroty(imolty,im),bsrotz(imolty,im)
             ! write(io_output,"(' ratio       ',6(f7.4,1x))") ratrax, ratray, ratraz,rarotx, raroty,rarotz
             write(io_output,"(' max.displ. ',6(1x,f11.4))") rmtrax(imolty,im),rmtray(imolty,im),rmtraz(imolty,im),rmrotx(imolty,im)&
              ,rmroty(imolty,im),rmrotz(imolty,im)
          end if

          acntrax(imolty,im) = acntrax(imolty,im)+bntrax(imolty,im)
          acntray(imolty,im) = acntray(imolty,im)+bntray(imolty,im)
          acntraz(imolty,im) = acntraz(imolty,im)+bntraz(imolty,im)
          acstrax(imolty,im) = acstrax(imolty,im)+bstrax(imolty,im)
          acstray(imolty,im) = acstray(imolty,im)+bstray(imolty,im)
          acstraz(imolty,im) = acstraz(imolty,im)+bstraz(imolty,im)

          acnrotx(imolty,im) = acnrotx(imolty,im)+bnrotx(imolty,im)
          acnroty(imolty,im) = acnroty(imolty,im)+bnroty(imolty,im)
          acnrotz(imolty,im) = acnrotz(imolty,im)+bnrotz(imolty,im)
          acsrotx(imolty,im) = acsrotx(imolty,im)+bsrotx(imolty,im)
          acsroty(imolty,im) = acsroty(imolty,im)+bsroty(imolty,im)
          acsrotz(imolty,im) = acsrotz(imolty,im)+bsrotz(imolty,im)

          bstrax(imolty,im) = 0.0E0_dp
          bstray(imolty,im) = 0.0E0_dp
          bstraz(imolty,im) = 0.0E0_dp
          bntrax(imolty,im) = 0.0E0_dp
          bntray(imolty,im) = 0.0E0_dp
          bntraz(imolty,im) = 0.0E0_dp
          bsrotx(imolty,im) = 0.0E0_dp
          bsroty(imolty,im) = 0.0E0_dp
          bsrotz(imolty,im) = 0.0E0_dp
          bnrotx(imolty,im) = 0.0E0_dp
          bnroty(imolty,im) = 0.0E0_dp
          bnrotz(imolty,im) = 0.0E0_dp
       end do
    end do
  end subroutine update_translation_rotation_max_displacement

  subroutine averageMaximumDisplacement(ibox,imol)
    integer,intent(in)::ibox,imol

    if (numberDimensionIsIsotropic(ibox).eq.3) then
       rmtrax(imol,ibox)=(rmtrax(imol,ibox)+rmtray(imol,ibox)+rmtraz(imol,ibox))/3
       rmtray(imol,ibox)=rmtrax(imol,ibox)
       rmtraz(imol,ibox)=rmtrax(imol,ibox)
       rmrotx(imol,ibox)=(rmrotx(imol,ibox)+rmroty(imol,ibox)+rmrotz(imol,ibox))/3
       rmroty(imol,ibox)=rmrotx(imol,ibox)
       rmrotz(imol,ibox)=rmrotx(imol,ibox)
    else if (numberDimensionIsIsotropic(ibox).eq.2) then
       rmtrax(imol,ibox)=(rmtrax(imol,ibox)+rmtray(imol,ibox))/2
       rmtray(imol,ibox)=rmtrax(imol,ibox)
       rmrotx(imol,ibox)=(rmrotx(imol,ibox)+rmroty(imol,ibox))/2
       rmroty(imol,ibox)=rmrotx(imol,ibox)
    end if
  end subroutine averageMaximumDisplacement

!> \brief Write some information about translations and rotations
  subroutine output_translation_rotation_stats(io_output)
    use util_string,only:format_n
    integer,intent(in)::io_output
    integer::ibox,i,ii
    real::ratvol
    character(LEN=default_path_length)::fmt(3)=(/"(' x-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"&
     ,"(' y-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"&
     ,"(' z-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"/)

    if (pm_atom_tra.ge.1E-6) then
        write(io_output,*)
        write(io_output,*) '### Atom Translations ###'
        write(io_output,*)
        if ( Abstrax .ne. 0.0E0_dp ) then
            ratvol = Abstrax / Abntrax
        else
            ratvol = 0.0E0_dp
        end if
        write(io_output,fmt(1)) Abntrax,ratvol,Armtrax
        if ( Abstrax .ne. 0.0E0_dp ) then
            ratvol = Abstray / Abntray
        else
            ratvol = 0.0E0_dp
        end if
        write(io_output,fmt(1)) Abntray,ratvol,Armtray
        if ( Abstrax .ne. 0.0E0_dp ) then
            ratvol = Abstraz / Abntraz
        else
            ratvol = 0.0E0_dp
        end if
        write(io_output,fmt(1)) Abntraz,ratvol,Armtraz
    end if

    write(io_output,*)
    write(io_output,*) '### Translations ###'
    write(io_output,*)
    do ibox = 1,nbox
       do i=1,nmolty
          write(io_output,"(A,I0,A)",advance='no') 'molecule typ = ',i,'   '
          write(io_output,'(A10)',advance='no')molecname(i)
          write(io_output,*)' in box',ibox
          acntrax(i,ibox) = acntrax(i,ibox) + bntrax(i,ibox)
          acstrax(i,ibox) = acstrax(i,ibox) + bstrax(i,ibox)
          if ( acntrax(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acstrax(i,ibox) / acntrax(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(1)) acntrax(i,ibox),ratvol,rmtrax(i,ibox)

          acntray(i,ibox) = acntray(i,ibox) + bntray(i,ibox)
          acstray(i,ibox) = acstray(i,ibox) + bstray(i,ibox)
          if ( acntray(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acstray(i,ibox) / acntray(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(2)) acntray(i,ibox),ratvol,rmtray(i,ibox)

          acntraz(i,ibox) = acntraz(i,ibox) + bntraz(i,ibox)
          acstraz(i,ibox) = acstraz(i,ibox) + bstraz(i,ibox)
          if ( acntraz(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acstraz(i,ibox) / acntraz(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(3)) acntraz(i,ibox),ratvol,rmtraz(i,ibox)
          write(io_output,*)
       end do
    end do

    write(io_output,*) '### Rotations ###'
    write(io_output,*)
    do ibox = 1,nbox
       do i=1,nmolty
          write(io_output,"(A,I0,A)",advance='no') 'molecule typ = ',i,'   '
          write(io_output,'(A10)',advance='no')molecname(i)
          write(io_output,*)' in box',ibox
          acnrotx(i,ibox) = acnrotx(i,ibox) + bnrotx(i,ibox)
          acsrotx(i,ibox) = acsrotx(i,ibox) + bsrotx(i,ibox)
          if ( acnrotx(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acsrotx(i,ibox) / acnrotx(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(1)) acnrotx(i,ibox),ratvol,rmrotx(i,ibox)

          acnroty(i,ibox) = acnroty(i,ibox) + bnroty(i,ibox)
          acsroty(i,ibox) = acsroty(i,ibox) + bsroty(i,ibox)
          if ( acnroty(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acsroty(i,ibox) / acnroty(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(2)) acnroty(i,ibox),ratvol,rmroty(i,ibox)

          acnrotz(i,ibox) = acnrotz(i,ibox) + bnrotz(i,ibox)
          acsrotz(i,ibox) = acsrotz(i,ibox) + bsrotz(i,ibox)
          if ( acnrotz(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acsrotz(i,ibox) / acnrotz(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(3)) acnrotz(i,ibox),ratvol,rmrotz(i,ibox)
          write(io_output,*)
       end do
    end do
  end subroutine output_translation_rotation_stats

  subroutine read_checkpoint_simple(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) Abntrax,Abntray,Abntraz,Abstrax,Abstray,Abstraz,bntrax,bntray,bntraz,bstrax,bstray,bstraz&
     ,acntrax,acntray,acntraz,acstrax,acstray,acstraz,bnrotx,bnroty,bnrotz,bsrotx,bsroty,bsrotz,acnrotx,acnroty,acnrotz,acsrotx&
     ,acsroty,acsrotz
    call mp_bcast(Abntrax,1,rootid,groupid)
    call mp_bcast(Abntray,1,rootid,groupid)
    call mp_bcast(Abntraz,1,rootid,groupid)
    call mp_bcast(Abstrax,1,rootid,groupid)
    call mp_bcast(Abstray,1,rootid,groupid)
    call mp_bcast(Abstraz,1,rootid,groupid)
    call mp_bcast(bntrax,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bntray,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bntraz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bstrax,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bstray,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bstraz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acntrax,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acntray,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acntraz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acstrax,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acstray,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acstraz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bnrotx,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bnroty,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bnrotz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsrotx,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsroty,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsrotz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acnrotx,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acnroty,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acnrotz,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acsrotx,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acsroty,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(acsrotz,ntmax*nbxmax,rootid,groupid)
  end subroutine read_checkpoint_simple

  subroutine write_checkpoint_simple(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) Abntrax,Abntray,Abntraz,Abstrax,Abstray,Abstraz,bntrax,bntray,bntraz,bstrax,bstray,bstraz,acntrax,acntray,acntraz,acstrax,acstray,acstraz,bnrotx,bnroty,bnrotz,bsrotx,bsroty,bsrotz,acnrotx,acnroty,acnrotz,acsrotx,acsroty,acsrotz
  end subroutine write_checkpoint_simple
end MODULE moves_simple
