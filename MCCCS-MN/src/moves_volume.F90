MODULE moves_volume
  use util_random,only:random
  use util_runtime,only:err_exit
  use util_kdtree,only:construct_kdtree, scale_kdtree
  use sim_system
  use sim_cell
  use energy_kspace,only:recip,calp,save_kvector,restore_kvector
  use energy_pairwise,only:sumup
  implicit none
  private
  save
  public::volume_1box,volume_2box,init_moves_volume,update_volume_max_displacement,output_volume_stats,read_checkpoint_volume&
   ,write_checkpoint_volume,allow_cutoff_failure,restore_displ_transl

  real,allocatable,public::acsvol(:),acnvol(:),acshmat(:,:),acnhmat(:,:),bsvol(:),bnvol(:),bshmat(:,:),bnhmat(:,:),acc_displ(:)
  real,allocatable::vboxn(:,:),vboxo(:,:),bxo(:),byo(:),bzo(:),xcmo(:),ycmo(:),zcmo(:),rxuo(:,:),ryuo(:,:),rzuo(:,:),qquo(:,:)&
   ,rcut_original(:)
  real::hmato(9),hmatio(9)
  integer::allow_cutoff_failure=-1 !< controls how volume move failures, the ones that will result in box lengths
  !< smaller than twice the cutoff, are handled: -1 = fetal error and program exits;
  !< 0 = simply rejects the move, which is equivalent to modifying the lower limit of integration of the partition function;
  !< 1 = allows the move and adjusts cutoff to be half of the new box lengths (will be restored if possible);
  !< 2 = allows the move but does not adjust the cutoff, which could be problematic because this results in a lower density
  !< in the cutoff radius (due to periodic boundary conditions)

contains
!> \brief Makes an isotropic volume change for NVT-Gibbs ensemble
!>
!> Perform change of the volume: random walk in ln(V1/V2) with V1+V2=const. \n
!> The maximum change is controlled by \b rmtrax and the
!> number of successful trial moves is stored in \b bsvol.
  subroutine volume_2box()
    real::rpair,rm,rbox,volo(nbxmax),volt,voln(nbxmax),rbcut(nbxmax),dfac(nbxmax),df,dx,dy,dz,expdv,min_boxl,v(nEnergy),dele,displ
    integer::ipair,ipairb,boxa,boxb,ibox,i,hbox,jbox,jhmat,imolty,j,ichoiq
    logical::lncubic,lx(nbxmax),ly(nbxmax),lz(nbxmax),ovrlap,ladjust
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start VOLUME_2BOX in ',myid
#endif

    ! select pair of boxes to do the volume move
    if ( nvolb .gt. 1 ) then
       rpair = random(-1)
       do ipair = 1, nvolb
          if ( rpair .lt. pmvolb(ipair) ) then
             ipairb = ipair
             exit
          end if
       end do
    else
       ipairb = 1
    end if
    boxa = box5(ipairb)
    boxb = box6(ipairb)

    bnvol(ipairb) = bnvol(ipairb) + 1.0E0_dp
    displ = 0.0E0_dp !< the accepted displacement

    call save_box(boxa)
    call save_box(boxb)
    call save_configuration((/boxa,boxb/))

    lncubic = .false.
    lx = .false.
    ly = .false.
    lz = .false.
    do ibox = 1, 2
       if (ibox .eq. 1) i = boxa
       if (ibox .eq. 2) i = boxb

       if (lsolid(i)) then
          ! volume move independently in x, y, z directions
          rm = random(-1)
          if ( rm .le. pmvolx ) then
             lx(i) = .true.
          else if ( rm .le. pmvoly ) then
             ly(i) = .true.
          else
             lz(i) = .true.
          end if

          if (.not. lrect(i)) then
             lncubic = .true.
             hbox = i
             volo(i) = cell_vol(i)

             ! select one of the cell edge
             if ( lx(i) ) then
                rbox = 3.0_dp*random(-1)
                if ( rbox .lt. 1.0E0_dp ) then
                   jhmat = 1
                else if (rbox .lt. 2.0E0_dp ) then
                   jhmat = 4
                else
                   jhmat = 7
                end if
             else if ( ly(i) ) then
                rbox = 2.0_dp*random(-1)
                if ( rbox .lt. 1.0E0_dp ) then
                   jhmat = 5
                else
                   jhmat = 8
                end if
             else
                jhmat = 9
             end if
          end if
       end if

       if (.not.lsolid(i).or.lrect(i)) then
          jbox = i
          if ( lpbcz ) then
             volo(i) = bxo(i)*byo(i)*bzo(i)
          else
             volo(i) = bxo(i)*byo(i)
          end if
       end if
    end do

    ! calculate total volume
    volt = volo(boxa) + volo(boxb)

    if ( lncubic ) then
       ! hbox is the non-orthorhombic box, jbox is the other box
       bnhmat(hbox,jhmat) = bnhmat(hbox,jhmat) + 1.0E0_dp
       hmat(hbox,jhmat) = hmat(hbox,jhmat) + rmhmat(hbox,jhmat)*( 2.0E0_dp*random(-1) - 1.0E0_dp )
       call matops(hbox)

       voln(hbox) = cell_vol(hbox)
       voln(jbox) = volt-voln(hbox)

       if (lsolid(jbox)) then
          ! volume move independently in x, y, z directions
          dfac(jbox)=voln(jbox)/volo(jbox)
          if (lx(jbox)) boxlx(jbox) = boxlx(jbox) * dfac(jbox)
          if (ly(jbox)) boxly(jbox) = boxly(jbox) * dfac(jbox)
          if (lz(jbox)) boxlz(jbox) = boxlz(jbox) * dfac(jbox)
       else
          if ( lpbcz ) then
             dfac(jbox)= (voln(jbox)/volo(jbox))**(1.0E0_dp/3.0E0_dp)
             boxlz(jbox) = boxlz(jbox)*dfac(jbox)
          else
             dfac(jbox)= sqrt(voln(jbox)/volo(jbox))
          end if
          boxlx(jbox) = boxlx(jbox)*dfac(jbox)
          boxly(jbox) = boxly(jbox)*dfac(jbox)
       end if

       if (allow_cutoff_failure.ne.2) then
          rbcut(hbox) = minval(min_width(hbox,:))/2.0_dp
          rbcut(jbox) = min(boxlx(jbox),boxly(jbox))/2.0_dp
          if (lpbcz) rbcut(jbox)=min(rbcut(jbox),boxlz(jbox)/2.0_dp)
          if (allow_cutoff_failure.eq.1) then
             ladjust=.false.
             if (rbcut(hbox).lt.rcut(hbox)) then
                ladjust=.true.
             else if (rcut(hbox).lt.rcut_original(hbox)) then
                ladjust=.true.
                if (rbcut(hbox).gt.rcut_original(hbox)) rbcut(hbox)=rcut_original(hbox)
             end if
             if (ladjust) then
                rbox=rbcut(hbox)
                rbcut(hbox)=rcut(hbox) ! set rbcut back to old cutoff
                rcut(hbox)=rbox ! update cutoff to half min boxlength
             else
                rbcut(hbox)=-1.0_dp
             end if

             ladjust=.false.
             if (rbcut(jbox).lt.rcut(jbox)) then
                ladjust=.true.
             else if (rcut(jbox).lt.rcut_original(jbox)) then
                ladjust=.true.
                if (rbcut(jbox).gt.rcut_original(jbox)) rbcut(jbox)=rcut_original(jbox)
             end if
             if (ladjust) then
                rbox=rbcut(jbox)
                rbcut(jbox)=rcut(jbox) ! set rbcut back to old cutoff
                rcut(jbox)=rbox        ! update cutoff to half min boxlength
             else
                rbcut(jbox)=-1.0_dp
             end if
          else if (rbcut(hbox).lt.rcut(hbox).or.rbcut(jbox).lt.rcut(jbox)) then
             hmat(hbox,jhmat) = hmato(jhmat)
             boxlx(jbox) = bxo(jbox)
             boxly(jbox) = byo(jbox)
             if ( lpbcz ) then
                boxlz(jbox) = bzo(jbox)
             end if
             if (allow_cutoff_failure.eq.-1) then
                call dump('final-config')
                write(io_output,*) 'w1:',min_width(hbox,1),'w2:',min_width(hbox,2),'w3:',min_width(hbox,3)
                call err_exit(__FILE__,__LINE__,'non-orthorhombic volume_2box move rejected. box width below cutoff size',myid+1)
             else if (allow_cutoff_failure.eq.0) then
                call matops(hbox)
                return
             end if
          end if
       end if

       ! determine the displacement of the COM
       df = dfac(jbox) - 1.0E0_dp
       do i = 1,nchain
          ibox = nboxi(i)
          imolty = moltyp(i)
          if (ibox .eq. hbox) then
             if ( lx(ibox) ) then
                dx = sxcm(i)*(hmat(hbox,1)-hmato(1))+sycm(i)*(hmat(hbox,4)-hmato(4))+szcm(i)*(hmat(hbox,7)-hmato(7))
                xcm(i) = xcm(i) + dx
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                end do
             else if ( ly(ibox) ) then
                dy = sycm(i)*(hmat(hbox,5)-hmato(5))+szcm(i)*(hmat(hbox,8)-hmato(8))
                ycm(i) = ycm(i) + dy
                do j = 1, nunit(imolty)
                   ryu(i,j) = ryu(i,j) + dy
                end do
             else
                dz = szcm(i)*(hmat(hbox,9)-hmato(9))
                zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          else if (ibox .eq. jbox) then
             if (lsolid(jbox)) then
                if ( lx(ibox) ) then
                   dx = xcm(i) * df
                   xcm(i) = xcm(i) + dx
                   do j = 1, nunit(imolty)
                      rxu(i,j) = rxu(i,j) + dx
                   end do
                else if ( ly(ibox) ) then
                   dy = ycm(i) * df
                   ycm(i) = ycm(i) + dy
                   do j = 1, nunit(imolty)
                      ryu(i,j) = ryu(i,j) + dy
                   end do
                else
                   dz = zcm(i) * df
                   zcm(i) = zcm(i) + dz
                   do j = 1, nunit(imolty)
                      rzu(i,j) = rzu(i,j) + dz
                   end do
                end if
             else
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
             end if
          end if
       end do
    else
       ! calculate new volume
       expdv=rmvol(ipairb)*(2.0E0_dp*random(-1)-1.0E0_dp)
       displ = abs(expdv) !< magnitude of the volume displacement
       expdv = volo(boxa)/volo(boxb)*exp(expdv)
       voln(boxa)= expdv*volt/(1+expdv)
       voln(boxb)= volt-voln(boxa)

       do i=1,2
          if (i.eq.1) ibox=boxa
          if (i.eq.2) ibox=boxb

          if (lsolid(ibox)) then
             ! volume move independently in x, y, z directions
             dfac(ibox)=voln(ibox)/volo(ibox)
             if (lx(ibox)) boxlx(ibox) = boxlx(ibox) * dfac(ibox)
             if (ly(ibox)) boxly(ibox) = boxly(ibox) * dfac(ibox)
             if (lz(ibox)) boxlz(ibox) = boxlz(ibox) * dfac(ibox)
          else
             if ( lpbcz ) then
                dfac(ibox)= (voln(ibox)/volo(ibox))**(1.0E0_dp/3.0E0_dp)
                boxlz(ibox) = boxlz(ibox) * dfac(ibox)
             else
                dfac(ibox)= sqrt(voln(ibox)/volo(ibox))
             end if
             boxlx(ibox) = boxlx(ibox) * dfac(ibox)
             boxly(ibox) = boxly(ibox) * dfac(ibox)
          end if
       end do

       if (allow_cutoff_failure.ne.2) then
          ! initialize rbcut to half minimum of new boxlength
          rbcut(boxa) = min(boxlx(boxa),boxly(boxa))/2.0_dp
          rbcut(boxb) = min(boxlx(boxb),boxly(boxb))/2.0_dp
          if (lpbcz) then
             rbcut(boxa)=min(rbcut(boxa),boxlz(boxa)/2.0_dp)
             rbcut(boxb)=min(rbcut(boxb),boxlz(boxb)/2.0_dp)
          end if
          ! determine what to do with rbcut
          if (allow_cutoff_failure.eq.1) then
             ladjust=.false.
             if (rbcut(boxa).lt.rcut(boxa)) then
                ladjust=.true.
             else if (rcut(boxa).lt.rcut_original(boxa)) then
                ladjust=.true.
                if (rbcut(boxa).gt.rcut_original(boxa)) rbcut(boxa)=rcut_original(boxa)
             end if
             if (ladjust) then
                rbox=rbcut(boxa)
                rbcut(boxa)=rcut(boxa) ! set rbcut back to old rcut
                rcut(boxa)=rbox ! update cutoff to half min boxlength
             else
                rbcut(boxa)=-1.0_dp
             end if

             ladjust=.false.
             if (rbcut(boxb).lt.rcut(boxb)) then
                ladjust=.true.
             else if (rcut(boxb).lt.rcut_original(boxb)) then
                ladjust=.true.
                if (rbcut(boxb).gt.rcut_original(boxb)) rbcut(boxb)=rcut_original(boxb)
             end if
             if (ladjust) then
                rbox=rbcut(boxb)
                rbcut(boxb)=rcut(boxb) ! set rbcut back to old cutoff
                rcut(boxb)=rbox        ! update cutoff to half min boxlength
             else
                rbcut(boxb)=-1.0_dp
             end if
          else if (rbcut(boxa).lt.rcut(boxa).or.rbcut(boxb).lt.rcut(boxb)) then
             boxlx(boxa) = bxo(boxa)
             boxlx(boxb) = bxo(boxb)
             boxly(boxa) = byo(boxa)
             boxly(boxb) = byo(boxb)
             if ( lpbcz ) then
                boxlz(boxa) = bzo(boxa)
                boxlz(boxb) = bzo(boxb)
             end if
             if (allow_cutoff_failure.eq.-1) then
                call dump('final-config')
                call err_exit(__FILE__,__LINE__,'A move was attempted that would lead to a boxlength less than twice rcut',myid+1)
             else if (allow_cutoff_failure.eq.0) then
                return
             end if
          end if
       end if

       ! determine new positions of the molecules
       ! calculate centre of mass and its displacement
       do i = 1, nchain
          ibox = nboxi(i)
          if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
             imolty = moltyp(i)
             df = dfac(ibox) - 1.0E0_dp
             if (lsolid(ibox)) then
                if ( lx(ibox) ) then
                   dx = xcm(i) * df
                   xcm(i) = xcm(i) + dx
                   do j = 1, nunit(imolty)
                      rxu(i,j) = rxu(i,j) + dx
                   end do
                else if ( ly(ibox) ) then
                   dy = ycm(i) * df
                   ycm(i) = ycm(i) + dy
                   do j = 1, nunit(imolty)
                      ryu(i,j) = ryu(i,j) + dy
                   end do
                else
                   dz = zcm(i) * df
                   zcm(i) = zcm(i) + dz
                   do j = 1, nunit(imolty)
                      rzu(i,j) = rzu(i,j) + dz
                   end do
                end if
             else
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
             end if
          end if
       end do
    end if

    if ( lchgall ) then
       if (lsolid(boxa).and.(.not.lrect(boxa))) then
          min_boxl = minval(min_width(boxa,:))
       else
          min_boxl = min(boxlx(boxa),boxly(boxa),boxlz(boxa))
       end if
       calp(boxa) = kalp(boxa)/min_boxl
       if (lsolid(boxb).and.(.not.lrect(boxb))) then
          min_boxl = minval(min_width(boxb,:))
       else
          min_boxl = min(boxlx(boxb),boxly(boxb),boxlz(boxb))
       end if
       calp(boxb) = kalp(boxb)/min_boxl
    end if

    ! create kdtree for the sumup
    if (lkdtree) then
        ! construct a kdtree for the fictitious box nbox+1
        ! when sumup is called from the volume move, the energy of nbox+1 and nbox+2 will be calculated instead of ibox
        ! when lcutcm (COM-kdtree), simply scale the coordinates
        if (lkdtree_box(boxa)) then
            if (lcutcm) then
                call scale_kdtree(boxa, dfac(boxa))
            else
                call construct_kdtree(boxa, nbox+1, .false.)
            end if
        end if

        if (lkdtree_box(boxb)) then
            if (lcutcm) then
                call scale_kdtree(boxb, dfac(boxb))
            else
                call construct_kdtree(boxb, nbox+2, .false.)
            end if
        end if
    end if

    do i = 1,2
       if ( i .eq. 1 ) ibox = boxa
       if ( i .eq. 2 ) ibox = boxb
       call sumup(ovrlap,v,ibox,.true.)
       if ( ovrlap ) goto 500
       vboxn(:,ibox) = v
       vboxn(ivTot,ibox) = vboxo(ivTot,ibox) + (vboxn(ivInterLJ,ibox)-vboxo(ivInterLJ,ibox)) + (vboxn(ivExt,ibox)-vboxo(ivExt,ibox)) + (vboxn(ivElect,ibox)-vboxo(ivElect,ibox)) + (vboxn(iv3body,ibox)-vboxo(iv3body,ibox)) ! inter, ext, elect, garo
    end do

    if ( lanes ) then
       ! for ANES algorithm, optimize the charge configuration
       ! on the new coordinates, continue to use the fluctuating charge
       ! algorithm to optimize the charge configurations, update the
       ! energy, coordinates and the ewald sum
       do i = 1,2
          if ( i .eq. 1 ) ibox = boxa
          if ( i .eq. 2 ) ibox = boxb
          vbox(ivTot,ibox) = vbox(ivTot,ibox) + (vboxn(ivTot,ibox) - vboxo(ivTot,ibox))
          vbox(ivInterLJ,ibox)  = vbox(ivInterLJ,ibox) +  (vboxn(ivInterLJ,ibox) - vboxo(ivInterLJ,ibox))
          vbox(ivTail,ibox) = vbox(ivTail,ibox) + (vboxn(ivTail,ibox) - vboxo(ivTail,ibox))
          vbox(ivExt,ibox) = vbox(ivExt,ibox) + (vboxn(ivExt,ibox) - vboxo(ivExt,ibox))
          vbox(ivElect,ibox) = vbox(ivElect,ibox) +  (vboxn(ivElect,ibox) - vboxo(ivElect,ibox))
          do ichoiq = 1,nchoiq(ibox)
             call flucq(0,ibox)
          end do
       end do
       dele = (vbox(ivTot,boxa) - vboxo(ivTot,boxa))+( vbox(ivTot,boxb)- vboxo(ivTot,boxb)) - ((nchbox(boxa)+1+ghost_particles(boxa)) *log(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+1+ghost_particles(boxb)) *log(voln(boxb)/volo(boxb))/beta)
    else if (lncubic) then
       dele = (vboxn(ivTot,boxa)-vboxo(ivTot,boxa)) + (vboxn(ivTot,boxb)-vboxo(ivTot,boxb)) - ((nchbox(boxa)+ghost_particles(boxa))*log(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+ghost_particles(boxb))*log(voln(boxb)/volo(boxb))/beta)
    else
       dele = (vboxn(ivTot,boxa)-vboxo(ivTot,boxa)) + (vboxn(ivTot,boxb)-vboxo(ivTot,boxb)) - ((nchbox(boxa)+1+ghost_particles(boxa))*log(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+1+ghost_particles(boxb))*log(voln(boxb)/volo(boxb))/beta)
    end if

    ! allows pressure difference (osmotic pressure)
    if (losmoticnvt) then
       dele = dele + (express(boxa)-express(boxb)) * (voln(boxa)-volo(boxa))
    end if

    ! acceptance test
    if (random(-1) .lt. exp(-beta*dele) ) then
       ! accepted
       bsvol(ipairb) = bsvol(ipairb) + 1.0E0_dp
       acc_displ(ipairb) = acc_displ(ipairb) + displ
       if ( lncubic ) then
          bshmat(hbox,jhmat) = bshmat(hbox,jhmat) + 1.0E0_dp
       end if
       call update_box(boxa)
       call update_box(boxb)
       if (allow_cutoff_failure.eq.1) then
          ! if rbcut>0, we needed to change the cutoff. In this case we need to
          ! make sure that the maximum displacements are not too large
          if (rbcut(boxa).gt.0) then
             call restore_displ_transl(boxa)
             if (L_Ewald_Auto) then
                ! new kvectors needed
                kalp(boxa) = 3.2_dp/rcut(boxa)
                calp(boxa) = kalp(boxa)
             end if
          end if
          if (rbcut(boxb).gt.0) then
             call restore_displ_transl(boxb)
             if (L_Ewald_Auto) then
                kalp(boxb) = 3.2_dp/rcut(boxb)
                calp(boxb) = kalp(boxb)
             end if
          end if
       end if
       return
    end if

    ! rejected
500 call restore_box(boxa)
    call restore_box(boxb)
    call restore_configuration((/boxa,boxb/))

    if (allow_cutoff_failure.eq.1) then
       ! if rejected, rbcut contains the cutoff before volume move.
       ! we need to restore that value
       if (rbcut(boxa).gt.0) rcut(boxa) = rbcut(boxa)
       if (rbcut(boxb).gt.0) rcut(boxb) = rbcut(boxb)
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end VOLUME_2BOX in ',myid,boxa,boxb
#endif
    return
  end subroutine volume_2box

!> \brief Makes an isotropic volume change under const. pressure.
!>
!> Perform change of the volume: random walk in V. \n
!> The maximum change is controlled by \b rmtrax and the
!> number of successful trial moves is stored in \b bsvol.
  subroutine volume_1box()
    real::rbox,volo,voln,rbcut,dx,dy,dz,dfac,df,v(nEnergy),dele,min_boxl,displ
    integer::ibox,boxvch,jhmat,i,imolty,j,ichoiq,nchain_boxvch
    logical::lx,ly,lz,ovrlap,ladjust,l_couple,l_consv
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start VOLUME_1BOX in ',myid
#endif
    ! Select a box at  random to change the volume of box
    rbox = random(-1)
    do ibox = 1,nbox
       if (rbox .lt. pmvlmt(ibox) ) then
          boxvch=ibox
          exit
       end if
    end do

    bnvol(boxvch) = bnvol(boxvch) + 1.0E0_dp
    displ = 0.0E0_dp

    call save_box(boxvch)
    call save_configuration((/boxvch/))

    nchain_boxvch = 0

    lx = .false.
    ly = .false.
    lz = .false.

    ! Is this going to be a bilayer move? If so, what kind?

    l_consv  = .false.
    l_couple = .false.

    if(boxvch.eq.1.and.l_bilayer) then
       if ((pm_consv.gt.0).and.(random(-1).lt.pm_consv)) then
             l_consv = .true.
       else
             l_couple = .true.
       end if
    end if


    if ( lsolid(boxvch) ) then
       ! volume move independently in x, y, z directions
       rbox = random(-1)
       if(l_couple) then
          if(rbox.le.pmvol_xy) then
             lx = .true.
             ly = .true.
          else
             lz = .true.
          end if
       else if (l_consv) then
             lx = .true.
             ly = .true.
             lz = .true.
       else
          if ( rbox .le. pmvolx ) then
             lx = .true.
          else if ( rbox .le. pmvoly ) then
             ly = .true.
          else
             lz = .true.
          end if
       end if

       if (.not.lrect(boxvch)) then
          volo = cell_vol(boxvch)
          ! select one of the cell edge
          if ( .not.l_bilayer ) then
             if ( lx ) then
                rbox = 3.0_dp*random(-1)
                if ( rbox .le. 1.0E0_dp ) then
                   jhmat = 1
                else if (rbox .le. 2.0E0_dp ) then
                   jhmat = 4
                else
                   jhmat = 7
                end if
             else if ( ly ) then
                rbox = 2.0_dp*random(-1)
                if ( rbox .le. 1.0E0_dp ) then
                   jhmat = 5
                else
                   jhmat = 8
                end if
             else
                jhmat = 9
             end if
          end if
       end if
    end if

    if (.not.lsolid(boxvch).or.lrect(boxvch)) then
       if ( lpbcz ) then
          volo = bxo(boxvch)*byo(boxvch)*bzo(boxvch)
       else
          volo = bxo(boxvch)*byo(boxvch)
       end if
    end if

    if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then

       if(l_couple) then ! Case 1: a coupled xy bilayer move
          if(lx.and.ly) then
             ! Update attempts
             bnhmat(boxvch,1) = bnhmat(boxvch,1) + 1.0E0_dp
             bnhmat(boxvch,5) = bnhmat(boxvch,5) + 1.0E0_dp

             ! Do a random move in boxlength using x and apply the equivalent change to y
             hmat(boxvch,1) = hmat(boxvch,1) + rmhmat(boxvch,1)* ( 2.0E0_dp*random(-1) - 1.0E0_dp )
             hmat(boxvch,5) = hmat(boxvch,1)
          else
             bnhmat(boxvch,9) = bnhmat(boxvch,9) + 1.0E0_dp
             hmat(boxvch,9) = hmat(boxvch,9) + rmhmat(boxvch,9)* ( 2.0E0_dp*random(-1) - 1.0E0_dp )
          end if
       else if(l_consv) then ! Case 2: a volume conserving bilayer move
          ! Update attempts
          bnhmat(boxvch,1) = bnhmat(boxvch,1) + 1.0E0_dp
          bnhmat(boxvch,5) = bnhmat(boxvch,5) + 1.0E0_dp
          bnhmat(boxvch,9) = bnhmat(boxvch,9) + 1.0E0_dp

          ! Do a random move in boxlength using x and apply the equivalent change to y
          hmat(boxvch,1) = hmat(boxvch,1) + rmhmat(boxvch,1)* ( 2.0E0_dp*random(-1) - 1.0E0_dp )
          hmat(boxvch,9) =  hmat(boxvch,5)* hmat(boxvch,5)/hmat(boxvch,1)/hmat(boxvch,1)*hmat(boxvch,9) ! since x and y are equivalent, this is saying that [(x^2_old)/(X^2_new)]*z_old = z_new
          hmat(boxvch,5) = hmat(boxvch,1)
       else ! Case 3: everything else
          bnhmat(boxvch,jhmat) = bnhmat(boxvch,jhmat) + 1.0E0_dp
          hmat(boxvch,jhmat) = hmat(boxvch,jhmat) + rmhmat(boxvch,jhmat)* ( 2.0E0_dp*random(-1) - 1.0E0_dp )
       end if

       call matops(boxvch)

       voln = cell_vol(boxvch)

       if (allow_cutoff_failure.ne.2) then
          rbcut = minval(min_width(boxvch,:))/2.0_dp
          if (allow_cutoff_failure.eq.1) then
             ladjust=.false.
             if (rbcut.lt.rcut(boxvch)) then
                ladjust=.true.
             else if (rcut(boxvch).lt.rcut_original(boxvch)) then
                ladjust=.true.
                if (rbcut.gt.rcut_original(boxvch)) rbcut=rcut_original(boxvch)
             end if
             if (ladjust) then
                rbox=rbcut
                rbcut=rcut(boxvch) ! set rbcut back to old cutoff
                rcut(boxvch)=rbox  ! update cutoff to half min boxlength
             else
                rbcut=-1.0_dp
             end if
          else if (rbcut.lt.rcut(boxvch)) then

             if(l_couple) then
                if(lx.and.ly) then
                   hmat(boxvch,1) = hmato(1)
                   hmat(boxvch,5) = hmato(5)
                else
                   hmat(boxvch,9) = hmato(9)
                end if
             else if (l_consv) then
                hmat(boxvch,1) = hmato(1)
                hmat(boxvch,4) = hmato(5) !possibly an error. Check before using l_consv
                hmat(boxvch,5) = hmato(9) ! same as line above
             else
                hmat(boxvch,jhmat) = hmato(jhmat)
             end if
             if (allow_cutoff_failure.eq.-1) then
                call dump('final-config')
                write(io_output,*) 'w1:',min_width(boxvch,1),'w2:',min_width(boxvch,2),'w3:',min_width(boxvch,3)
                call err_exit(__FILE__,__LINE__,'non-rectangular volume move rejected. box width below cutoff size',myid+1)
             else if (allow_cutoff_failure.eq.0) then
                call matops(boxvch)
                return
             end if
          end if
       end if

       ! determine the displacement of the COM
       do i = 1,nchain
          if (nboxi(i) .eq. boxvch) then
             nchain_boxvch = nchain_boxvch + 1
             imolty = moltyp(i)
             if ( lx ) then
                dx = sxcm(i)*(hmat(boxvch,1)-hmato(1))+sycm(i)*(hmat(boxvch,4)-hmato(4))+szcm(i)*(hmat(boxvch,7)-hmato(7))
                xcm(i) = xcm(i) + dx
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                end do
             end if
             if ( ly ) then
                dy = sycm(i)*(hmat(boxvch,5)-hmato(5))+szcm(i)*(hmat(boxvch,8)-hmato(8))
                ycm(i) = ycm(i) + dy
                do j = 1, nunit(imolty)
                   ryu(i,j) = ryu(i,j) + dy
                end do
             end if
             if ( lz ) then
                dz = szcm(i)*(hmat(boxvch,9)-hmato(9))
                zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          end if
       end do
    else
       ! calculate new volume
       displ = rmvol(boxvch) * ( 2.0E0_dp*random(-1) - 1.0E0_dp )
       voln = volo + displ
       displ = abs(displ)

       if (lsolid(boxvch)) then
          ! volume move independently in x, y, z directions
          dfac=voln/volo
          if (lx) boxlx(boxvch) = boxlx(boxvch) * dfac
          if (ly) boxly(boxvch) = boxly(boxvch) * dfac
          if (lz) boxlz(boxvch) = boxlz(boxvch) * dfac
       else
          if ( lpbcz ) then
             dfac = (voln/volo)**(1.0E0_dp/3.0E0_dp)
             boxlz(boxvch) = boxlz(boxvch) * dfac
          else
             dfac= sqrt(voln/volo)
          end if
          boxlx(boxvch) = boxlx(boxvch) * dfac
          boxly(boxvch) = boxly(boxvch) * dfac
       end if

       if (allow_cutoff_failure.ne.2) then
          ! determine rbcut: max possible rcut = 2*mininumBoxlength
          rbcut = min(boxlx(boxvch),boxly(boxvch))/2.0_dp
          if (lpbcz) rbcut=min(rbcut,boxlz(boxvch)/2.0_dp)
          if (allow_cutoff_failure.eq.1) then
             ladjust=.false.
             if (rbcut.lt.rcut(boxvch)) then
                ladjust=.true.
             else if (rcut(boxvch).lt.rcut_original(boxvch)) then
                ! must have made previous change to rcut
                ladjust=.true.
                if (rbcut.gt.rcut_original(boxvch)) then
                   ! box size has increased so new rbcut predicted to be >
                   ! original. Set rbcut back to original so rcut is then
                   ! changed back to original
                   rbcut=rcut_original(boxvch)
                end if
             end if
             if (ladjust) then
                ! we need to adjust rcut
                rbox=rbcut
                rbcut=rcut(boxvch) ! set rbcut back to old cutoff
                rcut(boxvch)=rbox  ! update cutoff to half min boxlength
             else
                rbcut=-1.0_dp !later used as flag to designate that there is no need
                              !to change rcut to rbcut
             end if
          else if (rbcut.lt.rcut(boxvch)) then
             boxlx(boxvch) = bxo(boxvch)
             boxly(boxvch) = byo(boxvch)
             if ( lpbcz ) then
                boxlz(boxvch) = bzo(boxvch)
             end if
             if (allow_cutoff_failure.eq.-1) then
                call dump('final-config')
                write(io_output,*) 'boxvch',boxvch
                call err_exit(__FILE__,__LINE__,'A move was attempted that would lead to a boxlength less than twice rcut',myid+1)
             else if (allow_cutoff_failure.eq.0) then
                return !just reject the move
             end if
          end if
       end if

       ! determine new positions of the molecules
       ! calculate centre of mass and its displacement
       df = dfac - 1.0E0_dp
       do i = 1, nchain
          ! Check if the chain i is in the correct box
          if (nboxi(i) .eq. boxvch) then
             nchain_boxvch = nchain_boxvch + 1
             imolty = moltyp(i)
             if (lsolid(boxvch)) then
                if ( lx ) then
                   dx = xcm(i) * df
                   xcm(i) = xcm(i) + dx
                   do j = 1, nunit(imolty)
                      rxu(i,j) = rxu(i,j) + dx
                   end do
                else if ( ly ) then
                   dy = ycm(i) * df
                   ycm(i) = ycm(i) + dy
                   do j = 1, nunit(imolty)
                      ryu(i,j) = ryu(i,j) + dy
                   end do
                else
                   dz = zcm(i) * df
                   zcm(i) = zcm(i) + dz
                   do j = 1, nunit(imolty)
                      rzu(i,j) = rzu(i,j) + dz
                   end do
                end if
             else
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
             end if
          end if
       end do
    end if

    if (nchain_boxvch.eq.0) then
       ! if no molecules in box
       if (allow_cutoff_failure.lt.0) then
          ! if default, error exit
          call err_exit(__FILE__,__LINE__,'Volume 1box attempted with no molec in box', myid+1)
       else
          ! allow_cutoff_failure was set to account for this
          ! reject move and do not count this as an attempt
          bnvol(boxvch) = bnvol(boxvch) - 1.0E0_dp
          goto 500
       end if
    end if


    if ( lchgall ) then
       if (lsolid(boxvch).and.(.not.lrect(boxvch))) then
          min_boxl = minval(min_width(boxvch,:))
       else
          min_boxl = min(boxlx(boxvch),boxly(boxvch),boxlz(boxvch))
       end if
       calp(boxvch) = kalp(boxvch)/boxlx(boxvch)
    end if

    ! create kdtree for the sumup
    if (lkdtree .and. lkdtree_box(ibox)) then
        ! construct a kdtree for the fictitious box nbox+1
        ! when sumup is called from the volume move, the energy of nbox+1 will be calculated instead of ibox
        if (lcutcm) then
            call scale_kdtree(ibox, dfac)
        else
            call construct_kdtree(ibox, nbox+1, .false.)
        end if
    end if

    call sumup(ovrlap,v,boxvch,.true.)
    if ( ovrlap ) goto 500
    vboxn(:,boxvch) = v
    vboxn(ivTot,boxvch) = vboxo(ivTot,boxvch) + (vboxn(ivInterLJ,boxvch)-vboxo(ivInterLJ,boxvch)) + (vboxn(ivExt,boxvch)-vboxo(ivExt,boxvch)) + (vboxn(ivElect,boxvch)-vboxo(ivElect,boxvch)) + (vboxn(iv3body,boxvch)-vboxo(iv3body,boxvch)) !inter, ext, elect, garo

    if ( lanes ) then
       ! for ANES algorithm, optimize the charge configuration
       ! on the new coordinates, continue to use the fluctuating charge
       ! algorithm to optimize the charge configurations, update the
       ! energy, coordinates and the ewald sum
       vbox(ivTot,boxvch)=vbox(ivTot,boxvch)+(vboxn(ivTot,boxvch)-vboxo(ivTot,boxvch))
       vbox(ivInterLJ,boxvch) = vbox(ivInterLJ,boxvch) + (vboxn(ivInterLJ,boxvch)-vboxo(ivInterLJ,boxvch))
       vbox(ivTail,boxvch)  = vbox(ivTail,boxvch) + (vboxn(ivTail,boxvch)-vboxo(ivTail,boxvch))
       vbox(ivExt,boxvch)   = vbox(ivExt,boxvch) + (vboxn(ivExt,boxvch)-vboxo(ivExt,boxvch))
       vbox(ivElect,boxvch) = vbox(ivElect,boxvch) + (vboxn(ivElect,boxvch)-vboxo(ivElect,boxvch))
       do ichoiq = 1,nchoiq(boxvch)
          call flucq(0,boxvch)
       end do
       dele = (vbox(ivTot,boxvch) - vboxo(ivTot,boxvch)) + express(boxvch)*(voln-volo) - ((nchbox(boxvch)+ghost_particles(boxvch)) * log(voln/volo) / beta )
    else
       dele = ( vboxn(ivTot,boxvch) - vboxo(ivTot,boxvch) ) + express(boxvch)*(voln-volo) - ((nchbox(boxvch)+ghost_particles(boxvch))*log(voln/volo)/beta)
    end if

    ! acceptance test
    if (random(-1) .lt. exp(-(beta*dele)) ) then
       ! accepted
       bsvol(boxvch) = bsvol(boxvch) + 1.0E0_dp
       acc_displ(boxvch) = acc_displ(boxvch) + displ
       if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
          if(l_couple) then
             if(lx.and.ly) then
                bshmat(boxvch,1) = bshmat(boxvch,1) + 1.0E0_dp
                bshmat(boxvch,5) = bshmat(boxvch,5) + 1.0E0_dp
             else
                bshmat(boxvch,9) = bshmat(boxvch,9) + 1.0E0_dp
             end if
          else if (l_consv) then
             bshmat(boxvch,1) = bshmat(boxvch,1) + 1.0E0_dp
             bshmat(boxvch,5) = bshmat(boxvch,5) + 1.0E0_dp
             bshmat(boxvch,9) = bshmat(boxvch,9) + 1.0E0_dp
          else
             bshmat(boxvch,jhmat) = bshmat(boxvch,jhmat) + 1.0E0_dp
          end if
       end if
       call update_box(boxvch)
       if (allow_cutoff_failure.eq.1) then
          ! if rbcut>0, we needed to change the cutoff. In this case we need to
          ! make sure that the maximum displacements are not too large
          ! and potentially update k_max
          if (rbcut.gt.0) then
             call restore_displ_transl(boxvch)
             if (L_Ewald_Auto) then
                kalp(boxvch) = 3.2_dp/rcut(boxvch)
                calp(boxvch) = kalp(boxvch)
             end if
          end if
       end if
       return
    end if

    ! rejected
500 call restore_box(boxvch)
    call restore_configuration((/boxvch/))

    ! if rejected, restore rcut
    if (allow_cutoff_failure.eq.1) then
       ! restore cutoff if needed.
       ! if rbcut > 0, a move was attempted that would make rcut < a
       ! boxlength/2. rbcut now contains original value, so we must change rcut
       ! back to the original value
       if (rbcut.gt.0) rcut(boxvch) = rbcut
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end VOLUME_1BOX in ',myid,boxvch
#endif
    return
  end subroutine volume_1box

  subroutine init_moves_volume(io_input,lprint)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::jerr,i,j,k
    real::rmvolume
    namelist /mc_volume/ tavol,iratv,pmvlmt,nvolb,pmvolb,box5,box6,pmvol,pmvolx,pmvoly,rmvolume,allow_cutoff_failure,l_bilayer,pm_consv, pmvol_xy

    if (allocated(acsvol)) deallocate(acsvol,acnvol,acshmat,acnhmat,bsvol,bnvol,acc_displ,bshmat,bnhmat,vboxn,vboxo,bxo,byo,bzo,xcmo,ycmo,zcmo,rxuo,ryuo,rzuo,qquo,rcut_original,stat=jerr)
    allocate(acsvol(nbxmax),acnvol(nbxmax),acshmat(nbxmax,9),acnhmat(nbxmax,9),bsvol(nbxmax),bnvol(nbxmax),bshmat(nbxmax,9)&
     ,bnhmat(nbxmax,9),vboxn(nEnergy,nbxmax),vboxo(nEnergy,nbxmax),bxo(nbxmax),byo(nbxmax),bzo(nbxmax),xcmo(nmax),ycmo(nmax)&
     ,zcmo(nmax),rxuo(nmax,numax),ryuo(nmax,numax),rzuo(nmax,numax),qquo(nmax,numax),rcut_original(nbxmax),acc_displ(nbxmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_moves_volume: allocation failed',jerr)

    ! defaults for namelist mc_volume
    acsvol = 0.E0_dp
    acnvol = 0.E0_dp
    acshmat = 0.0E0_dp
    acnhmat = 0.0E0_dp
    bsvol = 0.0E0_dp
    bnvol = 0.0E0_dp
    acc_displ = 0.0E0_dp
    bshmat = 0.0E0_dp
    bnhmat = 0.0E0_dp
    rcut_original = rcut
    l_bilayer = .false.
    pm_consv = 0.0d0
    pmvol_xy = 0.0d0

    nvolb=nbox*(nbox-1)/2
    do i=1,nvolb
       pmvolb(i)=real(i,dp)/nvolb
    end do
    k=0

    do i=1,nbox
       pmvlmt(i)=-1 !dummy integer, set default after read if needed
       do j=i+1,nbox
          k=k+1
          box5(k)=i
          box6(k)=j
       end do
    end do
    if (lnpt) then
       rmvolume=1.0E3_dp
    else
       rmvolume=1.0E-3_dp
    end if

    ! end defaults for namelist mc_volume
    !> read namelist mc_volume
    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=mc_volume,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_volume',jerr)
    end if


    if (ALL(pmvlmt(1:nbox).eq.-1)) then
       ! no pmvlmt provided; give defaults
       do i=1,nbox
          pmvlmt(i)=real(i,dp)/nbox
       end do
       if (lexzeo) pmvlmt(1) = 0.0E0_dp
    else if (ANY(pmvlmt(1:nbox).eq.-1)) then
       call err_exit(__FILE__,__LINE__,'error in pmvlmt: probabilities not &
                                        & specified for all boxes', myid+1)
    else if ((lexzeo) .and. (pmvlmt(1).gt.0.0E0_dp)) then
       call err_exit(__FILE__,__LINE__,'error in pmvlmt: you should not be &
                                        & doing volume moves on box 1 when lexzeo = T', myid+1)
    end if


    call mp_bcast(tavol,1,rootid,groupid)
    call mp_bcast(iratv,1,rootid,groupid)
    call mp_bcast(pmvlmt,nbox,rootid,groupid)
    call mp_bcast(nvolb,1,rootid,groupid)
    call mp_bcast(pmvolb,nvolb,rootid,groupid)
    call mp_bcast(box5,nvolb,rootid,groupid)
    call mp_bcast(box6,nvolb,rootid,groupid)
    call mp_bcast(pmvol,1,rootid,groupid)
    call mp_bcast(pmvolx,1,rootid,groupid)
    call mp_bcast(pmvoly,1,rootid,groupid)
    call mp_bcast(rmvolume,1,rootid,groupid)
    call mp_bcast(allow_cutoff_failure,1,rootid,groupid)

    rmvol=rmvolume

    if (.not.lfold.and.pmvol.gt.0)  call err_exit(__FILE__,__LINE__,'volume move only correct with folded coordinates',myid+1)
    if (ALL(allow_cutoff_failure.ne.(/-1,0,1,2/))) then
       call err_exit(__FILE__,__LINE__,'mc_volume: invalid value for allow_cutoff_failure = '//integer_to_string(allow_cutoff_failure),-1)
    end if

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_VOLUME','------------------------------------------'
       write(io_output,'(A,F4.2)') 'target volume acceptance ratio (tavol): ',tavol
       write(io_output,'(A,I0)') 'frequency to adjust maximum volume displacement: ',iratv
       write(io_output,'(A,F8.3)') 'initial maximum volume displacement (rmvol): ',rmvolume
       write(io_output,'(A,G16.9)') 'pmvol: ',pmvol
       if (ANY(lsolid(1:nbox))) write(io_output,'(2(A,G16.9))') 'pmvolx: ',pmvolx,', pmvoly: ',pmvoly
       do i = 1,nbox
          write(io_output,'(A,I0,A,G16.9)') '   pmvlmt for box ',i,': ',pmvlmt(i)
       end do
       write(io_output,'(A,I0)') 'nvolb: ',nvolb
       write(io_output,'(A,I0)') 'cutoff will be addressed with option:  ',allow_cutoff_failure
       do i = 1,nvolb
          write(io_output,'(3(A,I0),A,G16.9)') '   box pair ',i,': between ',box5(i),' and ',box6(i),',   pmvolb = ',pmvolb(i)
          if ((lsolid(box5(i)).and..not.lrect(box5(i))).and.(lsolid(box6(i)).and..not.lrect(box6(i)))) call err_exit(__FILE__,__LINE__,'can not perform volume move between two non-rectangular boxes',myid+1)
       end do
       if(l_bilayer) then
          write(io_output,*)'l_bilayer is true.'
          if(pm_consv.gt.0)write(io_output,*)'...Will perform volume conserving cell moves with a probability of: ',pm_consv
          if(pm_consv.lt.1) then
             write(io_output,*)'...Will perform coupled xy volume moves with a probability of:      ',1-pm_consv
             write(io_output,*)'...and the probability for coupled volume moves in xy, and z are:   ',pmvol_xy, 1-pmvol_xy
          end if

          if(.not.lsolid(1)) call err_exit(__FILE__,__LINE__,'mc_volume: If l_bilayer is true, lsolid must be true  for box 1.',-1)
          if(lrect(1))       call err_exit(__FILE__,__LINE__,'mc_volume: If l_bilayer is true, lrect  must be false for box 1.',-1)
       end if
    end if
  end subroutine init_moves_volume

!> \brief Adjust maximum volume displacement
  subroutine update_volume_max_displacement(io_output)
    integer,intent(in)::io_output
    integer::ibox,j
    real::ratvol

    do ibox = 1, nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             if ( bnhmat(ibox,j) .gt. 0.5E0_dp ) then
                ratvol = bshmat(ibox,j) / bnhmat(ibox,j)
                if (ratvol .eq. 0.0E0_dp) then
                   rmhmat(ibox,j) = rmhmat(ibox,j) * 0.1E0_dp
                else
                   rmhmat(ibox,j) = rmhmat(ibox,j)*ratvol/tavol
                end if
             end if
          end do
       else
          if ( bnvol(ibox) .gt. 0.5E0_dp ) then
             ratvol = bsvol(ibox) / bnvol(ibox)
             if ( ratvol .eq. 0.0E0_dp ) then
                rmvol(ibox) = rmvol(ibox) * 0.1E0_dp
             else
                rmvol(ibox) = rmvol(ibox) * ratvol / tavol
                if (rmvol(ibox).gt.(0.10E0_dp*boxlx(ibox)*boxly(ibox)*boxlz(ibox))) then
                   rmvol(ibox)=0.1E0_dp*(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
                end if
             end if
          end if
       end if
    end do

    if (myid.eq.rootid) then
       do ibox = 1, nbox
          if (lsolid(ibox) .and. .not. lrect(ibox)) then
             do j = 1,9
                write(io_output,"(' h-matrix change:  bn =',f8.1, '   bs =',f8.1,'   max.displ. =',e12.5)") bnhmat(ibox,j),bshmat(ibox,j), rmhmat(ibox,j)
             end do
          else
             write(io_output,"(' volume change:  bn =',f8.1, '   bs =',f8.1,'   max.displ. =',e12.5)") bnvol(ibox),bsvol(ibox),rmvol(ibox)
          end if
       end do
    end if

    do ibox = 1, nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
             acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
             bshmat(ibox,j) = 0.0E0_dp
             bnhmat(ibox,j) = 0.0E0_dp
          end do
       else
          acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
          acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
          bnvol(ibox) = 0.0E0_dp
          bsvol(ibox) = 0.0E0_dp
       end if
    end do
  end subroutine update_volume_max_displacement

  subroutine output_volume_stats(io_output)
    integer,intent(in)::io_output
    integer::ibox,j
    real::ratvol,acc_displ_avg
    character(LEN=default_path_length)::fmt="(' h-matrix attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"

    write(io_output,*)
    write(io_output,*) '### Volume change       ###'
    do ibox = 1,nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
             acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
             if ( acshmat(ibox,j) .gt. 0.5E0_dp) then
                write(io_output,fmt) acnhmat(ibox,j), acshmat(ibox,j)/acnhmat(ibox,j),rmhmat(ibox,j)
             else
                write(io_output,fmt) acnhmat(ibox,j),0.0E0_dp,rmhmat(ibox,j)
             end if
          end do
       else
          acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
          acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
          if ( acnvol(ibox) .ne. 0.0E0_dp ) then
             ratvol = acsvol(ibox) / acnvol(ibox)
             acc_displ_avg = acc_displ(ibox)/acsvol(ibox)
          else
             ratvol = 0.0E0_dp
             acc_displ_avg = 0.0E0_dp
          end if
          write(io_output,"(' attempts =',f8.1,'   ratio =',f6.3, 'max.displ. =',e11.4,'   avg.acc.displ. =',e11.4)") &
            acnvol(ibox),ratvol,rmvol(ibox),acc_displ_avg
       end if
    end do
  end subroutine output_volume_stats

!> \brief Store old box lengths and energy
  subroutine save_box(box)
    integer,intent(in)::box
    integer::j
    real::vdum

    vboxo(ivTot,box) = vbox(ivTot,box)
    vboxo(ivInterLJ,box) = vbox(ivInterLJ,box)
    vboxo(ivTail,box) = vbox(ivTail,box)
    vboxo(ivExt,box) = vbox(ivExt,box)
    vboxo(ivElect,box) = vbox(ivElect,box)
    vboxo(ivFlucq,box)= vbox(ivFlucq,box)
    vboxo(iv3body,box)= vbox(iv3body,box)

    bxo(box) = boxlx(box)
    byo(box) = boxly(box)
    if ( lpbcz ) bzo(box) = boxlz(box)

    if (lsolid(box) .and. .not. lrect(box)) then
       do j = 1,9
          hmato(j) = hmat(box,j)
          hmatio(j) = hmati(box,j)
       end do
    end if

    if ( lewald ) then
       call save_kvector(box)
       call recip(box,vdum,vdum,3)
    end if
  end subroutine save_box

!> \brief Store old chain configuration
  subroutine save_configuration(boxes)
    use sim_particle,only:save_neighbor_list
    integer,intent(in)::boxes(:)
    integer::i,j

    do i = 1, nchain
       if (ANY(nboxi(i).eq.boxes)) then
          xcmo(i) = xcm(i)
          ycmo(i) = ycm(i)
          if (lpbcz) zcmo(i) = zcm(i)
          do j = 1, nunit(moltyp(i))
             rxuo(i,j) = rxu(i,j)
             ryuo(i,j) = ryu(i,j)
             if ( lpbcz ) rzuo(i,j) = rzu(i,j)
             qquo(i,j) = qqu(i,j)
          end do
          if (lneighbor.or.lgaro) then
             call save_neighbor_list(i)
          end if
       end if
    end do
  end subroutine save_configuration

!> \brief Restore max displacement in translation
! this is only called if we are decreasing the cutoff to be half the boxlength,
! in which case we need to make sure the maximum displacement is not too large
! for this new cutoff
  subroutine restore_displ_transl(box)
    integer,intent(in)::box
    integer::imolty

    ! adjust maximum displacement in traslation if needed
    do imolty = 1,nmolty
       if (rmtrax(imolty,box) .gt. 2.0E0_dp*rcut(box)) rmtrax(imolty,box)=2.0E0_dp*rcut(box)
       if (rmtray(imolty,box) .gt. 2.0E0_dp*rcut(box)) rmtray(imolty,box)=2.0E0_dp*rcut(box)
       if (rmtraz(imolty,box) .gt. 2.0E0_dp*rcut(box)) rmtraz(imolty,box)=2.0E0_dp*rcut(box)
    end do
  end subroutine restore_displ_transl

!> \brief Restore old energy, box lengths
  subroutine restore_box(box)
    integer,intent(in)::box
    integer::j
    real::vdum, dfac

    vbox(ivTot,box) = vboxo(ivTot,box)
    vbox(ivInterLJ,box) = vboxo(ivInterLJ,box)
    vbox(ivTail,box) = vboxo(ivTail,box)
    vbox(ivExt,box) = vboxo(ivExt,box)
    vbox(ivElect,box) = vboxo(ivElect,box)
    vbox(ivFlucq,box)= vboxo(ivFlucq,box)
    vbox(iv3body,box)= vboxo(iv3body,box)

    dfac = bxo(box) / boxlx(box)

    ! rescale the COM-kdtree
    if (lcutcm .and. lkdtree .and. lkdtree_box(box)) call scale_kdtree(box, dfac)

    boxlx(box)   = bxo(box)
    boxly(box)   = byo(box)
    if ( lpbcz ) boxlz(box)   = bzo(box)

    if (lsolid(box) .and. .not. lrect(box)) then
       do j = 1,9
          hmat(box,j) = hmato(j)
          hmati(box,j) = hmatio(j)
       end do
       call matops(box)
    end if

    if ( lewald ) then
       call restore_kvector(box)
       call recip(box,vdum,vdum,4)
    end if
  end subroutine restore_box

!> \brief Restore old energy, box lengths
  subroutine restore_configuration(boxes)
    use sim_particle,only:restore_neighbor_list
    integer,intent(in)::boxes(:)
    integer::i,j

    do i = 1, nchain
       if (ANY(nboxi(i).eq.boxes)) then
          xcm(i) = xcmo(i)
          ycm(i) = ycmo(i)
          if ( lpbcz ) zcm(i) = zcmo(i)
          do j = 1, nunit(moltyp(i))
             rxu(i,j) = rxuo(i,j)
             ryu(i,j) = ryuo(i,j)
             if ( lpbcz ) rzu(i,j) = rzuo(i,j)
             qqu(i,j) = qquo(i,j)
          end do
          if (lneighbor.or.lgaro) then
             call restore_neighbor_list(i)
          end if
       end if
    end do

  end subroutine restore_configuration

  subroutine update_box(box)
    use sim_particle,only:ctrmas
    use util_kdtree,only:update_box_kdtree
    integer,intent(in)::box

    if ( .not. lanes ) then
       vbox(ivTot,box)    = vbox(ivTot,box) + (vboxn(ivTot,box) - vboxo(ivTot,box))
       vbox(ivInterLJ,box) = vbox(ivInterLJ,box) + (vboxn(ivInterLJ,box) - vboxo(ivInterLJ,box))
       vbox(ivTail,box)  = vbox(ivTail,box) + (vboxn(ivTail,box) - vboxo(ivTail,box))
       vbox(ivExt,box)   = vbox(ivExt,box) + (vboxn(ivExt,box) - vboxo(ivExt,box))
       vbox(ivElect,box) = vbox(ivElect,box) + (vboxn(ivElect,box) - vboxo(ivElect,box))
       vbox(iv3body,box) = vbox(iv3body,box) + (vboxn(iv3body,box)-vboxo(iv3body,box))
    end if

    ! update centers of mass
    call ctrmas(.true.,box,0,5)

    ! update coordinates in bead-kdtree
    if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(box)) then
        call update_box_kdtree(box)
    end if

    ! update linkcell, if applicable
    if (licell .and. (box .eq. boxlink)) then
       call build_linked_cell()
    end if
  end subroutine update_box

  subroutine read_checkpoint_volume(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) bnvol,bsvol,acnvol,acsvol,bnhmat,bshmat,acnhmat,acshmat
    call mp_bcast(bnvol,nbxmax,rootid,groupid)
    call mp_bcast(bsvol,nbxmax,rootid,groupid)
    call mp_bcast(acnvol,nbxmax,rootid,groupid)
    call mp_bcast(acsvol,nbxmax,rootid,groupid)
    call mp_bcast(bnhmat,nbxmax*9,rootid,groupid)
    call mp_bcast(bshmat,nbxmax*9,rootid,groupid)
    call mp_bcast(acnhmat,nbxmax*9,rootid,groupid)
    call mp_bcast(acshmat,nbxmax*9,rootid,groupid)
  end subroutine read_checkpoint_volume

  subroutine write_checkpoint_volume(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) bnvol,bsvol,acnvol,acsvol,bnhmat,bshmat,acnhmat,acshmat
  end subroutine write_checkpoint_volume
end MODULE moves_volume
