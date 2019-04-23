MODULE prop_pressure
  use var_type,only:dp
  use const_phys,only:qqfact,N_Avogadro
  use util_math,only:erfunc
  use util_mp,only:mp_sum
  use sim_system
  use sim_cell
  use energy_kspace,only:recippress,calp
  use energy_pairwise
  use util_kdtree,only:range_search
  implicit none
  private
  save
  public::pressure

contains
!> \brief Calculates the pressure for a configuration, unit [kPa]
!>
!> \warning Tabulated potential and Feuston-Garofalini
!> are NOT supported.
!> \note New potenial functional form needs to be added here,
!> and in U2 (energy calculation), coru (tail corrections to
!> energy), corp (tail corrections to pressure)
  subroutine pressure(press,surf,comp,ibox)
    use const_math,only:sqrtpi
    use const_phys,only:k_B
    use util_kdtree,only:range_search
    real,intent(out)::press,surf,comp
    integer,intent(in)::ibox

    real::rbcut,rcutsq,calpi,calpisq,pxx,pyy,pzz,xcmi,ycmi,zcmi,rcmi,fxcmi,fycmi,fzcmi,rxuij,ryuij,rzuij,rijsq,rcm,rcmsq&
     ,rxui,ryui,rzui,rij,fij,repress,rpxx,rpyy,rpzz&
     ,rpxy,rpyx,rpxz,rpzx,rpyz,rpzy,vol,volsq,rhosq,pwell
    integer::i,imolty,j,jmolty,ii,jj,ntii,ntjj,ntij,iii,jjj,k
    logical::lqimol,lexplt,lij2,lqjmol,lcoulo(numax,numax)

    ! Kdtree
    real, allocatable :: lij_list(:,:)
    logical :: ovrlap
! --------------------------------------------------------------------
    if (lsolid(ibox).and.(.not.lrect(ibox))) then
       vol = cell_vol(ibox)
    else
       vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
    end if

    if (lideal(ibox)) then
       press = k_B*1.0E27_dp * (nchbox(ibox)+ghost_particles(ibox))/beta/vol
       surf = 0.0_dp
       comp = 1.0_dp
       return
    end if

    if (lpbc) call setpbc(ibox)

    rbcut  = rcut(ibox)
    rcutsq = rbcut*rbcut
    calpi=calp(ibox)
    calpisq=calpi*calpi

    press = 0.0E0_dp
    pxx = 0.0E0_dp
    pyy = 0.0E0_dp
    pzz = 0.0E0_dp
    pips = 0.0E0_dp

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
    ! if(LSOLPAR.and.(ibox.eq.2)) then
    !    press= 1.380662E4_dp * ( ( nchbox(ibox) / beta) - ( press/3.0E0_dp ) ) / ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
    !    surf = 0.0E0_dp
    !    return
    ! end if

    ! loop over all chains i
    ! RP added for MPI
    do i = myid+1, nchain - 1,numprocs
       ! check if i is in relevant box ###
       if ( nboxi(i) .eq. ibox ) then
          if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
             call press_calc_kdtree(ibox, i, press, pxx, pyy, pzz)
          else
             imolty = moltyp(i)
             lqimol = lelect(imolty)
             if (nugrow(imolty).eq.nunit(imolty)) then
                 lexplt = .false.
             else
                 lexplt = .true.
             end if
             xcmi = xcm(i)
             ycmi = ycm(i)
             zcmi = zcm(i)
             if (lcutcm) then
                 rcmi = rcmu(i)
             else
                 lij2 = .true.
             end if

             if (lcutcm .and. lkdtree .and. lkdtree_box(ibox)) then
                  call energy_inter_com_kd_tree(ibox, i, xcmi, ycmi, zcmi, rbcut, 0.0 &
                       , .false., lij_list, .true., ovrlap)
             end if

             ! loop over all chains j with j>i
             do j = i + 1, nchain
                ! check for simulation box ###
                if ( nboxi(j) .eq. ibox ) then
                   jmolty = moltyp(j)
                   lqjmol = lelect(jmolty)
                   fxcmi = 0.0E0_dp
                   fycmi = 0.0E0_dp
                   fzcmi = 0.0E0_dp

                   ! Ctrmas cutoff
                   if (lcutcm .and. lkdtree .and. lkdtree_box(ibox)) then
                       if (lij_list(1, j) .lt. 0.5d0) then
                           cycle
                       else
                           lij2 = .true.
                           goto 1025
                       end if
                   end if

                   if ( lcutcm ) then
                      ! check if ctrmas within rcmsq
                      rxuij = xcmi - xcm(j)
                      ryuij = ycmi - ycm(j)
                      rzuij = zcmi - zcm(j)
                      ! minimum image the ctrmas pair separations
                      if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rcm = rbcut + rcmi + rcmu(j)
                      rcmsq = rcm*rcm
                      if (rijsq .le. rcmsq) then
                         lij2 = .true.
                      else if (lqimol .and. lqjmol .and. lchgall) then
                         lij2 = .false.
                      else
                         cycle
                      end if
                   end if

1025               continue

                   ! loop over all beads ii of chain i
                   do ii = 1, nunit(imolty)
                      ntii = ntype(imolty,ii)
                      rxui = rxu(i,ii)
                      ryui = ryu(i,ii)
                      rzui = rzu(i,ii)

                      ! loop over all beads jj of chain j
                      bead2: do jj = 1, nunit(jmolty)
                         ! check exclusion table
                         if (lexclu(imolty,ii,jmolty,jj)) cycle bead2

                         ntjj = ntype(jmolty,jj)
                         if ( lij2 ) then
                            if ((.not.(lij(ntii).and.lij(ntjj))) .and.(.not.(lqchg(ntii).and.lqchg(ntjj)))) cycle bead2
                         else
                            if (.not.(lqchg(ntii).and.lqchg(ntjj))) cycle bead2
                         end if

                         ntij=type_2body(ntii,ntjj)

                         rxuij = rxui - rxu(j,jj)
                         ryuij = ryui - ryu(j,jj)
                         rzuij = rzui - rzu(j,jj)
                         ! minimum image the pair separations ***
                         if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                         rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                         rij=sqrt(rijsq)

                         ! compute whether the charged groups will interact & fij
                         fij = 0.0E0_dp
                         if (lqimol.and.lqjmol.and.lqchg(ntii).and.lqchg(ntjj) ) then
                            if (.not.lewald) then
                               if (.not.lchgall) then
                                  iii = leaderq(imolty,ii)
                                  jjj = leaderq(jmolty,jj)
                                  if (iii.eq.ii.and.jjj.eq.jj) then
                                     ! set up the charge-interaction table
                                     if (rijsq.lt.rcutsq) then
                                        lcoulo(iii,jjj) = .true.
                                     else
                                        lcoulo(iii,jjj) = .false.
                                     end if
                                  end if
                               end if
                               !> \todo tabulated potential
                               if (lchgall.or.lcoulo(iii,jjj) ) then
                                  fij = -qqfact*qqu(i,ii)*qqu(j,jj)/rijsq/rij
                               end if
                            else if (lchgall.or.rijsq.lt.rcutsq) then
                                  fij = -qqfact*qqu(i,ii)*qqu(j,jj) &
                                        /rijsq*(2.0_dp*calpi*exp(-calpisq*rijsq)/sqrtpi+erfunc(calpi*rij)/rij)
                            end if
                         end if

                         if (rijsq.lt.rcutsq.or.lijall) then
                             fij = fij + fij_calculation(i, imolty, ii, j, jmolty, jj, rijsq, rij, ntii, ntjj, ntij)
                         end if
                         fxcmi = fxcmi + fij * rxuij
                         fycmi = fycmi + fij * ryuij
                         fzcmi = fzcmi + fij * rzuij
                      end do bead2
                   end do

                   ! calculate distance between c-o-m ---
                   rxuij = xcmi - xcm(j)
                   ryuij = ycmi - ycm(j)
                   rzuij = zcmi - zcm(j)
                   ! minimum image the pair separations ***
                   if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)

                   press = press +  fxcmi*rxuij + fycmi*ryuij + fzcmi*rzuij

                   ! for surface tension
                   ! this is correct for the coulombic part and for LJ.  Note sign difference!
                   pxx = pxx - fxcmi*rxuij
                   pyy = pyy - fycmi*ryuij
                   pzz = pzz - fzcmi*rzuij
                   pips(1,2) = pips(1,2) - rxuij*fycmi
                   pips(1,3) = pips(1,3) - rxuij*fzcmi
                   pips(2,1) = pips(2,1) - ryuij*fxcmi
                   pips(2,3) = pips(2,3) - ryuij*fzcmi
                   pips(3,1) = pips(3,1) - rzuij*fxcmi
                   pips(3,2) = pips(3,2) - rzuij*fycmi
                end if
             end do !< loop over j>i
          end if !< if lkdtree
       end if    !< if particle i is in this box
    end do       !< loop over particle i

    call mp_sum(press,1,groupid)
    call mp_sum(pxx,1,groupid)
    call mp_sum(pyy,1,groupid)
    call mp_sum(pzz,1,groupid)
    call mp_sum(pips,size(pips),groupid)
! ################################################################

    if (lewald) then
       ! Compute the reciprocal space contribution
       ! by using the thermodynamic definition
       call recippress(ibox,repress,rpxx,rpyy,rpzz,rpxy,rpyx,rpxz,rpzx,rpyz,rpzy)
       press = press - repress
       pxx = pxx + rpxx
       pyy = pyy + rpyy
       pzz = pzz + rpzz
       pips(1,2) = pips(1,2) + qqfact*rpxy
       pips(1,3) = pips(1,3) + qqfact*rpxz
       pips(2,1) = pips(2,1) + qqfact*rpyx
       pips(2,3) = pips(2,3) + qqfact*rpyz
       pips(3,1) = pips(3,1) + qqfact*rpzx
       pips(3,2) = pips(3,2) + qqfact*rpzy
    end if

    pips(1,1) = pxx
    pips(2,2) = pyy
    pips(3,3) = pzz

    ! Compute the Gaussian well contribution to the pressure
    pwell = 0.0E0_dp
    pwellips = 0.0E0_dp

    ! KM for MPI - comment:
    ! this could likely be parallelized in the future
    if (lmipsw) then
       do i = 1,nchain
          imolty = moltyp(i)
          if (lwell(imolty)) then
             rxui = xcm(i)
             ryui = ycm(i)
             rzui = zcm(i)
             do j = 1, nwell(imolty)*nunit(imolty)
                k = j - int(j/nunit(imolty))*nunit(imolty)
                if (k.eq.0) k = nunit(imolty)
                rxuij = rxui-rxwell(j,imolty)
                ryuij = ryui-rywell(j,imolty)
                rzuij = rzui-rzwell(j,imolty)
                call mimage (rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                rcm = rcut(ibox)+rcmu(i)
                rcmsq = rcm*rcm
                if (rijsq.lt.rcmsq) then
                   do ii = 1, nunit(imolty)
                      if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                      rxui = rxu(i,ii)
                      ryui = ryu(i,ii)
                      rzui = rzu(i,ii)
                      rxuij = rxui-rxwell(j,imolty)
                      ryuij = ryui-rywell(j,imolty)
                      rzuij = rzui-rzwell(j,imolty)
                      call mimage (rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                      fij = 2.0E0_dp*awell(ii,k,imolty)*bwell* exp(-bwell*rijsq)
                      pwell = pwell+fij*rijsq
                      pwellips(1,1) = pwellips(1,1)+fij*rxuij*rxuij
                      pwellips(2,2) = pwellips(2,2)+fij*ryuij*ryuij
                      pwellips(3,3) = pwellips(3,3)+fij*rzuij*rzuij
                      pwellips(1,2) = pwellips(1,2)+fij*rxuij*ryuij
                      pwellips(1,3) = pwellips(1,3)+fij*rxuij*rzuij
                      pwellips(2,3) = pwellips(2,3)+fij*ryuij*rzuij
                   end do
                end if
             end do
          end if
       end do
       pwellips(2,1) = pwellips(1,2)
       pwellips(3,1) = pwellips(1,3)
       pwellips(3,2) = pwellips(2,3)

       pipsw = -(press/3.0E0_dp)/vol
       pwellipsw = -(pwell/3.0E0_dp)/vol

       do i = 1, 3
          do j = 1, 3
             pips(i,j) = pips(i,j)/vol
             pwellips(i,j) = -pwellips(i,j)/vol
          end do
       end do

       if (lstagea) then
          press = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*press
       else if (lstageb) then
          press = etais*press+lambdais*pwell
       else if (lstagec) then
          press = (etais+(1.0E0_dp-etais)*lambdais)*press +(1.0E0_dp-lambdais)*pwell
       end if
    end if

    press = ((nchbox(ibox)+ghost_particles(ibox))/beta - press/3.0_dp)/vol
    surf = pzz - 0.5_dp*(pxx + pyy)
    ! divide by surface area and convert from K to put surf in mN/m
    surf = k_B*1.0E23_dp*surf/(2.0_dp*boxlx(ibox)*boxly(ibox))

    ! tail correction and impulsive force contribution to pressure
    if (ltailc.or.(.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3)) then
       ! add tail corrections for the Lennard-Jones energy
       ! Not adding tail correction for the ghost particles
       ! as they are ideal (no interaction) Neeraj.
       volsq = vol*vol
       do imolty=1, nmolty
          do jmolty=1, nmolty
             rhosq = ncmt(ibox,imolty)*ncmt(ibox,jmolty)/volsq
             press=press + corp(imolty,jmolty,rhosq,ibox)
          end do
       end do
    end if

    comp = press * vol * beta/(nchbox(ibox)+ghost_particles(ibox))
    press = k_B*1.0E27_dp * press
   return
  end subroutine pressure

!> \brief Impulsive force and tail corrections to pressure
!>
!> \warning Both corrections for tabulated potential and
!> Feuston-Garofalini, and tail corrections for MMFF94 are
!> NOT supported.
!> \note New potenial functional form needs to be added here,
!> and in U2 (energy calculation), coru (tail corrections to
!> energy), prop_pressure::pressure (force calculation)
  function corp(imolty,jmolty,rhosq,ibox)
    use const_math,only:twopi
    real::corp
    integer,intent(in)::imolty,jmolty,ibox
    real,intent(in)::rhosq
    real::rbcut,rbcut2,rbcut3,tmp,rci1,rci3,rci6,sigma2,epsilon2
    integer::ii,jj,ntii,ntjj,ntij

    rbcut=rcut(ibox)
    rbcut2=rbcut*rbcut
    rbcut3=rbcut2*rbcut

    corp = 0.0_dp
    do ii = 1, nunit(imolty)
       ntii = ntype(imolty,ii)
       do jj = 1, nunit(jmolty)
          ntjj = ntype(jmolty,jj)
          ntij=type_2body(ntii,ntjj)
          if (lgaro.or.lsami.or.lmuir) then
          else if (nonbond_type(ntij).eq.1) then
             ! LJ 12-6
             if (lexpand(imolty).and.lexpand(jmolty)) then
                sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))/2.0_dp
                epsilon2 = 4.0_dp*sqrt(epsilon_f(imolty,ii)*epsilon_f(jmolty,jj))
             else if (lexpand(imolty)) then
                sigma2 = (sigma_f(imolty,ii)+vvdW_b(2,ntjj))/2.0_dp
                epsilon2 = 4.0_dp*sqrt(epsilon_f(imolty,ii)*vvdW_b(1,ntjj))
             else if (lexpand(jmolty)) then
                sigma2 = (vvdW_b(2,ntii)+sigma_f(jmolty,jj))/2.0_dp
                epsilon2 = 4.0_dp*sqrt(vvdW_b(1,ntii)*epsilon_f(jmolty,jj))
             else
                sigma2 = vvdW(2,ntij)
                epsilon2 = vvdW(1,ntij)
             end if

             tmp  = sigma2**3
             rci3 = tmp/rbcut3
             rci6 = rci3*rci3
             if (ltailc) then ! tail corrections
                corp = corp + epsilon2*tmp*rci3*(rci6*4.0_dp/3.0_dp-2.0_dp)
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                corp = corp + epsilon2*tmp*rci3*(rci6-1.0_dp)
             end if
          else if (nonbond_type(ntij).eq.2) then
             ! Buckingham exp-6
             if (ltailc) then ! tail corrections
                if (vvdW(2,ntij).ne.0) then
                   tmp = vvdW(2,ntij)*rbcut
                   corp = corp - 2.0_dp*vvdW(3,ntij)/rbcut3 + (-6.0_dp+tmp*(6.0_dp+tmp*(tmp-3.0_dp)))*vvdW(1,ntij)*exp(tmp)/(vvdW(2,ntij)**3)
                end if
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                corp = corp + vvdW(1,ntij)*exp(vvdW(2,ntij)*rbcut)*rbcut3 - vvdW(3,ntij)/rbcut3
             end if
          else if (nonbond_type(ntij).eq.3) then
             ! Mie
             tmp=vvdW(2,ntij)/rbcut
             if (ltailc) then ! tail corrections
                rci1=tmp**(vvdW(3,ntij)-3)
                rci3=tmp**(vvdW(4,ntij)-3)
                corp = corp + vvdW(1,ntij)*(vvdW(2,ntij)**3)*(rci1*vvdW(3,ntij)/(vvdW(3,ntij)-3)-rci3*vvdW(4,ntij)/(vvdW(4,ntij)-3))
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                corp = corp + vvdW(1,ntij)*(tmp**vvdW(3,ntij)-tmp**vvdW(4,ntij))*rbcut3
             end if
          else if (nonbond_type(ntij).eq.4) then
             ! MMFF94 or DPD
             if (ltailc) then ! tail corrections for MMFF94 are not implemented yet
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                if (vvdW(2,ntij).ne.0) then
                   rci1 = rbcut/vvdW(2,ntij)
                   tmp = rci1**7
                   corp = corp + vvdW(1,ntij) * ((1.07_dp/(rci1+0.07_dp))**7) * (1.12_dp/(tmp+0.12_dp)-2.0_dp) * rbcut3
                end if
             end if
          else if (nonbond_type(ntij).eq.5) then
             ! LJ 9-6
             rci3=(vvdW(2,ntij)/rbcut)**3
             if (ltailc) then ! tail corrections
                corp = corp + vvdW(1,ntij)*(vvdW(2,ntij)**3)*rci3*(3.0_dp*rci3-6.0_dp)
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                corp = corp + vvdW(1,ntij)*rci3*rci3*(2.0_dp*rci3 - 3.0_dp)
             end if
          else if (nonbond_type(ntij).eq.6) then
             ! Generalized LJ
             if (ltailc) then ! tail corrections
                tmp  = 2.0_dp*vvdW(4,ntij)-3.0_dp
                rci1 = (vvdW(2,ntij)/rbcut)**tmp
                rci3 = (vvdW(2,ntij)/rbcut)**(vvdW(4,ntij)-3.0_dp)
                corp = corp + vvdW(1,ntij)*(vvdW(2,ntij)**3)*2.0_dp*vvdW(4,ntij)*(rci1/tmp-rci3/(vvdW(4,ntij)-3.0_dp))
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                tmp = (vvdW(2,ntij)/rbcut)**vvdW(4,ntij)
                corp = corp + vvdW(1,ntij)*tmp*(tmp-2.0_dp)
             end if
          else if (nonbond_type(ntij).eq.7) then
             ! LJ 12-6-8
             if (ltailc) then ! tail corrections
                corp = corp + 4.0_dp/3.0_dp*vvdW(1,ntij)/(rbcut3**3)-2.0_dp*vvdW(2,ntij)/rbcut3-8.0_dp/5.0_dp*vvdW(3,ntij)/rbcut3/rbcut2
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                tmp=rbcut2**3
                corp = corp + vvdW(1,ntij)/tmp/tmp-vvdW(2,ntij)/tmp-vvdW(3,ntij)/tmp/rbcut2
             end if
          else if (nonbond_type(ntij).eq.10) then
             ! LJ 12-6-8-10
             if (ltailc) then ! tail corrections
                corp = corp + (4.0_dp*vvdW(1,ntij))/(3.0_dp*rbcut3**3)-(2.0_dp*vvdW(2,ntij))/rbcut3-(8.0_dp*vvdW(3,ntij))/(5.0_dp*rbcut3*rbcut2)-(10.0_dp*vvdW(4,ntij))/(7.0_dp*rbcut3*rbcut3*rbcut)
             else if (.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3) then ! impulsive force corrections
                tmp=rbcut2**3
                corp = corp + vvdW(1,ntij)/(tmp*tmp)-vvdW(2,ntij)/tmp-vvdW(3,ntij)/(tmp*rbcut2)-vvdW(4,ntij)/(tmp*rbcut2*rbcut2)
             end if
          else if (ALL(nonbond_type(ntij).ne.(/-1,0,8,9/))) then
             call err_exit(__FILE__,__LINE__,'corp: undefined nonbond type',myid+1)
          end if
       end do
    end do

    corp=twopi/3.0_dp*rhosq*corp

    return
  end function corp

  ! calculate fij in the pressure calculation, given the molecule number i and j
  ! the molecule type imolty and jmolty
  ! the bead number ii and jj
  ! the distance squared between two beads, rijsq
  ! and the interaction type ntii, ntjj and ntij
  function fij_calculation(i, imolty, ii, j, jmolty, jj, rijsq, rij, ntii, ntjj, ntij) result(fij)
      real :: fij
      integer, intent(in) :: i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij
      real, intent(in) :: rijsq, rij
      real::rs1,rs2,rs4,rs6,rs7,rs8,rs14,sr1,sr2,sr3,sr6,sr7,sigma2,epsilon2,qave,tmp

      fij = 0.0E0_dp
      !> \todo when lsami, lmuir, lgaro is .true.
      if (lgaro.or.lsami.or.lmuir) then
      else if (nonbond_type(ntij).eq.1) then
          ! LJ 12-6
          if (lexpand(imolty).and.lexpand(jmolty)) then
              sigma2=(sigma_f(imolty,ii)+sigma_f(jmolty,jj))/2.0_dp
              sr2=sigma2*sigma2/rijsq
              epsilon2=4.0_dp*sqrt(epsilon_f(imolty,ii)*epsilon_f(jmolty,jj))
          else if (lexpand(imolty)) then
              sigma2=(sigma_f(imolty,ii)+vvdW_b(2,ntjj))/2.0_dp
              sr2=sigma2*sigma2/rijsq
              epsilon2=4.0_dp*sqrt(epsilon_f(imolty,ii)*vvdW_b(1,ntjj))
          else if (lexpand(jmolty)) then
              sigma2=(vvdW_b(2,ntii)+sigma_f(jmolty,jj))/2.0_dp
              sr2=sigma2*sigma2/rijsq
              epsilon2=4.0_dp*sqrt(vvdW_b(1,ntii)*epsilon_f(jmolty,jj))
          else
              sr2=vvdW(3,ntij)/rijsq
              epsilon2=vvdW(1,ntij)
          end if

          if (lfepsi) then
              sr6 = (1.0_dp/rijsq)**3
              if ((.not.lqchg(ntii)).and.(.not.lqchg(ntjj))) then
                  if (nunit(imolty).eq.4) then
                      !> \bug TIP-4P structure (temperary use?)
                      qave = (qqu(i,4)+qqu(j,4))/2.0_dp
                  else
                      qave = (qqu(i,4)+qqu(i,5)+qqu(j,4)+qqu(j,5))*0.85_dp
                  end if
              else
                  qave=(qqu(i,ii)+qqu(j,jj))/2.0_dp
              end if
              fij = fij + 12.0_dp*epsilon2*sr6*(-sr6*(aslope*(qave-a0)*(qave-a0)+ashift) &
                    +0.5_dp*(bslope*(qave-b0)*(qave-b0)+ bshift))/rijsq
          else
              sr6 = sr2**3
              fij = fij - 12.0_dp*epsilon2*sr6*(sr6-0.5_dp)/rijsq
          end if
      else if (nonbond_type(ntij).eq.2) then
          ! Buckingham exp-6
          rs1=vvdW(2,ntij)*rij
          fij = fij + (vvdW(1,ntij)*rs1*exp(rs1)+6.0_dp*vvdW(3,ntij)/(rijsq**3))/rijsq
      else if (nonbond_type(ntij).eq.3) then
          ! Mie
          sr1 = vvdW(2,ntij) / rij
          fij = fij - vvdW(1,ntij)*(vvdW(3,ntij)*sr1**vvdW(3,ntij)-vvdW(4,ntij)*sr1**vvdW(4,ntij))/rijsq
      else if (nonbond_type(ntij).eq.4) then
          ! MMFF94
          if (vvdW(2,ntij).ne.0) then
              rs1 = rij/vvdW(2,ntij)
              rs2 = rs1*rs1
              rs6 = rs2**3
              rs7 = rs1*rs6
              sr1 = 1.07_dp/(rs1+0.07_dp)
              sr7 = sr1**7.0_dp
              fij = fij - 7.0_dp*vvdW(1,ntij)*rs1*sr7*( sr1*(1.12_dp/(rs7+0.12_dp)-2.0_dp)/1.07_dp &
                    + 1.12_dp*rs6/((rs7+0.12_dp)*(rs7+0.12_dp)) )/rijsq
          end if
      else if (nonbond_type(ntij).eq.5) then
          ! LJ 9-6
          sr3=(vvdW(2,ntij)/rij)**3
          sr6=sr3*sr3
          fij = fij - 18.0_dp*vvdW(1,ntij)*sr6*(sr3-1.0_dp)/rijsq
      else if (nonbond_type(ntij).eq.6) then
          ! Generalized LJ
          sr1 = vvdW(2,ntij)/rij
          if (rij.le.vvdW(2,ntij)) then
              tmp = sr1**(vvdW(3,ntij)/2.0_dp)
              fij = fij - vvdW(1,ntij)*vvdW(3,ntij)*tmp*(tmp-1.0_dp)/rijsq
          else
              tmp = sr1**vvdW(4,ntij)
              fij = fij - vvdW(1,ntij)*2.0_dp*vvdW(4,ntij)*tmp*(tmp-1.0_dp)/rijsq
          end if
      else if (nonbond_type(ntij).eq.7) then
          ! LJ 12-6-8
          rs4=rijsq*rijsq
          rs8=rs4*rs4
          tmp=rs8*rs4*rijsq !^14
          fij = fij - 12.0_dp*vvdW(1,ntij)/tmp + 6.0_dp*vvdW(2,ntij)/rs8 + 8.0_dp*vvdW(3,ntij)/rs8/rijsq
      else if (nonbond_type(ntij).eq.8) then
          ! DPD
          ! In general the force between two DPD particles is
          ! given by
          ! F(r_ij) = -\del_{ij} U(r_{ij})
          !         = -\del_{ij} (a_{ij}/2 (1-r_{ij}/rmin)^2
          ! applying the chain rule this gives:
          !         = -(a_{ij}/2)*(2(1-r_{ij}/rmin))*(-1/rmin)*\hat{r}_{ij}
          ! simplifying to:
          !         =(a_{ij}/rmin)*(1-r_{ij}/rmin)*\hat{r}_{ij}
          ! the form here looks totally bizarrre due to the
          ! way the code handles the different pieces of
          ! information.
          if(rij<=vvdW(2,ntij).and.rij>0.0_dp) then ! if rij<rmin
              fij=fij-2.0_dp*vvdW(1,ntij)*(1.0_dp-rij/vvdW(2,ntij))/(rij*vvdW(2,ntij))
          end if
      else if (nonbond_type(ntij).eq.10) then
          ! LJ 12-6-8-10
          rs4=rijsq*rijsq
          rs8=rs4*rs4
          rs14=rs8*rs4*rijsq !^14
          fij = fij - (12.0_dp*vvdW(1,ntij))/rs14 + (6.0_dp*vvdW(2,ntij))/rs8 + (8.0_dp*vvdW(3,ntij))/(rs8*rijsq) + (10.0_dp*vvdW(4,ntij))/(rs8*rs4)
      else if (ALL(nonbond_type(ntij).ne.(/-1,0,9/))) then
          call err_exit(__FILE__,__LINE__,'pressure: undefined nonbond type',myid+1)
      end if

      return

  end function fij_calculation

  ! Pressure calculation when bead-kdtree is used
  ! ibox, i: box and molecule of interest
  ! press: accumulated pressure
  ! pxx, pyy, pzz: accumulated surface tension?
  subroutine press_calc_kdtree(ibox, i, press, pxx, pyy, pzz)
     use const_math, only : sqrtpi
     integer, intent(in) :: ibox, i
     real :: rbcut, press, pxx, pyy, pzz

     type(tree), pointer :: kd_tree
     real :: coord(3), fxcm(nchain), fycm(nchain), fzcm(nchain)
     logical :: ovrlap, lqimol, lqjmol, lcoulo(numax,numax)
     real, allocatable :: inter_list(:,:)
     integer :: iInter, inter_list_dim, dist_calc_num, dist_calc_num_temp
     real :: rij, rijsq, rcutsq, rxuij, ryuij, rzuij, xcmi, ycmi, zcmi, fij, calpi, calpisq
     integer :: ii, iii, j, jj, jjj, imolty, jmolty, ntii, ntjj, ntij

     kd_tree => mol_tree(ibox)%tree
     rbcut = rcut(ibox)
     rcutsq = rbcut * rbcut
     imolty = moltyp(i)
     lqimol = lelect(imolty)
     xcmi = xcm(i)
     ycmi = ycm(i)
     zcmi = zcm(i)
     calpi=calp(ibox)
     calpisq=calpi*calpi

     fxcm = 0.0E0_dp !< these are 1d arrays
     fycm = 0.0E0_dp
     fzcm = 0.0E0_dp

     ! find all bead pairs and accumulate forces
     do ii = 1, nunit(imolty)
         coord(1) = rxu(i, ii)
         coord(2) = ryu(i, ii)
         coord(3) = rzu(i, ii)
         ntii = ntype(imolty, ii)

         ! find beads that are within the cutoff from bead ii
         call range_search(kd_tree, coord, nmax*numax, 0.0, rbcut, i, ii, .true., ovrlap, inter_list_dim,&
                            inter_list, dist_calc_num_temp, .true.)

         ! loop over all beads found
         bead: do iInter = 1, inter_list_dim
             rijsq = inter_list(1, iInter)
             j = inter_list(2, iInter)
             jj = inter_list(3, iInter)
             jmolty = moltyp(j)

             ! check exclusion table
             if (lexclu(imolty,ii,jmolty,jj)) cycle bead

             rxuij = inter_list(4, iInter)
             ryuij = inter_list(5, iInter)
             rzuij = inter_list(6, iInter)
             ntjj = ntype(jmolty, jj)
             ntij = type_2body(ntii, ntjj)
             lqjmol = lelect(jmolty)
             rij = sqrt(rijsq)
             fij = 0.0E0_dp

             ! calculate the charge interactions
             if (lqimol .and. lqjmol .and. lqchg(ntii) .and. lqchg(ntjj)) then
                 if (.not. lewald) then
                      if (.not. lchgall) then
                         iii = leaderq(imolty, ii)
                         jjj = leaderq(jmolty, jj)
                         if ((iii .eq. ii) .and. (jjj .eq. jj)) then
                             ! set up the charge-interaction table
                             if (rijsq.lt.rcutsq) then
                                lcoulo(iii,jjj) = .true.
                             else
                                lcoulo(iii,jjj) = .false.
                             end if
                          end if
                      end if
                      !> \todo tabulated potential
                      if (lchgall.or.lcoulo(iii,jjj) ) then
                          fij = -qqfact*qqu(i,ii)*qqu(j,jj)/rijsq/rij
                      end if
                 else
                      fij = -qqfact*qqu(i,ii)*qqu(j,jj)/rijsq&
                            *(2.0_dp*calpi*exp(-calpisq*rijsq)/sqrtpi+erfunc(calpi*rij)/rij)
                 end if
             end if

             ! LJ part of the fij
             fij = fij + fij_calculation(i, imolty, ii, j, jmolty, jj, rijsq, rij, ntii, ntjj, ntij)

             ! Accumulate the forces in f*cm arrays
             fxcm(j) = fxcm(j) + rxuij * fij
             fycm(j) = fycm(j) + ryuij * fij
             fzcm(j) = fzcm(j) + rzuij * fij

          end do bead
    end do

    ! loop over jchain and calculate pressure accumulants
    do j = i + 1, nchain
         if ((nboxi(j) .eq. ibox) .and. &
              ((abs(fxcm(j)) .gt. 0.0) .or. (abs(fycm(j)) .gt. 0.0) .or. (abs(fzcm(j)) .gt. 0.0))) then
             rxuij = xcmi - xcm(j)
             ryuij = ycmi - ycm(j)
             rzuij = zcmi - zcm(j)
             if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)

             press = press + fxcm(j)*rxuij + fycm(j)*ryuij + fzcm(j)*rzuij

             ! for surface tension
             ! this is correct for the coulombic part and for LJ.  Note sign difference!
             pxx = pxx - fxcm(j) * rxuij
             pyy = pyy - fycm(j) * ryuij
             pzz = pzz - fzcm(j) * rzuij
             pips(1,2) = pips(1,2) - rxuij * fycm(j)
             pips(1,3) = pips(1,3) - rxuij * fzcm(j)
             pips(2,1) = pips(2,1) - ryuij * fxcm(j)
             pips(2,3) = pips(2,3) - ryuij * fzcm(j)
             pips(3,1) = pips(3,1) - rzuij * fxcm(j)
             pips(3,2) = pips(3,2) - rzuij * fycm(j)
         end if
    end do

  end subroutine press_calc_kdtree

end MODULE prop_pressure
