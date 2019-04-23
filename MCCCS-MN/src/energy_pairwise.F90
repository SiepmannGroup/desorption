MODULE energy_pairwise
  use var_type,only:dp,default_string_length
  use const_math,only:onepi,twopi,sqrtpi,degrad
  use const_phys,only:qqfact
  use util_math,only:erfunc
  use util_runtime,only:err_exit
  use util_string,only:uppercase
  use util_files,only:readLine
  use util_mp,only:mp_sum,mp_lor,mp_allgather
  use util_kdtree,only:update_tree_height,range_search
  use sim_system
  use sim_cell
  use energy_kspace,only:calp,sself,correct
  use energy_intramolecular,only:lininter_bend
  use energy_external,only:U_ext
  use energy_sami,only:ljsami,ljmuir
  use energy_garofalini
  use energy_3body,only:hasThreeBody
  use energy_4body,only:hasFourBody
  implicit none
  private
  save
  public::sumup,energy,boltz,coru,read_ff,init_ff,U2,type_2body,vdW_nParameter,nonbond_type,energy_inter_com_kd_tree

  real,parameter::overlapValue=1.0E+20_dp,a15(2)=(/4.0E7_dp,7.5E7_dp/) !< 1-5 correction term for unprotected hydrogen-oxygen interaction; 1 for ether oxygens, 2 for alcohol oxygens
  !< OLD VALUES: a15(2)=/17.0_dp**6,16.0_dp**6/)
  integer,parameter::vdW_nParameter(-1:10)=(/0,0,2,3,4,2,2,4,3,2,3,4/)
  integer,allocatable::atom_type(:),nonbond_type(:) !< type -1: tabulated potential
  !< type 1: Lennard-Jones 12-6, U(r) = 4*epsilon*[(sigma/r)^12-(sigma/r)^6]
  !< vvdW_b_1 = epsilon, vvdW_b_2 = sigma
  !< vvdW_1 = 4*epsilon, vvdW_2 = sigma, vvdW_3 = sigma^2
  !< type 2: Buckingham exponential-6, U(r) = A*exp(-B*r) - C/r^6
  !< vvdW_b_1 = A, vvdW_b_2 = B, vvdW_b_3 = C
  !< vvdW_1 = A, vvdW_2 = -B, vvdW_3 = C
  !< type 3: Mie, U(r) = C*epsilon*[(sigma/r)^n0-(sigma/r)^n1], where C = n0/(n0-n1) * (n0/n1)^[n1/(n0-n1)]
  !< vvdW_b_1 = epsilon, vvdW_b_2 = sigma, vvdW_b_3 = n0, vvdW_b_4 = n1
  !< vvdW_1 = C*epsilon, vvdW_2 = sigma, vvdW_3 = n0, vvdW_4 = n1
  !< type 4: MMFF94 (Merck Molecular Force Field: Thomas A. Halgren, J Am Chem Soc 1992,114:7827-7843) buffered 14-7, U(r) = 4*epsilon*{1.07/[(r/r0)+0.07]}^7 * {1.12/[(r/r0)^7+0.12]-2}
  !< vvdW_b_1 = epsilon, vvdW_b_2 = r0
  !< vvdW_1 = 4*epsilon, vvdW_2 = r0, vvdW_3 = r0^2
  !< type 5: Lennard-Jones 9-6, U(r) = 4*epsilon*[2*(r0/r)^9-3*(r0/r)^6],
  !< vvdW_b_1 = epsilon, vvdW_b_2 = r0
  !< vvdW_1 = 4*epsilon, vvdW_2 = r0
  !< type 6: Generalized Lennard-Jones (Ref: J Chem Phys 2004,120:4994),
  !< U(r) = 4*epsilon*[(r0/r)^n0-2*(r0/r)^(n0/2)] if r<=r0
  !<      = 4*epsilon*[(r0/r)^(2*n1)-2*(r0/r)^n1] if r>r0
  !< vvdW_b_1 = epsilon, vvdW_b_2 = r_0, vvdW_b_3 = n0, vvdW_b_4 = n1
  !< vvdW_1 = 4*epsilon, vvdW_2 = r_0, vvdW_3 = n0, vvdW_4 = n1
  !< type 7: Lennard-Jones 12-6-8, U(r) = A/r^12 - B/r^6 - C/r^8
  !< vvdW_b_1 = A, vvdW_b_2 = B, vvdW_b_3 = C
  !< vvdW_1 = A, vvdW_2 = B, vvdW_3 = C
  !< type 8: DPD potential
  !< U(r) = a/2 (1-r/rmin)^2 if rij <= rmin
  !<      = 0.0                 if rij > rmin
  !< vvdW_b_1 = a, vvdW_b_2 = rmin
  !< vvdW_1 = a, vvdW_2 = rmin, vvdW_3 = rmin^2
  !< type 9: Hard-core square-well,
  !< U(r) = +inf     if r < sigma
  !<      = -epsilon if sigma <= r < lambda*sigma
  !<      = 0        if r >= lambda*sigma
  !< vvdW_b_1 = epsilon, vvdW_b_2 = sigma, vvdW_b_3 = lambda
  !< vvdW_1 = epsilon, vvdW_2 = sigma, vvdW_3 = lambda*sigma
  !< type 10: Lennard-Jones 12-6-8-10, U(r) = A/r^12 - C/r^6 - D/r^8 - E/r^10
  !< vvdW_b_1 = A, vvdW_b_2 = C, vvdW_b_3 = D, vvdW_b_4 = E

  integer,allocatable::vdWsplits(:,:),electsplits(:,:)
  real,allocatable::rvdW(:,:,:),tabvdW(:,:,:),relect(:,:,:),tabelect(:,:,:)
  integer::ntabvdW,ntabelect
contains
!*****************************************************************
!> \brief Calculates the total potential energy for a configuration.
!>
!> \param ovrlap logical, true for substantial atom overlap
!> \param v* energies
!> \param ibox box number
!> \param lvol true if called from moves_volume.F90, no output of summary infomation
!******************************************************************
  subroutine sumup(ovrlap,v,ibox,lvol)
    use sim_particle,only:init_neighbor_list,lnn,add_neighbor_list
    use energy_kspace,only:recipsum
    use energy_intramolecular,only:U_bonded
    use energy_3body,only:U3System
    use energy_4body,only:U4System

    real::v(nEnergy),vrecipsum,vwell,my_velect
    logical::ovrlap,lvol
    logical::lexplt,lqimol,lqjmol,lcoulo(numax,numax),lij2,liji,lqchgi
    integer::i,imolty,ii,j,jmolty,jj,ntii,ntjj,ntij,iunit,ibox,nmcount,ntj,k,mmm,jcell(nmax)
    integer::mole_start,mole_end
    real::rcutsq,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij,rijsq,rho,rij,rbcut,calpi,qqii
    ! real::vtemp
    real::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,vol
    ! Neeraj & RP for MPI
    real::sum_vvib,sum_vbend,sum_vtg
    ! Kdtree
    real, allocatable :: lij_list(:,:)
    integer :: i_run
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start SUMUP in ',myid,' for box ', ibox
#endif

    if ( lpbc ) call setpbc(ibox)

    rbcut = rcut(ibox)
    rcutsq = rbcut*rbcut
    calpi=calp(ibox)
    rminsq = rmin * rmin

    ! KM for MPI
    my_velect = 0.0E0_dp

    ovrlap = .false.
    v = 0.0E0_dp
    vwell = 0.0E0_dp

    ! check the molecule count ***
    nmcount = 0
    do i = 1, nchain
       if ( nboxi(i) .eq. ibox ) then
          nmcount=nmcount+1
       end if
    end do
    if ( nmcount .ne. nchbox(ibox) ) then
       call err_exit(__FILE__,__LINE__,'SUMUP: nmcount ne nchbox'//integer_to_string(nmcount)//integer_to_string(nchbox(ibox))&
        ,myid+1)
    end if

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
    if (lneigh.or.lneighbor.or.lgaro) call init_neighbor_list(ibox)

    ! loop over all chains i
    if (.not.lideal(ibox)) then
        if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
            ! computing inter-molecular interactions using kd-tree
            call energy_inter_kd_tree_sumup(ibox, lvol, v, ovrlap)
            if (ovrlap) return
        else
           ! MPI
           do i = 1, nchain - 1
              ! check if i is in relevant box ###
              if ( nboxi(i) .eq. ibox ) then
                 imolty = moltyp(i)
                 lqimol = lelect(imolty)

                 if ( lcutcm .and. lvol ) then
                    xcmi = xcm(i)
                    ycmi = ycm(i)
                    zcmi = zcm(i)
                    rcmi = rcmu(i)
                 else
                    lij2 = .true.
                 end if

                 if (licell.and.(ibox.eq.boxlink)) then
                    ii=1
                    rxui = xcm(i)
                    ryui = ycm(i)
                    rzui = zcm(i)
                    call get_cell_neighbors(rxui,ryui,rzui,ibox,jcell,mole_end)
                    mole_start = 1
                 else
                    mole_start = i + 1 + myid
                    mole_end = nchain
                 end if

                 if ( nugrow(imolty) .eq. nunit(imolty) ) then
                    lexplt = .false.
                 else
                    lexplt = .true.
                 end if

                 if (lcutcm .and. lkdtree .and. lkdtree_box(ibox)) then
                    call energy_inter_com_kd_tree(ibox, i, xcm(i), ycm(i), zcm(i), rbcut, rmin &
                        , .false., lij_list, .true., ovrlap)
                    if (ovrlap) goto 199
                 end if

                 ! loop over all chains j with j>i
                 molecule2: do k = mole_start, mole_end, numprocs

                    ! if cell list is used, skip the COM distance calculation part
                    if (licell.and.(ibox.eq.boxlink)) then
                        j = jcell(k)
                        if (j .le. i) then
                            cycle molecule2
                        else
                            jmolty = moltyp(j)
                            lij2 = .true.
                            goto 1024
                        end if
                    else
                        j = k
                    end if

                    ! check for simulation box ###
                    if ( nboxi(j) .eq. ibox ) then
                       jmolty = moltyp(j)
                       lqjmol = lelect(jmolty)

                       if (lcutcm .and. lkdtree .and. lkdtree_box(ibox)) then
                           if (lij_list(1, j) .lt. 0.5d0) then
                               cycle molecule2
                           else
                               lij2 = .true.
                               goto 1024
                           end if
                       end if

                       if (lcutcm .and. lvol ) then
                          ! check if ctrmas within rcmsq
                          rxuij = xcmi - xcm(j)
                          ryuij = ycmi - ycm(j)
                          rzuij = zcmi - zcm(j)
                          ! minimum image the ctrmas pair separations
                          if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)

                          rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                          rcm = rbcut + rcmi + rcmu(j)
                          rcmsq = rcm*rcm

                          ! if ( lneighbor .and. rcmsq .lt. rbsmax**2 .and. rcmsq .gt. rbsmin**2 ) then
                          ! neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
                          ! neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
                          ! neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
                          ! neighbor(neigh_cnt(j,imolty),j,imolty)=i
                          ! end if

                          if ( rijsq .gt. rcmsq ) then
                             if ( lqimol .and. lqjmol .and. lchgall ) then
                                lij2 = .false.
                             else
                                cycle molecule2
                             end if
                          else
                             lij2 = .true.
                          end if
                       end if

1024                   continue

                       do ii = 1,nunit(imolty)
                          ntii = ntype(imolty,ii)
                          liji = lij(ntii)
                          lqchgi = lqchg(ntii)
                          rxui = rxu(i,ii)
                          ryui = ryu(i,ii)
                          rzui = rzu(i,ii)

                          ! loop over all beads jj of chain j
                          bead2: do jj = 1, nunit(jmolty)
                             ! check exclusion table
                             if (lexclu(imolty,ii,jmolty,jj)) cycle bead2

                             ntjj = ntype(jmolty,jj)
                             if ( lij2 ) then
                                if ((.not.(liji.and.lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj)))) cycle bead2
                             else
                                if (.not.(lqchgi.and.lqchg(ntjj))) cycle bead2
                             end if

                             ntij=type_2body(ntii,ntjj)

                             if (lexpee) rminsq=rminee(ntij)*rminee(ntij)

                             rxuij = rxui - rxu(j,jj)
                             ryuij = ryui - ryu(j,jj)
                             rzuij = rzui - rzu(j,jj)
                             ! minimum image the pair separations ***
                             if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)

                             rijsq=(rxuij*rxuij)+(ryuij*ryuij)+(rzuij*rzuij)
                             if (rijsq.lt.rminsq .and. .not.(lexpand(imolty).or.lexpand(jmolty))) then
                                if ( .not. lvol .and.myid.eq.rootid) then
                                   write(io_output,*) 'overlap inter'
                                   write(io_output,*)'rijsq rminsq',rijsq,rminsq
                                   write(io_output,*) 'i ii', i, ii
                                   write(io_output,*) 'i-pos', rxui,ryui,rzui
                                   write(io_output,*) 'j jj', j, jj
                                   write(io_output,*) 'j-pos',  rxu(j,jj),ryu(j,jj),rzu(j,jj)
                                end if
                                ovrlap = .true.
                                ! RP added for MPI to compensate ovrlap
                                ! return
                                goto 199
                             else if (rijsq.lt.rcutsq .or. lijall) then
                                v(ivInterLJ)=v(ivInterLJ)+U2(rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                             end if

                             v(ivElect)=v(ivElect)+Q2(rijsq,rcutsq,i,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo,ibox)

                             if (lneigh.and.rijsq.le.rcutnn(ibox)**2) then
                                lnn(i,j) = .true.
                                lnn(j,i) = .true.
                             end if

!cc  KM for MPI
!cc  all processors need to know neighbor information
!cc  lneighbor and lgaro will not work in parallel
!cc  calculation of neighbors assumes everything is sequential
                             !> \bug Wouldn't testing for center-of-mass distance be a better criteria?
                             if (lneighbor.and.ii.eq.1.and.jj.eq.1.and.rijsq.lt.rbsmax**2.and.rijsq.gt.rbsmin**2) then
                                call add_neighbor_list(i,j,sqrt(rijsq),rxuij,ryuij,rzuij)
                             else if (lgaro) then
                                if (isNeighbor(rijsq,ntii,ntjj)) then
                                   call add_neighbor_list(i,j,sqrt(rijsq),rxuij,ryuij,rzuij)
                                end if
                             end if
                          end do bead2
                       end do !do ii = 1,nunit(imolty)
                    end if
                 end do molecule2
              end if
           end do !do i = 1, nchain - 1
           ! Returning from ovrlap--------------
199        continue

! KM don't check overlap until allreduce is finished
! if(ovrlap .eq. .true.)then
! write(io_output,*)'521: in sumup ovrlap=',ovrlap,'myid=',myid
! end if
! -----------------------------------------
           call mp_lor(ovrlap,1,groupid)
           if(ovrlap)then
               ! write(io_output,*)'530: in sumup ovrlap=',ovrlap,'myid=',myid
               return
           end if

           call mp_sum(v(ivInterLJ),1,groupid)
           call mp_sum(v(ivElect),1,groupid)

       end if ! if lkdtree

! KEA garofalini 3 body potential
       if (lgaro.and..not.lideal(ibox)) then
          call triad()
          v(iv3body)=vthreebody()
       end if

       if (hasThreeBody.and..not.lideal(ibox)) v(ivInterLJ)=v(ivInterLJ)+U3System(ibox)
       if (hasFourBody.and..not.lideal(ibox)) v(ivInterLJ)=v(ivInterLJ)+U4System(ibox)

       if (ltailc) then
          ! add tail corrections for the Lennard-Jones energy
          if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
             vol = cell_vol(ibox)
          else
             vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
          end if
          do jmolty = 1, nmolty
             rho = ncmt(ibox,jmolty) / vol
             do imolty = 1, nmolty
                v(ivTail) = v(ivTail) +  ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
             end do
          end do
          v(ivInterLJ) = v(ivInterLJ) + v(ivTail)
       end if
    end if
!$$$c      write(io_output,*)
!$$$c      write(io_output,*) '+++++++'
!$$$c      vtemp = v(ivElect)
!$$$c      write(io_output,*) 'direct space part:',v(ivElect)*qqfact

    if ( ldielect ) then
       call dipole(ibox,0)
    end if

    if (lewald.and..not.lideal(ibox)) then
       call recipsum(ibox,vrecipsum)
       ! update self terms and correction terms
       sself = 0.0E0_dp
       correct = 0.0E0_dp
       ! combine to reduce numerical error
       ! vsc = 0.0E0_dp

       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             do ii = 1,nunit(imolty)
                sself = sself + qqu(i,ii)*qqu(i,ii)
                ! 1.772.. is the square root of pi
                ! vsc = vsc - qqu(i,ii)*qqu(i,ii)*calpi/1.772453851E0_dp
                do jj = ii+1,nunit(imolty)
                   rxuij = rxu(i,ii) - rxu(i,jj)
                   ryuij = ryu(i,ii) - ryu(i,jj)
                   rzuij = rzu(i,ii) - rzu(i,jj)
                   ! lpbc is not called here because it's intra-chain interaction
                   rij = sqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)

                   ! correct should only be calculated if ii and jj should NOT interact,
                   ! so only calculating it if lqinclu is false
                   ! this part is 1,2 and 1,3
                   if (.not. lqinclu(imolty,ii,jj)) then
                      correct = correct + qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                      ! vsc = vsc + qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                   else
                      correct=correct+(1.0E0_dp - qscale2(imolty,ii,jj))*qqu(i,ii)*qqu(i,jj)* (erfunc(calpi*rij)-1.0E0_dp)/rij
                      ! vsc = vsc + (1.0E0_dp - qscale2(imolty,ii,jj))*qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                   end if
                end do
             end do
          end do
       end do

       call mp_sum(correct,1,groupid)
       call mp_sum(sself,1,groupid)

! vdipole = (dipolex*dipolex+dipoley*dipoley+dipolez*dipolez)*(2.0E0_dp*onepi)/(3.0E0_dp*boxlx(ibox)**3.0E0_dp)
! write(io_output,*) dipolex,dipoley,dipolez
       sself = -sself*calpi/sqrtpi
       v(ivElect) = v(ivElect) + sself + correct + vrecipsum
! v(ivElect) = v(ivElect) + vsc + vrecipsum
    end if

!$$$c at this point velect contains all intermolecular charge interactions,
!$$$c plus the ewald self term and intramolecular corrections
!$$$
! write(io_output,*)
! write(io_output,*) '== After Inter === velect is:',v(ivElect)*qqfact
!$$$
!$$$       vtemp = v(ivElect)

! ################################################################

! have to recalculate ewald terms if volume changes
    if (.not.lvol.or.(lvol.and.ANY(lelect(1:nmolty)))) then
! *******************************
! INTRACHAIN INTERACTIONS ***
! *******************************
       ! write(io_output,*) 'starting intrachain'
       ! loop over all chains i
       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             do ii = 1, nunit(imolty)-1
! write(io_output,*) 'ntype(imolty,ii),ii',ntype(imolty,ii),ii
                ntii = ntype(imolty,ii)
                rxui = rxu(i,ii)
                ryui = ryu(i,ii)
                rzui = rzu(i,ii)
                do jj = ii+1, nunit(imolty)
                   if (linclu(imolty,ii,jj).or.lqinclu(imolty,ii,jj)) then
                      ntjj = ntype(imolty,jj)

                      ntij=type_2body(ntii,ntjj)

                      if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
! KNS 11-23-15 insert rxui, ryui rzui from above
                      rxuij = rxui - rxu(i,jj)
                      ryuij = ryui - ryu(i,jj)
                      rzuij = rzui - rzu(i,jj)
                      ! lpbc is not called here because it's intra-chain interaction
                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

                      if (linclu(imolty,ii,jj)) then
                         if (rijsq.lt.rminsq.and..not.lexpand(imolty)) then
                            if ( .not. lvol ) then
                               write(io_output,*) 'overlap intra'
                               write(io_output,*) 'rijsq rminsq', rijsq,  rminsq
                               write(io_output,*) 'i ii', i, ii
                               write(io_output,*) 'i-pos', rxui,ryui,rzui
                               write(io_output,*) 'jj', jj
                               write(io_output,*) 'j-pos' ,rxu(i,jj),ryu(i,jj),rzu(i,jj)
                            end if
                            ovrlap = .true.
                            ! RP added for MPI to compensate ovrlap
                            ! return
                            goto 299
                            ! -------------------------------
                         else
                            ! there's no cutoff when calculating the intra-chain interaction
                            !> \bug skip intra if it is bending 1-3 and using a table?
                            if (L_bend_table) then
                               do mmm=1,inben(imolty,ii)
                                  if (ijben3(imolty,ii,mmm).eq.jj)then
                                     v(ivIntraLJ) = v(ivIntraLJ) + lininter_bend(sqrt(rijsq),itben(imolty,ii,mmm))
                                     goto 94
                                  end if
                               end do
                            end if

                            v(ivIntraLJ)=v(ivIntraLJ)+U2(rijsq,i,imolty,ii,ntii,i,imolty,jj,ntjj,ntij)
                         end if
                      end if !if (linclu(imolty,ii,jj))

94                    if (lqinclu(imolty,ii,jj)) then
                         ! calculate intramolecular charge interaction
                         my_velect=my_velect+qscale2(imolty,ii,jj)*Q2(rijsq,rcutsq,i,imolty,ii,ntii,lqchg(ntii),i,imolty,jj&
                          ,ntjj,calpi,lcoulo,ibox)
                      end if
                   end if
                end do
             end do
          end do
       end do

! RP added for MPI----- Returning from ovrlap--------------
299    continue

       call mp_lor(ovrlap,1,groupid)
       if(ovrlap)then
          ! write(io_output,*)'941: in sumup ovrlap=',ovrlap,'myid=',myid
          return
       end if
! -----------------------------------------
       call mp_sum(v(ivIntraLJ),1,groupid)
       call mp_sum(my_velect,1,groupid)
       v(ivElect) = v(ivElect) + my_velect

!$$$       vtemp = v(ivElect) - vtemp
!$$$
!$$$c       write(io_output,*) '== Intra Velect ===',vtemp*qqfact
! write(io_output,*) '== After Intra  === velect is:',v(ivElect)*qqfact
! write(io_output,*) 'vintra ', v(ivIntraLJ)
! write(io_output,*) 'vinter ', v(ivInterLJ)
! write(io_output,*) 'test', v(ivIntraLJ)
! ################################################################

! *************************************
! INTRACHAIN FLUCQ INTERACTIONS ***
! *************************************
!c RP added for MPI
!> \bug removed at some point?
       do i = 1,nchain
! do i = myid+1,nchain,numprocs
! calculate intramolecular flucq energy for chain i
          if( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             if ( lelect(imolty) ) then

                if ( lflucq(imolty) ) then
                   iunit = nunit(imolty)
                   do ii = 1, iunit

                      ntii = ntype(imolty,ii)
                      qqii = qqu(i,ii)
                      do jj = ii, iunit

                         if ( ii .eq. jj) then
                            v(ivFlucq) = v(ivFlucq) + xiq(ntii)*qqii + jayself(ntii)*qqii*qqii
                         else
                            ntjj = ntype(imolty,jj)
                            ntij = type_2body(ntii,ntjj)

                            v(ivFlucq) = v(ivFlucq)  + jayq(ntij)*qqii*qqu(i,jj)
                         end if
                      end do
                   end do
                   v(ivFlucq) = v(ivFlucq) - fqegp(imolty)
!> removed by BLE 01/15; don't want to zero out this energy overall just because one molecule type is not fluctuating
                !else
                   !v(ivFlucq) = 0.0E0_dp
                end if
             end if
          end if
       end do
! ------------------------------------------

! **************************************************
! CALCULATION OF VIB. + BEND. + TORS. ENERGY ***
! **************************************************

       ! NOTE here virtual coordinates can be used!!!
       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             call U_bonded(i,imolty,sum_vvib,sum_vbend,sum_vtg)
             v(ivStretching)=v(ivStretching)+sum_vvib
             v(ivBending)=v(ivBending)+sum_vbend
             v(ivTorsion)=v(ivTorsion)+sum_vtg

          end do
       end do
       call mp_sum(v(ivStretching),1,groupid)
       call mp_sum(v(ivBending),1,groupid)
       call mp_sum(v(ivTorsion),1,groupid)
    end if !if ( .not. lvol .or. (lvol .and. lewald) )

! ################################################################

! ***************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

! for adsorption isotherms, don't calculate energy w/surface
! in box 2
    if ((lelect_field).or.((ibox.eq.1).and.(lexzeo.or.lslit.or.lgraphite.or.lsami.or.lmuir))) then
       ! MPI
       do imolty=1,nmolty
          if (.not.lexclu_zeo(imolty)) then
             do nmcount=myid+1,ncmt(ibox,imolty),numprocs
                i=parbox(nmcount,ibox,imolty)
                do j = 1, nunit(imolty)
                   ntj = ntype(imolty,j)
                   v(ivExt)=v(ivExt)+U_ext(ibox,i,j,ntj)
                end do
             end do
          end if
       end do

       call mp_sum(v(ivExt),1,groupid)
    end if

! --------------------------------------------------------------------
! calculation of additional gaussian potential needed in thermodynamic
! integration in stages b and c
! --------------------------------------------------------------------
    if (lmipsw) then
       ! MPI
       do imolty=1,nmolty
          !> \bug most likely only works when nbox = 1
       do nmcount=myid+1,ncmt(ibox,imolty),numprocs
          i=parbox(nmcount,ibox,imolty)
          if (lwell(imolty)) then
             rxui = xcm(i)
             ryui = ycm(i)
             rzui = zcm(i)
             do j = 1, nwell(imolty)*nunit(imolty)
                k = j-int(j/nunit(imolty))*nunit(imolty)
                if (k.eq.0) k = nunit(imolty)
                rxuij = rxui-rxwell(j,imolty)
                ryuij = ryui-rywell(j,imolty)
                rzuij = rzui-rzwell(j,imolty)
                call mimage(rxuij,ryuij,rzuij,ibox)
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
                      call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                      vwell = vwell-awell(ii,k,imolty)*exp(-bwell*rijsq)
                   end do
                end if
             end do
          end if
       end do
       end do

       call mp_sum(vwell,1,groupid)
    end if

! ----------------------------------------------------------------------------

! write(io_output,*) 'self,corr:',(v(ivElect)-vrecipsum)*qqfact
! write(io_output,*) 'vsc, new self cor:',vsc*qqfact
! write(io_output,*) 'recip space part :',vrecipsum*qqfact
! write(io_output,*) 'sc and recip:',(vsc+vrecipsum)*qqfact

    if (.not.L_elect_table) then
       v(ivElect) = v(ivElect)*qqfact
    end if

    v(ivTot) = v(ivInterLJ) + v(ivIntraLJ) + v(ivExt) + v(ivElect) + v(ivFlucq) + v(iv3body)

! write(io_output,*) 'v in sumup',v(ivTot)

    vipsw = v(ivTot)
    vwellipsw = vwell

    if (lstagea) then
       v(ivTot) = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*v(ivTot)
    else if (lstageb) then
       v(ivTot) = etais*v(ivTot)+lambdais*vwell
    else if (lstagec) then
       v(ivTot) = (etais+(1.0E0_dp-etais)*lambdais)*v(ivTot)+(1.0E0_dp-lambdais)*vwell
    end if

    v(ivTot) = v(ivTot) + v(ivStretching) + v(ivBending) + v(ivTorsion)

    if ( .not. lvol.and.myid.eq.rootid ) then
       write(io_output,*)
       write(io_output,*) 'sumup control'
       write(io_output,*) 'number of chains', nchbox(ibox)
       do i = 1, nmolty
          write(io_output,"(A,I2,A,A10,I8)") 'number of chains of type   ',i,' ', molecname(i),ncmt(ibox,i)
       end do
       write(io_output,*) 'inter lj energy ', v(ivInterLJ)
       write(io_output,*) 'intra lj energy ', v(ivIntraLJ)
       if (ltailc) write(io_output,*) 'Tail correction ', v(ivTail)
       write(io_output,*) 'bond vibration  ', v(ivStretching)
       write(io_output,*) 'bond bending    ', v(ivBending)
       write(io_output,*) 'torsional       ', v(ivTorsion)
       write(io_output,*) 'external        ', v(ivExt)
       write(io_output,*) 'coulombic energy', v(ivElect)
       ! write(io_output,*) 'exact energy    ', 1.74756*1.67*831.441/3.292796
       write(io_output,*) 'fluc Q energy   ', v(ivFlucq)
       write(io_output,*) 'well energy     ', vwellipsw
       if(lgaro) write(io_output,*) '3-body garo     ', v(iv3body)
       write(io_output,*) 'total energy    ', v(ivTot)
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end SUMUP in ',myid
#endif

    return
  end subroutine sumup

!*****************************************************************
!> \brief Calculates the total potential energy for a configuration.
!>
!> \param i calculate the energies associated with chain i
!> \param imolty molecule type of chain i
!> \param v* energies
!> \param flagon flag for old(flagon=1)/new(flagon=2) configurations
!> \param ibox box number of chain i
!> \param istart, iuend calculate energies from bead istart to bead iuend for chain i
!> \param lljii whether to include intramolecular LJ interactions
!> \param ovrlap atom overlap
!> \param ltors whether to calculate torsional energy
!> \param lcharge_table whether need to set up charge interaction table; true if called from CBMC
!> \param lfavor
!*****************************************************************
  subroutine energy(i,imolty,v,flagon,ibox,istart,iuend,lljii,ovrlap,ltors,lcharge_table,lfavor,lAtom_traxyz)
    use sim_particle,only:lnn,lnn_t,neighbor,neigh_cnt,ndij,nxij,nyij,nzij,ctrmas,add_neighbor_list_molecule,neighi,neigh_icnt&
     ,ndiji,nxiji,nyiji,nziji
    use energy_intramolecular,only:U_torsion
    use energy_3body,only:U3MolSys
    use energy_4body,only:U4MolSys
    use energy_garofalini,only:triad_en

    logical::lqimol,lqjmol,lexplt,lcoulo(numax,numax),lfavor,lij2,liji,lqchgi
    logical::lljii,ovrlap,ltors,lcharge_table,lAtom_traxyz

    integer::growii,growjj,k,jcell(nmax),nmole
    integer::i,ibox,istart,iuend,ii,ntii,flagon,jjj,iii,mmm,j,jj,ntjj,ntij,ntj,imolty,jmolty,jjend
    integer::nchp2

    real::v(nEnergy),rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq,ryuij,rzuij,rij,rijsq,rbcut,calpi
    real::vwell
    real::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq

    real, allocatable :: lij_list(:,:)
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start ENERGY in ',myid,' for molecule ',i,' in box ',ibox
#endif

    if ( lpbc ) call setpbc(ibox)

    rbcut = rcut(ibox)
    rcutsq = rbcut * rbcut
    calpi = calp(ibox)
    if (ldual) rcinsq = rcutin*rcutin
    rminsq = rmin * rmin

    ovrlap = .false.
    !> changed by BLE (01/15) to be specific to the parts of energy that are calculated in this subroutine
    !> this avoids zeroing out the FQ energy for polarizable models anytime energy is called
    !v = 0.0E0_dp
    v(ivTot) = 0.0E0_dp
    v(ivInterLJ) = 0.0E0_dp
    v(ivIntraLJ) = 0.0E0_dp
    v(ivExt) = 0.0E0_dp
    v(ivElect) = 0.0E0_dp
    v(ivEwald) = 0.0E0_dp
    v(ivTorsion) = 0.0E0_dp
    v(iv3body) = 0.0E0_dp
    sself  = 0.0E0_dp
    correct = 0.0E0_dp

    if (flagon.eq.2) then
       if (lneigh) lnn_t(:,i)=.false.
       if (lneighbor.or.lgaro) neigh_icnt = 0
    end if

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
    lqimol = lelect(imolty)

    if (nugrow(imolty) .eq. nunit(imolty)) then
       lexplt = .false.
    else
       lexplt = .true.
       growii = nugrow(imolty)
    end if

    nchp2 = nchain + 2
    do ii = 1,nunit(imolty)
       rxu(nchp2,ii) = rxuion(ii,flagon)
       ryu(nchp2,ii) = ryuion(ii,flagon)
       rzu(nchp2,ii) = rzuion(ii,flagon)
       qqu(nchp2,ii) = qquion(ii,flagon)
    end do
    nboxi(nchp2) = ibox
    moltyp(nchp2) = imolty

    if ( lcutcm .or. lfavor ) then
       ! calculate the center of mass of chain i and give it a dummy #
       call ctrmas(.false.,ibox,nchp2,9)
       xcmi = xcm(nchp2)
       ycmi = ycm(nchp2)
       zcmi = zcm(nchp2)
       rcmi = rcmu(nchp2)
       ! write(io_output,*) 'rcmi:',rcmi
    else
       lij2 = .true.
    end if

    if (licell.and.(ibox.eq.boxlink)) then
       ii=1
       rxui = rxuion(ii,flagon)
       ryui = ryuion(ii,flagon)
       rzui = rzuion(ii,flagon)
       call get_cell_neighbors(rxui,ryui,rzui,ibox,jcell,nmole)
    else
       nmole = nchain
    end if

! loop over all chains except i
! JLR 11-24-09 don't loop if box is ideal gas
    if (.not.(lideal(ibox))) then
! END JLR 11-24-09
! RP added for MPI
       if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
           call energy_inter_kd_tree_energy(i, imolty, v, flagon, ibox, istart, iuend, ovrlap, lcoulo)
           if (ovrlap) return
       else
           ! compute the COM distances if COM-based kdtree is used
           if (lcutcm .and. lkdtree .and. lkdtree_box(ibox)) then
               call energy_inter_com_kd_tree(ibox, i, xcmi, ycmi, zcmi, rbcut, rmin &
                    , .false., lij_list, .false., ovrlap)
               if (ovrlap) return
           end if

           do k = myid+1,nmole,numprocs
! do k = 1, nmole
              if (licell.and.(ibox.eq.boxlink)) then
                 j = jcell(k)
              else
                 j = k
              end if

              jmolty = moltyp(j)
              lqjmol = lelect(jmolty)
              growjj = nugrow(jmolty)

! check for simulation box ###
              if ( ( ibox .eq. nboxi(j) ) .and. (i .ne. j )) then
                 if ( lneigh ) then
                    if ( .not. lnn(j,i) ) cycle
                 end if

                 if (lcutcm .and. lkdtree .and. lkdtree_box(ibox)) then
                     if (lij_list(1, j) .lt. 0.5d0) then !< if outside com-cutoff, cycle the do loop on mlcls
                         cycle
                     else
                         lij2 = .true.
                         goto 2047 !< skip the following brute-force com-cutoff calculation
                     end if
                 end if

                 if (lcutcm .or. lfavor) then
                    ! check if ctrmas within rcmsq
                    rxuij = xcmi - xcm(j)
                    ryuij = ycmi - ycm(j)
                    rzuij = zcmi - zcm(j)
                    ! minimum image the ctrmas pair separations ***
                    if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                    rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                    rcm = rbcut + rcmi + rcmu(j)
                    rcmsq = rcm*rcm
                    ! write(io_output,*) rcm,rcmi,rcmu(j)
                    if ( lfavor ) then
                        favor(j) = (rminsq/rijsq)**2*5.0E0_dp
                        favor2(j) = rminsq/rijsq
                    end if
                    if ( rijsq .gt. rcmsq .and. lcutcm) then
                        if ( lqimol .and. lqjmol .and. lchgall ) then
                            lij2 = .false.
                        else
                            cycle
                        end if
                    else
                        lij2 = .true.
                    end if
                end if

2047            continue

                if ( lcharge_table .and. (.not. lchgall) ) then
                    ! called from CBMC and must set up charge-interaction table ---
                    do ii = 1,nugrow(imolty)
                        do jj = 1,nugrow(jmolty)
                            iii = leaderq(imolty,ii)
                            jjj = leaderq(jmolty,jj)
                            if ( iii .eq. ii .and. jjj .eq. jj ) then
                                rxuij = rxuion(ii,flagon) - rxu(j,jj)
                                ryuij = ryuion(ii,flagon) - ryu(j,jj)
                                rzuij = rzuion(ii,flagon) - rzu(j,jj)
                                if ( lpbc )  call mimage(rxuij,ryuij,rzuij,ibox)
                                rijsq = rxuij*rxuij + ryuij*ryuij  + rzuij*rzuij
                                if ((rijsq .lt. rcutsq) .or. lijall) then
                                    lcoulo(ii,jj) = .true.
                                else
                                    lcoulo(ii,jj) = .false.
                                end if
                            end if
                        end do
                    end do
                end if

                ! loop over all beads ii of chain i
                do ii = istart, iuend
                    ntii = ntype(imolty,ii)
                    liji = lij(ntii)
                    lqchgi = lqchg(ntii)
                    rxui = rxuion(ii,flagon)
                    ryui = ryuion(ii,flagon)
                    rzui = rzuion(ii,flagon)

                    ! loop over all beads jj of chain j
                    do jj = 1, nunit(jmolty)
                        ! check exclusion table
                        if ( lexclu(imolty,ii,jmolty,jj) ) cycle

                        ntjj = ntype(jmolty,jj)
                        if ( lij2 ) then
                            if ((.not.(liji.and.lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj))))  cycle
                        else
                            if (.not.(lqchgi.and.lqchg(ntjj))) cycle
                        end if

                        ntij=type_2body(ntii,ntjj)

                        if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                        rxuij = rxui - rxu(j,jj)
                        ryuij = ryui - ryu(j,jj)
                        rzuij = rzui - rzu(j,jj)
                        ! minimum image the pair separations ***
                        if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                        rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                        if (rijsq.lt.rminsq .and. .not.(lexpand(imolty).or.lexpand(jmolty))) then
                            ovrlap = .true.
                        ! RP added for MPI
                        ! return
                            goto 99
                        else if ((rijsq .lt. rcutsq) .or. lijall) then
                            v(ivInterLJ)=v(ivInterLJ)+U2(rijsq,nchp2,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                        end if

                        v(ivElect)=v(ivElect)+Q2(rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo,ibox)

                        if (lneigh.and.flagon.eq.2.and.rijsq.le.rcutnn(ibox)**2) lnn_t(j,i)=.true.

                        ! KM lneighbor and lgaro does not work in parallel
                        !> \bug Wouldn't testing for center-of-mass distance be a better criteria?
                        if (lneighbor.and.ii.eq.1.and.jj.eq.1.and.flagon.eq.2.and.rijsq.lt.rbsmax**2.and.rijsq.gt.rbsmin**2) then
                            call add_neighbor_list_molecule(j,sqrt(rijsq),rxuij,ryuij,rzuij)
                        else if (lgaro.and.flagon.eq.2) then
                            if (isNeighbor(rijsq,ntii,ntjj)) then
                                call add_neighbor_list_molecule(j,sqrt(rijsq),rxuij,ryuij,rzuij)
                            end if
                        end if
                    end do
                end do
            end if
          end do

! RP added for MPI
! Returning from ovrlap--------------
99        continue

          call mp_lor(ovrlap,1,groupid)
          if(ovrlap)then
              ! write(io_output,*)'630 in energy ovrlap=',ovrlap,'myid=',myid
              return
          end if

          call mp_sum(v(ivInterLJ),1,groupid)
          call mp_sum(v(ivElect),1,groupid)

        end if !< if kdtree
    end if !< if box is ideal

! -----------------------------------------------

!kea - garo: add three body loop for intermolecular interactions
    if (lgaro.and..not.lideal(ibox)) then
       if(flagon.eq.2) then
          v(iv3body)=triad_en(i,neigh_icnt,neighi,ndiji,nxiji,nyiji,nziji,.true.)
       else if (flagon.eq.1) then
          v(iv3body)=triad_en(i,neigh_cnt(i),neighbor(:,i),ndij(:,i),nxij(:,i),nyij(:,i),nzij(:,i),.false.)
       end if
    end if

    if (hasThreeBody.and..not.lideal(ibox)) v(ivInterLJ)=v(ivInterLJ)+U3MolSys(i,istart,iuend,flagon)
    if (hasFourBody.and..not.lideal(ibox)) v(ivInterLJ)=v(ivInterLJ)+U4MolSys(i,istart,iuend,flagon)

! ################################################################

! the intramolecular van der waals and ewald terms have to be calculated
! for the explicit atom placement models
! *******************************
! INTRACHAIN INTERACTIONS ***
! *******************************

! JLR 11-19-09 commenting this out, alway do mimage for intrachain
! for expanded ensemble
! lmim = .false.
! mlen2 = rcmu(nchp2)*2E0_dp
! if ( mlen2>boxlx(ibox) .or. mlen2>boxly(ibox) .or. mlen2>boxlz(ibox)) lmim = .true.
! END JLR 11-19-09

! calculate intramolecular energy correction for chain i
    do ii = istart, iuend
       ntii = ntype(imolty,ii)
       rxui = rxuion(ii,flagon)
       ryui = ryuion(ii,flagon)
       rzui = rzuion(ii,flagon)
       if (lAtom_traxyz) then
          jjend=nunit(imolty)
       else
          jjend=ii-1
       end if
       do jj = 1,jjend
          if (jj.eq.ii) cycle
          ntjj = ntype(imolty,jj)

          ntij=type_2body(ntii,ntjj)

          if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
!KNS 11-23-15 insert rxui, ryui, rzui determined above
          rxuij = rxui - rxuion(jj,flagon)
          ryuij = ryui - ryuion(jj,flagon)
          rzuij = rzui - rzuion(jj,flagon)

          ! lpbc is not called here because it's intra-chain interaction
          rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

          if (lqinclu(imolty,ii,jj) ) then
             ! calculation of intramolecular electrostatics
             v(ivElect)=v(ivElect)+qscale2(imolty,ii,jj)*Q2(rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchg(ntii),nchp2,imolty,jj&
              ,ntjj,calpi,lcoulo,ibox)
          end if

          ! calculation of other non-bonded interactions
          if ( linclu(imolty,ii,jj) ) then
             if (lljii) then
                if (rijsq.lt.rminsq .and. .not.lexpand(imolty)) then
                   ovrlap = .true.
                   ! write(io_output,*) 'intra ovrlap:',ii,jj
                   return
                else
                   ! there's no cutoff when calculating the intra-chain interaction
                   !> \bug skip intra if it is bending 1-3 and using a table?
                   if (L_bend_table) then
                      do mmm=1,inben(imolty,ii)
                         if (ijben3(imolty,ii,mmm).eq.jj) then
                            v(ivIntraLJ) = v(ivIntraLJ) + lininter_bend(sqrt(rijsq),itben(imolty,ii,mmm))
                            goto 96
                         end if
                      end do
                   end if

                   v(ivIntraLJ)=v(ivIntraLJ)+U2(rijsq,nchp2,imolty,ii,ntii,nchp2,imolty,jj,ntjj,ntij)
                end if
             end if
          end if

96        if (lewald.and..not.lideal(ibox)) then
             ! compute the ewald intramolecular (self and correction) terms for
             ! the interactions of the placed atoms with themselves, and with the
             ! rest of their own molecule, if there's no interaction
             ! these are 1,2 and 1,3
             rij = sqrt(rijsq)
             if (.not. lqinclu(imolty,ii,jj)) then
                correct=correct+qquion(ii,flagon)*qquion(jj,flagon)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                ! 1,4 interaction which we scale by qscale
             else
                correct=correct+(1.0E0_dp-qscale2(imolty,ii,jj))*qquion(ii,flagon)*qquion(jj,flagon)&
                 *(erfunc(calpi*rij)-1.0E0_dp)/rij
             end if
          end if
       end do
       if (lewald.and..not.lideal(ibox)) then
          sself = sself + qquion(ii,flagon)*qquion(ii,flagon)
       end if
    end do
    if (lewald.and..not.lideal(ibox)) then
       sself = -sself * calpi/sqrtpi
       v(ivEwald) = sself + correct
    end if
! ################################################################

! ***************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

    if ((lelect_field.and.lqimol).or.((ibox.eq.1).and.(lexzeo.or.lslit.or.lgraphite.or.lsami.or.lmuir))) then
       if (.not.lexclu_zeo(imolty)) then
          do j = istart,iuend
             ntj = ntype(imolty,j)
             v(ivExt)=v(ivExt)+U_ext(ibox,nchp2,j,ntj)
          end do
       end if
    end if

! *********************************************************************
! calculation of torsion energy for explicit atom methyl groups ****
! *********************************************************************
    if ( ltors ) then
       v(ivTorsion)=U_torsion(nchp2,imolty,nugrow(imolty)+1,.true.)
    end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
    vwell = 0.0E0_dp
    if (lwell(imolty).and.lmipsw) then
       rxui = xcmi
       ryui = ycmi
       rzui = zcmi
       do j = 1, nwell(imolty)*nunit(imolty)
          k = j - int(j/nunit(imolty))*nunit(imolty)
          if (k.eq.0) k = nunit(imolty)
          rxuij = rxui-rxwell(j,imolty)
          ryuij = ryui-rywell(j,imolty)
          rzuij = rzui-rzwell(j,imolty)
          call mimage(rxuij,ryuij,rzuij,ibox)
          rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
          rcm = rcut(ibox)+rcmi
          rcmsq = rcm*rcm
          if (rijsq.lt.rcmsq) then
             do ii = 1, nunit(imolty)
                if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                rxui = rxuion(ii,flagon)
                ryui = ryuion(ii,flagon)
                rzui = rzuion(ii,flagon)
                rxuij = rxui-rxwell(j,imolty)
                ryuij = ryui-rywell(j,imolty)
                rzuij = rzui-rzwell(j,imolty)
                call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                vwell = vwell-awell(ii,k,imolty)*exp(-bwell*rijsq)
             end do
          end if
       end do
    end if

! ----------------------------------------------------------------------------

    if (.not.L_elect_table) then
       v(ivElect) = v(ivElect)*qqfact
       v(ivEwald) = v(ivEwald)*qqfact
    end if

! note that vintra is only computed when the flag lljii is true
    v(ivTot) = v(ivInterLJ) + v(ivExt) + v(ivIntraLJ) + v(ivElect) + v(ivEwald) + v(iv3body)
! write(io_output,*) 'vinter:',v(ivInterLJ),'vext:',v(ivExt),'vintra:',v(ivIntraLJ),'velect',v(ivElect),'vewald:',v(ivEwald),'v'

    if (flagon.eq.1) then
       vipswo = v(ivTot)
       vwellipswo = vwell
    else
       vipswn = v(ivTot)
       vwellipswn = vwell
    end if

    if (lmipsw) then
       if (lstagea) then
          v(ivTot) = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*v(ivTot)
       else if (lstageb) then
          v(ivTot) = etais*v(ivTot)+lambdais*vwell
       else if (lstagec) then
          v(ivTot) = (etais+(1.0E0_dp-etais)*lambdais)*v(ivTot)+(1.0E0_dp-lambdais)*vwell
       end if
    end if

#ifdef __DEBUG__
    ! write(io_output,*) 'v :', v(ivTot)
    write(io_output,*) 'end ENERGY in ',myid
#endif
    return
  end subroutine energy

!*****************************************************************
!> \brief Calculates the potential energy and the boltzmann factor
!>       for ichoi trial positions.
!>
!> \param lnew true for new configurations
!> \param lfirst true for insertion of the first bead in swap moves
!> \param ovrlap logical variable, true for walk termination
!> \param i calculates Boltzmann weights for the newly-grown beads in chain i
!> \param icharge usually identical to i
!> \param imolty molecule type of chain i
!> \param ibox box number of chain i
!> \param ichoi number of trial positions
!> \param iufrom the bead from which the new beads were grown
!> \param ntogrow number of new beads that have been grown
!> \param glist the list of new beads that have been grown; ntogrow entries
!> \param maxlen maximum possible distance of the newly-grown beads from iufrom
!*****************************************************************
  subroutine boltz(lnew,lfirst,ovrlap,i,icharge,imolty,ibox,ichoi,iufrom,ntogrow,glist,maxlen)
    use sim_particle,only:lnn
    use util_mp,only:mp_set_displs

    logical::lnew,ovrlap,lcmno(nmax),lfirst
    logical::lqimol,lqjmol,liji,lqchgi
    integer::ichoi,growjj,igrow,count,glist(numax),icharge,cnt,jcell(nmax)
    integer::i,imolty,ibox,ntogrow,itrial,ntii,j,jj,ntjj,ntij,iu,jmolty,iufrom,ii,k,nmole
    ! integer::NRtype
    real::rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij,rij,rijsq,maxlen,rcm,rcmsq,corr,rcutmax,rbcut
    real::v(nEnergy),vwell,rcutsq,rcinsq
    integer::mmm

    ! RP added for MPI
    integer::rcounts(numprocs),displs(numprocs),my_start,my_end,blocksize,my_itrial
    real::my_vtry(nchmax),my_vtrintra(nchmax),my_vtrext(nchmax),my_vtrinter(nchmax),my_vtrelect(nchmax),my_vtrewald(nchmax)&
     ,my_bfac(nchmax),my_vipswot(nchmax),my_vwellipswot(nchmax),my_vipswnt(nchmax),my_vwellipswnt(nchmax)
    logical::my_lovr(nchmax)

    real, allocatable :: lij_list(:,:)
! ------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start BOLTZ in ',myid
#endif

    if ( lpbc ) call setpbc(ibox)

    ! determine the potential cutoffs
    rcutsq = rcut(ibox)*rcut(ibox)
    rbcut = rcut(ibox)

    ! KM initialize variables
    do j=1,ichoi
       my_lovr(j) = .false.
       lovr(j) = .false.
       my_vtry(j) = 0.0E0_dp
       my_vtrintra(j) = 0.0E0_dp
       my_vtrelect(j) = 0.0E0_dp
       my_vtrext(j) = 0.0E0_dp
       my_vtrinter(j) = 0.0E0_dp
       my_vtrewald(j) = 0.0E0_dp
       my_bfac(j) = 0.0E0_dp
       my_vipswot(j) = 0.0E0_dp
       my_vwellipswot(j) = 0.0E0_dp
       my_vipswnt(j) = 0.0E0_dp
       my_vwellipswnt(j) = 0.0E0_dp
    end do

    if ( ldual ) then
       ! use rcutin for both types of interactions (except intra)
       rbcut = rcutin
       rcinsq = rcutin*rcutin
    else
       ! compute the cutoffs squared for each interaction
       rcinsq =  rcutsq
       if ( lcutcm ) then
          ! not needed when ldual is true since will use rcutin then
          rcutmax = rcut(ibox)
       end if
    end if

    ! compute minimum cutoff squared
    rminsq = rmin * rmin

    lqimol = lelect(imolty)
    igrow = nugrow(imolty)

    if ( lcutcm .and. (.not. lfirst) ) then
       ! check previous bead (iufrom) COM for each molecule in the box
       if ( lnew ) then
          ! for trial chain ###
          rxui  = rxnew(iufrom)
          ryui  = rynew(iufrom)
          rzui  = rznew(iufrom)
       else
          ! for old chain ###
          rxui  = rxu(i,iufrom)
          ryui  = ryu(i,iufrom)
          rzui  = rzu(i,iufrom)
       end if

!> \todo to be MPI parallelized
!> COM distance calculation
       if (lkdtree .and. lkdtree_box(ibox)) then
           call energy_inter_com_kd_tree(ibox, i, rxui, ryui, rzui, rbcut, 0.0, lfirst, lij_list, .false., ovrlap)

           ! update lij_list to lcmno
           do j = 1, nchain
               if (lij_list(1, j) .lt. 0.5d0) then
                   lcmno(j) = .true.
               else
                   lcmno(j) = .false.
               end if
           end do
       else
           do j = 1,nchain
              lcmno(j) = .false.
              if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
                 rxuij = rxui-xcm(j)
                 ryuij = ryui-ycm(j)
                 rzuij = rzui-zcm(j)
                 ! minimum image the pseudo-ctrmas pair separation
                 if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)

                 rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

                 if ( ldual ) then
                    rcm = rcutin + rcmu(j) + maxlen
                    rcmsq = rcm*rcm
                 else
                    rcm = rcutmax + rcmu(j) + maxlen
                    rcmsq = rcm*rcm
                 end if
                 if (rijsq .gt. rcmsq ) lcmno(j) = .true.
              end if
           end do
       end if !< if not lkdtree
    end if

    ! RP added for MPI
    blocksize = ichoi/numprocs
    rcounts = blocksize
    blocksize = ichoi - blocksize * numprocs
    if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
    call mp_set_displs(rcounts,displs,blocksize,numprocs)
    my_start = displs(myid+1) + 1
    my_end = my_start + rcounts(myid+1) - 1
    !if (ldebug) write(myid+100,*)'boltz: my_start=',my_start,'; my_end=',my_end,'; ichoi=',ichoi,'; rcounts=',rcounts,'; displs=',displs

    my_itrial = 0
    ! do itrial = 1, ichoi
    do itrial = my_start,my_end
       my_itrial = my_itrial + 1
       ! lovr(itrial) = .false.
       my_lovr(my_itrial) = .false.

       v(ivInterLJ) = 0.0E0_dp
       v(ivIntraLJ) = 0.0E0_dp
       v(ivExt) = 0.0E0_dp
       v(ivElect) = 0.0E0_dp
       v(ivEwald) = 0.0E0_dp

       ! Only if L_Coul_CBMC is true, then compute electrostatic interactions/corrections
       if(L_Coul_CBMC.and.lewald.and..not.lideal(ibox)) then
          do count = 1,ntogrow
             ii = glist(count)
             ! This part does not change for fixed charge moves, but is
             ! used in the swap rosenbluth weight. - ewald self term
             v(ivEwald) = v(ivEwald) - qqu(icharge,ii)*qqu(icharge,ii)*calp(ibox)/sqrtpi
          end do
       end if

       ! no intramolecular interactions if this is the first bead
       if ( .not. lfirst ) then
! *****************************************
! INTRACHAIN BEAD-BEAD INTERACTIONS ***
! *****************************************
          ! cycle through molecule and check bead by bead
          do iu = 1, igrow
             ! see if iu exists in the new chain yet
             if (.not. lexist(iu)) cycle
             ntjj = ntype(imolty,iu)

             ! loop over all the grown beads
             do count = 1,ntogrow
                ii = glist(count)
                ! assign bead type for ii,iu, and the cross term
                ntii = ntype(imolty,ii)
                ntij=type_2body(ntii,ntjj)

                ! see if iu has nonbonded intramolecular interaction with ii
                if (linclu(imolty,ii,iu).or.(L_Coul_CBMC.and.(lqinclu(imolty,ii,iu).or.lewald))) then

                   if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
                   ! determine distances
                   if ( lnew ) then
                      ! use new trial chain coordinates
                      rxuij  = rxnew(iu) - rxp(count,itrial)
                      ryuij  = rynew(iu) - ryp(count,itrial)
                      rzuij  = rznew(iu) - rzp(count,itrial)
                   else
                      ! use old chain coordinates
                      rxuij  = rxu(i,iu) - rxp(count,itrial)
                      ryuij  = ryu(i,iu) - ryp(count,itrial)
                      rzuij  = rzu(i,iu) - rzp(count,itrial)
                   end if
                   ! lpbc is not called here because it's intra-chain interaction
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                end if
                if ( linclu(imolty,ii,iu) .or. lqinclu(imolty,ii,iu)) then
                   if ( linclu(imolty,ii,iu) ) then
                      if (rijsq.lt.rminsq.and..not.lexpand(imolty)) then
                         ! RP added for MPI
                         my_lovr(my_itrial) = .true.
                         ! write(io_output,*) 'intra overlap'
                         goto 19
                      else
                         ! there's no cutoff when calculating the intra-chain interaction
                         !> \bug skip intra if it is bending 1-3 and using a table?
                         if (L_bend_table) then
                            do mmm=1,inben(imolty,ii)
                               if (ijben3(imolty,ii,mmm).eq.iu) then
                                  v(ivIntraLJ) = v(ivIntraLJ) + lininter_bend(sqrt(rijsq),itben(imolty,ii,mmm))
                                  goto 96
                               end if
                            end do
                         end if

                         v(ivIntraLJ)=v(ivIntraLJ)+U2(rijsq,i,imolty,ii,ntii,i,imolty,iu,ntjj,ntij)
                      end if
                   end if

                   ! intramolecular charge interaction
                   ! compute velect (coulomb and ewald)
96                 if (L_Coul_CBMC.and.lqinclu(imolty,ii,iu).and.lqchg(ntii).and.lqchg(ntjj).and.rijsq.lt.rcinsq) then
                      ! boltz.f has problem to compute the electrostatic interactions
                      ! in a group-based way because the leader q might not be grown at
                      ! present, so it calculates electrostatic interaction not based on
                      ! group but on its own distance in SC, but should be corrected
                      ! later by calling energy subroutine.
                      rij = sqrt(rijsq)
                      if (L_elect_table) then
                         v(ivElect) = v(ivElect) + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(icharge,iu)&
                          *lininter_elect(rij,ntii,ntjj)
                      else if (lewald.and.(.not.lideal(ibox))) then
                         ! compute real space term of vewald
                         v(ivElect) = v(ivElect) + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(icharge,iu)*erfunc(calp(ibox)*rij)&
                          / rij
                         ! ewald sum correction term
                         corr = (1.0E0_dp - qscale2(imolty,ii,iu))*qqu(icharge,ii)*qqu(icharge,iu)&
                          *(erfunc(calp(ibox)*rij)-1.0E0_dp)/rij
                         v(ivEwald) = v(ivEwald) + corr
                      else
                         v(ivElect) = v(ivElect) + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(i,iu)/rij
                      end if
                   end if
                   ! end charge calculation
                else if (L_Coul_CBMC.and.lewald) then
                   ! will only add correction here if lqinclu is false.
                   ! ewald sum correction term
                   rij = sqrt(rijsq)
                   corr = qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0E0_dp) /rij
                   v(ivEwald) = v(ivEwald) + corr
                end if
             end do
          end do
          !> \bug double Ewald correction?
          if (L_Coul_CBMC.and.lewald.and.ntogrow.gt.1.and..not.lideal(ibox)) then
             ! ewald sum correction term for interactions of the
             ! growing beads with each other
             ! this is 1-3, so don't need to consult lqinclu
             ! should change this since it corrects for all currently grown beads in this
             ! step, which could at somepoint be further than 1-3 apart!!! (say for rigrot...)
             do cnt = 1,ntogrow-1
                iu = glist(cnt)
                do count = cnt+1,ntogrow
                   ii = glist(count)
                   ! determine distances - use trial chain coordinates
                   rxuij  = rxp(cnt,itrial) - rxp(count,itrial)
                   ryuij  = ryp(cnt,itrial) - ryp(count,itrial)
                   rzuij  = rzp(cnt,itrial) - rzp(count,itrial)
                   !> \bug no mimage convertion?
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   rij   = sqrt(rijsq)
                   ! ewald sum correction term
                   corr = qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0E0_dp)/rij
                   v(ivEwald) = v(ivEwald) + corr
                end do
             end do
          end if
       end if

! JLR 11-24-09 don't compute if lideal
       if (.not.lideal(ibox)) then
! END JLR 11-24-09
          if (licell.and.(ibox.eq.boxlink)) then
! we only use count = 1, the rest should be taken care of
! with rintramax
             count = 1
             rxui = rxp(count,itrial)
             ryui = ryp(count,itrial)
             rzui = rzp(count,itrial)
             call get_cell_neighbors(rxui,ryui,rzui,ibox,jcell,nmole)
          else
             nmole = nchain
          end if

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
! loop over all chains except i
          if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
              call energy_inter_kd_tree_boltz(i, icharge, imolty, v, ibox, ntogrow, ovrlap, glist, itrial)

              if (ovrlap) then
                  my_lovr(my_itrial) = .true.
                  goto 19
              end if
          else !< if kdtree

              ! compute COM distances if lfirst (where COM distances have not been computed) and COM-kdtree is used
              if (lcutcm .and. lfirst .and. lkdtree .and. lkdtree_box(ibox)) then
                 call energy_inter_com_kd_tree(ibox, i, rxp(1,itrial), ryp(1,itrial), rzp(1,itrial) &
                    , rbcut, 0.0, lfirst, lij_list, .false., ovrlap)
              end if

              do_nmole:do k = 1, nmole
                 if (licell.and.(ibox.eq.boxlink)) then
                    j = jcell(k)
                 else
                    j = k
                 end if

                 ! check for simulation box
                 if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
                    ! check neighbor list
                    if (lneigh.and..not.lnew) then
                       if (.not.lnn(j,i)) cycle do_nmole
                    end if
                    ! check COM table calculated above
                    if (.not.lfirst.and.lcmno(j).and.lcutcm) cycle do_nmole

                    jmolty = moltyp(j)
                    lqjmol = lelect(jmolty)
                    growjj = nugrow(jmolty)

                    ! loop over all beads of molecule i grown this step
108                 do count = 1,ntogrow
                       ! assign bead type for ii
                       ii = glist(count)
                       ntii = ntype(imolty,ii)
                       liji = lij(ntii)
                       lqchgi = lqchg(ntii)

                       ! assign positions to r*ui
                       rxui = rxp(count,itrial)
                       ryui = ryp(count,itrial)
                       rzui = rzp(count,itrial)

                       if ( lfirst .and. lcutcm ) then
                          ! check if ctrmas within rcmsq
                          if (lkdtree .and. lkdtree_box(ibox)) then
                             if (lij_list(1, j) .lt. 0.5d0) cycle do_nmole
                          else
                             rxuij = rxui-xcm(j)
                             ryuij = ryui-ycm(j)
                             rzuij = rzui-zcm(j)
                             ! minimum image the ctrmas pair separations
                             if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                             rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                             ! determine cutoff
                             if ( ldual ) then
                                ! must be lfirst so no previous bead
                                rcm = rcutin + rcmu(j)
                             else
                                ! standard lcutcm cutoff
                                rcm = rcutmax + rcmu(j)
                             end if
                             ! check if interaction distance is greater than cutoff
                             rcmsq = rcm*rcm
                             if ( rijsq .gt. rcmsq ) cycle do_nmole
                          end if
                       end if

                       ! loop over all beads jj of chain j
                       do jj = 1,nunit(jmolty)
                          ! check exclusion table
                          if ( lexclu(imolty,ii,jmolty,jj) ) cycle
                          ! start iswatch add-on ***
                          !> \todo is there a way to pull this out of the loops?
                          ! if (liswatch.and.j.eq.other.and.(.not.liswinc(jj,jmolty))) then
                          !    cycle
                          ! end if
                          ! end iswatch add-on ***

                          ntjj = ntype(jmolty,jj)
                          if ((.not.(liji.and.lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj)))) cycle

                          ntij=type_2body(ntii,ntjj)

                          if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                          rxuij = rxui - rxu(j,jj)
                          ryuij = ryui - ryu(j,jj)
                          rzuij = rzui - rzu(j,jj)
                          ! minimum image the pair separations ***
                          if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                          rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                          ! compute vinter (eg. lennard-jones)
                          if (rijsq.lt.rminsq.and..not.(lexpand(imolty).or.lexpand(jmolty))) then
                             my_lovr(my_itrial) = .true.
                             ! write(io_output,*) 'j:',j,jj
                             ! write(io_output,*) 'rjsq:',rijsq,rminsq
                             goto 19
                          else if (rijsq.lt.rcinsq.or.lijall) then
                             v(ivInterLJ)=v(ivInterLJ)+U2(rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                          end if

                          ! compute velect (coulomb and ewald)
                          if (L_Coul_CBMC.and.lqchg(ntii).and.lqchg(ntjj).and.rijsq.lt.rcinsq) then
                             ! boltz.f has problem to compute the electrostatic interactions
                             ! in a group-based way because the leader q might not be grown at
                             ! present, so it calculates electrostatic interaction not based on
                             ! group but on its own distance in SC, but should be corrected
                             ! later by calling energy subroutine.
                             rij = sqrt(rijsq)
                             if (L_elect_table) then
                                v(ivElect) = v(ivElect) + qqu(icharge,ii)*qqu(j,jj)*lininter_elect(rij,ntii,ntjj)
                             else if (lewald) then
                                ! compute real space term of velect
                                v(ivElect) = v(ivElect) + qqu(icharge,ii)*qqu(j,jj)*erfunc(calp(ibox)*rij)/rij
                             else
                                ! compute all electrostatic interactions
                                v(ivElect) = v(ivElect) + qqu(icharge,ii)*qqu(j,jj)/ rij
                             end if
                          end if
                       end do
                    end do
                 end if
              end do do_nmole
          end if !< lkdtree
       end if
! ################################################################

! **************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

! not for grand can. with ibox=2 !
! required for histogram reweighting to work for monolayer
! phase diagrams.
! not used for adsorption isotherms
       if ((ibox .eq. 1).and.(.not.lexclu_zeo(imolty))) then
          if ((lelect_field.and.lqimol).or.((ibox.eq.1).and.(lexzeo.or.lslit.or.lgraphite.or.lsami.or.lmuir))) then
             do count = 1,ntogrow
                ! assign bead type for ii
                ii = glist(count)
                ntii = ntype(imolty,ii)
                rxu(nchain+2,ii) = rxp(count,itrial)
                ryu(nchain+2,ii) = ryp(count,itrial)
                rzu(nchain+2,ii) = rzp(count,itrial)
                v(ivExt)=v(ivExt)+U_ext(ibox,nchain+2,ii,ntii)
             end do
          end if
       end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
       vwell = 0.0E0_dp
       if (lwell(imolty).and.lmipsw) then
          rxui = xcm(i)
          ryui = ycm(i)
          rzui = zcm(i)
          do j = 1, nwell(imolty)*nunit(imolty)
             k = j - int(j/nunit(imolty))*nunit(imolty)
             if (k.eq.0) k = nunit(imolty)
             rxuij = rxui-rxwell(j,imolty)
             ryuij = ryui-rywell(j,imolty)
             rzuij = rzui-rzwell(j,imolty)
             call mimage(rxuij,ryuij,rzuij,ibox)
             rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
             rcm = rcut(ibox)+rcmu(i)+maxlen
             rcmsq = rcm*rcm
             if (rijsq.lt.rcmsq) then
                do count = 1, ntogrow
                   ii = glist(count)
                   if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                   rxui = rxp(count,itrial)
                   ryui = ryp(count,itrial)
                   rzui = rzp(count,itrial)
                   rxuij = rxui-rxwell(j,imolty)
                   ryuij = ryui-rywell(j,imolty)
                   rzuij = rzui-rzwell(j,imolty)
                   call mimage(rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                   vwell = vwell-awell(ii,k,imolty)*exp(-bwell*rijsq)
                end do
             end if
          end do
       end if
! ----------------------------------------------------------------------------

! *********************************************
! CALCULATION OF TOTAL POTENTIAL ENERGY ***
! *********************************************
! write(23,*)
! write(23,*) 'wirting out total energy'
! write(23,*) 'Lnew', lnew
! if (NRtype.eq.1) then
! write(23,*) 'called from swap'
! else
! write(23,*) 'called from rigrot'
! end if

19     if ( my_lovr(my_itrial) ) then
          my_bfac(my_itrial) = 0.0E0_dp
       else
          if (.not.L_elect_table) then
             v(ivElect) = v(ivElect)*qqfact
             v(ivEwald) = v(ivEwald)*qqfact
          end if
          v(ivTot) = v(ivInterLJ)+v(ivIntraLJ)+v(ivExt)+v(ivElect)+v(ivEwald)
          if (.not.lnew) then
             my_vipswot(my_itrial) = v(ivTot)
             my_vwellipswot(my_itrial) = vwell
          else
             my_vipswnt(my_itrial) = v(ivTot)
             my_vwellipswnt(my_itrial) = vwell
          end if

          if (lstagea) then
             v(ivTot) = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*v(ivTot)
          else if (lstageb) then
             v(ivTot) = etais*v(ivTot)+lambdais*vwell
          else if (lstagec) then
             v(ivTot) = (etais+(1.0E0_dp-etais)*lambdais)*v(ivTot)+ (1.0E0_dp-lambdais)*vwell
          end if

          my_vtry(my_itrial) = v(ivTot)
          my_vtrintra(my_itrial) = v(ivIntraLJ)
          my_vtrext(my_itrial)   = v(ivExt)
          my_vtrinter(my_itrial) = v(ivInterLJ)
          my_vtrelect(my_itrial) = v(ivElect)
          my_vtrewald(my_itrial) = v(ivEwald)
          ! write(23,*) 'itrial' ,itrial
          ! write(23,*) vtr(ivTot,itrial), vtr(ivIntraLJ,itrial), vtr(ivExt,itrial),
          !   &       vtr(ivInterLJ,itrial),vtr(ivElect,itrial), vtr(ivEwald,itrial)

          if ((my_vtry(my_itrial)*beta).gt.(2.3E0_dp*softcut))then
             ! write(io_output,*) 'caught by softcut',vtr(ivTot,itrial)*beta
             my_lovr(my_itrial) = .true.
             my_bfac(my_itrial) = 0.0E0_dp
          else if((my_vtry(my_itrial)*beta).lt.-2.303E0_dp*308)then
             ! write(io_output,*) '### warning: weight too big out of range'
             my_lovr(my_itrial) = .true.
             my_bfac(my_itrial) = 0.0E0_dp
          else
             my_bfac(my_itrial) = exp ( -(my_vtry(my_itrial)*beta) )
          end if
       end if
    end do

    call mp_allgather(my_vtry,vtr(ivTot,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrintra,vtr(ivIntraLJ,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrext,vtr(ivExt,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrinter,vtr(ivInterLJ,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrelect,vtr(ivElect,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrewald,vtr(ivEwald,:),rcounts,displs,groupid)
    call mp_allgather(my_bfac,bfac,rcounts,displs,groupid)
    call mp_allgather(my_vipswot,vipswot,rcounts,displs,groupid)
    call mp_allgather(my_vwellipswot,vwellipswot,rcounts,displs,groupid)
    call mp_allgather(my_vipswnt,vipswnt,rcounts,displs,groupid)
    call mp_allgather(my_vwellipswnt,vwellipswnt,rcounts,displs,groupid)
    call mp_allgather(my_lovr,lovr,rcounts,displs,groupid)
    ovrlap = .true.
    if (ANY(.not.lovr(1:ichoi))) ovrlap=.false.
! ----------------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'end BOLTZ in ',myid
#endif
    return
  end subroutine boltz

!> \brief Tail corrections to energy
!>
!> \warning Tabulated potential, MMFF94 and Feuston-Garofalini
!> are NOT supported.
!> \note New potenial functional form needs to be added here,
!> and in U2 (energy calculation), prop_pressure::pressure
!> (force calculation), corp (tail corrections to pressure)
  function coru(imolty,jmolty,rho,ibox)
    real::coru
    integer,intent(in)::imolty,jmolty,ibox
    real,intent(in)::rho
    real::rbcut,rbcut2,rbcut3,tmp,rci1,rci3,rci6,sigma2,epsilon2
    integer::ii,jj,ntii,ntjj,ntij

    rbcut=rcut(ibox)
    rbcut2=rbcut*rbcut
    rbcut3=rbcut2*rbcut

    coru = 0.0_dp
    do ii = 1, nunit(imolty)
       ntii = ntype(imolty,ii)
       do jj = 1, nunit(jmolty)
          ntjj = ntype(jmolty,jj)
          ntij = type_2body(ntii,ntjj)
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
             coru = coru + epsilon2*tmp*rci3*(rci6/9.0_dp-1/3.0_dp)
          else if (nonbond_type(ntij).eq.2) then
             ! Buckingham exp-6
             if (vvdW(2,ntij).ne.0) then
                tmp = vvdW(2,ntij)*rbcut
                coru = coru - (2.0_dp+tmp*(tmp-2.0_dp))*vvdW(1,ntij)*exp(tmp)/(vvdW(2,ntij)**3) - vvdW(3,ntij)/3.0_dp/rbcut3
             end if
          else if (nonbond_type(ntij).eq.3) then
             ! Mie
             tmp=vvdW(2,ntij)/rbcut
             rci1=tmp**(vvdW(3,ntij)-3)
             rci3=tmp**(vvdW(4,ntij)-3)
             coru = coru + vvdW(1,ntij)*(vvdW(2,ntij)**3)*(rci1/(vvdW(3,ntij)-3)-rci3/(vvdW(4,ntij)-3))
          else if (nonbond_type(ntij).eq.5) then
             ! LJ 9-6
             rci3=(vvdW(2,ntij)/rbcut)**3
             coru = coru + vvdW(1,ntij)*(vvdW(2,ntij)**3)*rci3*(rci3/3.0_dp-1.0_dp)
          else if (nonbond_type(ntij).eq.6) then
             ! Generalized LJ
             tmp  = 2.0_dp*vvdW(4,ntij)-3.0_dp
             rci1 = (vvdW(2,ntij)/rbcut)**tmp
             rci3 = (vvdW(2,ntij)/rbcut)**(vvdW(4,ntij)-3.0_dp)
             coru = coru + vvdW(1,ntij)*(vvdW(2,ntij)**3)*(rci1/tmp-rci3*2.0_dp/(vvdW(4,ntij)-3.0_dp))
          else if (nonbond_type(ntij).eq.7) then
             ! LJ 12-6-8
             coru = coru + vvdW(1,ntij)/9.0_dp/(rbcut3**3)-vvdW(2,ntij)/3.0_dp/rbcut3-vvdW(3,ntij)/5.0_dp/rbcut3/rbcut2
          else if (nonbond_type(ntij).eq.10) then
             ! LJ 12-6-8-10
             coru = coru + vvdW(1,ntij)/9.0_dp/(rbcut3**3)-vvdW(2,ntij)/3.0_dp/rbcut3-vvdW(3,ntij)/5.0_dp/rbcut3/rbcut2-vvdW(4,ntij)/(7.0_dp*rbcut3*rbcut3*rbcut)
          else if (ALL(nonbond_type(ntij).ne.(/-1,0,4,8,9/))) then
             call err_exit(__FILE__,__LINE__,'coru: undefined nonbond type',myid+1)
          end if
       end do
    end do

    coru=twopi*rho*coru

  end function coru

!> \brief Read in non-bonded potentials
  subroutine read_ff(io_ff,lmixlb,lmixjo,lmixwh,lmixkong)
    use util_search,only:initiateTable,destroyTable,addToTable,indexOf,tightenTable
    use util_memory,only:reallocate
    use util_mp,only:mp_bcast
    use energy_intramolecular,only:read_ff_bonded
    use energy_sami,only:susami,sumuir
    use energy_garofalini,only:init_garofalini

    integer,INTENT(IN)::io_ff
    logical,intent(in)::lmixlb,lmixjo,lmixwh,lmixkong
    integer,parameter::initial_size=20
    character(LEN=default_string_length)::line_in
    integer::jerr,i,j,ij,ji,nmix,itmp

    !> Looking for section ATOMS
    if (allocated(atoms%list)) then
       call destroyTable(atoms)
       deallocate(atom_type,vvdW_b,qelect,mass,lij,lqchg,chemid,stat=jerr)
    end if
    nntype=0
    if (myid.eq.rootid) then
       rewind(io_ff)
       CYCLE_READ_ATOMS:DO
          call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_atoms

          if (UPPERCASE(line_in(1:5)).eq.'ATOMS') then
             call initiateTable(atoms,initial_size)
             allocate(atom_type(1:initial_size),vvdW_b(1:4,1:initial_size),qelect(1:initial_size),mass(1:initial_size)&
              ,lij(1:initial_size),lqchg(1:initial_size),chemid(1:initial_size),stat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_pairwise: atoms allocation failed',jerr)
             atom_type = 0
             vvdW_b = 0.0_dp
             qelect = 0.0_dp
             mass = 0.0_dp
             lij = .true.
             lqchg = .false.
             do
                call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section ATOMS',myid)
                if (UPPERCASE(line_in(1:9)).eq.'END ATOMS') exit
                nntype=nntype+1
                read(line_in,*) i
                i=addToTable(atoms,i,expand=.true.)
                if (i.gt.ubound(atom_type,1)) then
                   call reallocate(atom_type,1,2*ubound(atom_type,1))
                   call reallocate(vvdW_b,1,4,1,2*ubound(vvdW_b,2))
                   call reallocate(qelect,1,2*ubound(qelect,1))
                   call reallocate(mass,1,2*ubound(mass,1))
                   call reallocate(lij,1,2*ubound(lij,1))
                   call reallocate(lqchg,1,2*ubound(lqchg,1))
                   call reallocate(chemid,1,2*ubound(chemid,1))
                end if
                read(line_in,*) j,atom_type(i),vvdW_b(1:vdW_nParameter(atom_type(i)),i),qelect(i),mass(i),chemid(i)
                if (qelect(i).ne.0) then
                   lqchg(i)=.true.
                else
                   lqchg(i)=.false.
                end if
                if (ALL(vvdW_b(1:vdW_nParameter(atom_type(i)),i).eq.0)) then
                   lij(i)=.false.
                else
                   lij(i)=.true.
                end if
             end do
             exit cycle_read_atoms
          end if
       END DO CYCLE_READ_ATOMS
    end if

    call mp_bcast(nntype,1,rootid,groupid)

    if (nntype.gt.0) then
       if (myid.eq.rootid) then
          call tightenTable(atoms)
          call reallocate(atom_type,1,nntype)
          call reallocate(vvdW_b,1,4,1,nntype)
          call reallocate(qelect,1,nntype)
          call reallocate(mass,1,nntype)
          call reallocate(lij,1,nntype)
          call reallocate(lqchg,1,nntype)
          call reallocate(chemid,1,nntype)
       else
          call initiateTable(atoms,nntype)
          allocate(atom_type(1:nntype),vvdW_b(1:4,1:nntype),qelect(1:nntype),mass(1:nntype),lij(1:nntype),lqchg(1:nntype)&
           ,chemid(1:nntype),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_pairwise: atoms allocation failed',myid)
       end if

       call mp_bcast(atoms%size,1,rootid,groupid)
       call mp_bcast(atoms%list,atoms%size,rootid,groupid)
       call mp_bcast(atom_type,nntype,rootid,groupid)
       call mp_bcast(vvdW_b,4*nntype,rootid,groupid)
       call mp_bcast(qelect,nntype,rootid,groupid)
       call mp_bcast(mass,nntype,rootid,groupid)
       call mp_bcast(lij,nntype,rootid,groupid)
       call mp_bcast(lqchg,nntype,rootid,groupid)
       call mp_bcast(chemid,rootid,groupid)

       nmix=nntype*nntype

       if (allocated(lpl)) deallocate(lpl,xiq,jayself,nonbond_type,vvdW,rminee,ecut,jayq,stat=jerr)
       allocate(lpl(1:nntype),xiq(1:nntype),jayself(1:nntype),nonbond_type(1:nmix),vvdW(1:4,1:nmix),rminee(1:nmix),ecut(1:nmix)&
        ,jayq(1:nmix),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_pairwise: nonbond allocation failed',jerr)

       lpl = .false.
       xiq = 0.0_dp
       nonbond_type=0
       vvdW=0.0_dp
       jayq=0.0_dp
    end if

    ! Computation of un-like interactions
    ! convert input data to program units
    ! calculate square sigmas and epsilons for lj-energy subroutines
    do i = 1, nntype
       do j = 1, nntype
          if (atom_type(i).ne.atom_type(j)) cycle

          ij = type_2body(i,j)
          nonbond_type(ij)=atom_type(i)

          if (ANY(nonbond_type(ij).eq.(/1,3,4,5,6,8/))) then
             ! LJ 12-6 or Mie or MMFF94 or LJ 9-6 or Generalized LJ
             if (lmixlb) then
                ! Lorentz-Berthelot rules --- sig_ij = 0.5 [ sig_i + sig_j ]
                vvdW(2:4,ij) = 0.5_dp*(vvdW_b(2:4,i)+vvdW_b(2:4,j))
             else if (lmixwh) then
                ! Waldman-Hagler rules --- sig_ij = [0.5 [sig_i^6 + sig_j^6] ]^(1/6)
                if (vvdW_b(2,i).eq.0.0_dp.or.vvdW_b(2,j).eq.0.0_dp) then
                   vvdW(2:3,ij)=0.0_dp
                else
                   vvdW(2,ij)=(0.5_dp*(vvdW_b(2,i)**6+vvdW_b(2,j)**6))**(1/6.0_dp)
                   vvdW(3,ij)=vvdW(2,ij)**2
                   vvdW(1,ij)=4.0_dp*sqrt(vvdW_b(1,i)*vvdW_b(1,j)) &
                            /(0.5_dp*(vvdW_b(2,i)**6+vvdW_b(2,j)**6))&
                            *vvdW_b(2,i)**3*vvdW_b(2,j)**3
                end if
             else if (lmixkong) then
                ! Kong rules --- sig_ij = f(sig_ij, eps_ij)
                if (vvdW_b(2,i).eq.0.0_dp.or.vvdW_b(2,j).eq.0.0_dp) then
                   vvdW(2:3,ij)=0.0_dp
                else
                   vvdW(2,ij)=(0.5_dp*((vvdW_b(1,i)*vvdW_b(2,i)**12)**(1/13.0_dp)&
                            +(vvdW_b(1,j)*vvdW_b(2,j)**12)**(1/13.0_dp)))**(13/6.0_dp)&
                            /(vvdW_b(1,i)*vvdW_b(2,i)**6*vvdW_b(1,j)*vvdW_b(2,j)**6)**(1/12.0_dp)
                   vvdW(3,ij)=vvdW(2,ij)**2
                   vvdW(1,ij)=4.0_dp*vvdW_b(1,i)*vvdW_b(2,i)**6*vvdW_b(1,j)*vvdW_b(2,j)**6&
                            /(0.5_dp*((vvdW_b(1,i)*vvdW_b(2,i)**12)**(1/13.0_dp)&
                            +(vvdW_b(1,j)*vvdW_b(2,j)**12)**(1/13.0_dp)))**13
                end if
             else if (lmixjo) then
                ! Jorgensen mixing rules --- sig_ij = [ sig_i * sig_j ]^(1/2)
                vvdW(2:4,ij) = sqrt(vvdW_b(2:4,i)*vvdW_b(2:4,j))
             end if

             if (nonbond_type(ij).eq.3) then
                vvdW(1,ij) = vvdW(3,ij)/(vvdW(3,ij)-vvdW(4,ij))*((vvdW(3,ij)/vvdW(4,ij))**(vvdW(4,ij)/(vvdW(3,ij)-vvdW(4,ij))))&
                 *sqrt(vvdW_b(1,i)*vvdW_b(1,j))
             else if (nonbond_type(ij).eq.8) then
                vvdW(1,ij) = sqrt(vvdW_b(1,i)*vvdW_b(1,j))
                vvdW(3,ij) = vvdW(2,ij)*vvdW(2,ij)
             else
                if (lmixlb .or. lmixjo) then
                   vvdW(1,ij) = 4.0_dp*sqrt(vvdW_b(1,i)*vvdW_b(1,j))
                   if (nonbond_type(ij).eq.1.or.nonbond_type(ij).eq.4) vvdW(3,ij) = vvdW(2,ij)*vvdW(2,ij)
                   if ((nonbond_type(ij).eq.1).and.(vvdW_b(2,i).eq.0.0_dp.or.vvdW_b(2,j).eq.0.0_dp)) vvdW(2:3,ij) = 0.0_dp
                end if
             end if
          else if (nonbond_type(ij).eq.9) then
             ! Hard-core square-well
             vvdW(1,ij) = sqrt(vvdW_b(1,i)*vvdW_b(1,j))
             if (lmixlb) then
                ! Lorentz-Berthelot rules --- sig_ij = 0.5 [ sig_i + sig_j ]
                vvdW(2,ij) = 0.5_dp*(vvdW_b(2,i)+vvdW_b(2,j))
                vvdW(3,ij) = 0.5_dp*(vvdW_b(2,i)*vvdW_b(3,i)+vvdW_b(2,j)*vvdW_b(3,j))
             else if (lmixjo) then
                ! Jorgensen mixing rules --- sig_ij = [ sig_i * sig_j ]^(1/2)
                vvdW(2:3,ij) = sqrt(vvdW_b(2:3,i)*vvdW_b(2:3,j))
                vvdW(3,ij) = vvdW(2,ij)*vvdW(3,ij)
             end if
          else if (nonbond_type(ij).eq.2) then
             ! Buckingham exp-6
             vvdW(1,ij)=sqrt(vvdW_b(1,i)*vvdW_b(1,j))
             vvdW(2,ij)=-0.5_dp*(vvdW_b(2,i)+vvdW_b(2,j))
             vvdW(3,ij)=sqrt(vvdW_b(3,i)*vvdW_b(3,j))
          else if (nonbond_type(ij).eq.7) then
             ! LJ 12-6-8
             vvdW(1:3,ij)=sqrt(vvdW_b(1:3,i)*vvdW_b(1:3,j))
          else if (nonbond_type(ij).eq.10) then
             ! LJ 12-6-8-10
             vvdW(1:4,ij)=sqrt(vvdW_b(1:4,i)*vvdW_b(1:4,j))
          end if
       end do
    end do

    !> Looking for section NONBOND
    if (myid.eq.rootid) then
       REWIND(io_ff)
       CYCLE_READ_NONBOND:DO
          call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_nonbond
          if (UPPERCASE(line_in(1:7)).eq.'NONBOND') then
             jerr=0
             exit cycle_read_nonbond
          end if
       END DO CYCLE_READ_NONBOND
    end if

    call mp_bcast(jerr,1,rootid,groupid)

    if (jerr.eq.0) then
       nmix=0
       do
          if (myid.eq.rootid) then
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section NONBOND',jerr)
          end if

          call mp_bcast(line_in,rootid,groupid)

          if (UPPERCASE(line_in(1:11)).eq.'END NONBOND') exit
          nmix=nmix+1
          read(line_in,*) i,j
          i=indexOf(atoms,i)
          j=indexOf(atoms,j)
          if (i.eq.0.or.j.eq.0) call err_exit(__FILE__,__LINE__,'read_ff: undefined bead in section NONBOND',jerr)
          ij=type_2body(i,j)
          read(line_in,*) itmp,itmp,nonbond_type(ij),vvdW(1:vdW_nParameter(nonbond_type(ij)),ij)

          if (nonbond_type(ij).eq.1) then
             ! LJ 12-6
             vvdW(1,ij)=4.0_dp*vvdW(1,ij)
             vvdW(3,ij)=vvdW(2,ij)**2
          else if (nonbond_type(ij).eq.2) then
             ! Buckingham exp-6
             vvdW(2,ij)=-vvdW(2,ij)
          else if (nonbond_type(ij).eq.3) then
             ! Mie
             vvdW(1,ij)=vvdW(3,ij)/(vvdW(3,ij)-vvdW(4,ij))*((vvdW(3,ij)/vvdW(4,ij))**(vvdW(4,ij)/(vvdW(3,ij)-vvdW(4,ij))))&
              *vvdW(1,ij)
          else if (nonbond_type(ij).eq.4) then
             ! MMFF94
             vvdW(3,ij)=vvdW(2,ij)**2
          else if (nonbond_type(ij).eq.5) then
             ! LJ 9-6
             vvdW(1,ij)=4.0_dp*vvdW(1,ij)
          else if (nonbond_type(ij).eq.6) then
             ! Generalized LJ
             vvdW(1,ij)=4.0_dp*vvdW(1,ij)
          else if (nonbond_type(ij).eq.8) then
             ! DPD potential
             vvdW(3,ij)=vvdW(2,ij)**2
          else if (nonbond_type(ij).eq.9) then
             ! Hard-core square-well
             vvdW(3,ij)=vvdW(2,ij)*vvdW(3,ij)
          else if (ALL(nonbond_type(ij).ne.(/-1,0,7/))) then
             call err_exit(__FILE__,__LINE__,'read_ff: undefined nonbond type',myid+1)
          end if

          ji=type_2body(j,i)
          nonbond_type(ji)=nonbond_type(ij)
          vvdW(:,ji)=vvdW(:,ij)

          if (ANY(vvdW(1:vdw_nParameter(nonbond_type(ij)),ij).ne.0)) then
             lij(i)=.true.
             lij(j)=.true.
          end if
       end do
    end if

    call read_tabulated_ff_pair()

    call read_ff_bonded(io_ff)

    if (lgaro) then
       call init_garofalini()
    else
       if (lsami) then
          call susami()
       end if
       if (lmuir) then
          call sumuir()
       end if
    end if
  end subroutine read_ff

  subroutine init_ff(io_input,lprint)
    use energy_external,only:init_energy_external
    INTEGER,INTENT(IN)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::i,j,ij
    real::rbcut,rbcutsq

    !> read external potentials
    call init_energy_external(io_input,lprint)

    if (lshift) then
       ecut=0.0_dp
       rbcut=rcut(1)
       rbcutsq=rbcut*rbcut
       do i = 1, nntype
          do j = i, nntype
             ij = type_2body(i,j)
             ecut(ij) = U2(rbcutsq,1,1,1,i,2,1,1,j,ij)
             ecut(type_2body(j,i)) = ecut(ij)
          end do
       end do
    end if

! ! --- TraPPE-UA? Methane [CH4] sp3 charged with polarizability
! sigi(28) = 3.73E0_dp
! epsi(28) = 148.0E0_dp
! ! is this correct?
! mass(28) = 16.043E0_dp
! qelect(28) = -0.572E0_dp
! lqchg(28) = .true.
! jayself(28) = 0.5E0_dp*117403E0_dp
! xiq(28) = 9449.3E0_dp
! chname(28) = 'Tr C CH4 chg pol '
! chemid(28)  = 'C  '

! ! --- Methane hydrogen charged with polarizibility
! sigi(29) = 0.0E0_dp
! epsi(29) = 0.0E0_dp
! mass(29) = 1.0078E0_dp
! qelect(29) = 0.143E0_dp
! lqchg(29) = .true.
! jayself(29) = 0.5E0_dp*177700E0_dp
! xiq(29) = 0.0E0_dp
! lij(29) = .false.
! chname(29) = 'Tr H CH4 chg pol '
! chemid(29)  = 'H  '

! ! --- SPC-FQ oxygen [O]   S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(109) = 3.176
! epsi(109) = 148.0E0_dp
! mass(109) = 15.999E0_dp
! qelect(109) = -0.672123708
! lqchg(109) = .true.
! xiq(109) = 36899.0E0_dp
! jayself(109) = (0.5E0_dp)*(503.2E0_dp)*(367.0E0_dp)
! chname(109) = 'SPC-FQ O water   '
! chemid(109)  = '0  '

! ! --- SPC-FQ hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(110) = 0.0E0_dp
! epsi(110) = 0.0E0_dp
! mass(110) = 1.0079E0_dp
! qelect(110) = 0.336061854
! lij(110) = .false.
! lqchg(110) = .true.
! xiq(110) = 0.0E0_dp
! jayself(110) = (0.5E0_dp)*(503.2E0_dp)*(392.2E0_dp)
! chname(110) = 'SPC-FQ H water   '
! chemid(110)  = 'H  '

! ! --- TIP4P-FQ Oxygen [O] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(111) = 3.159E0_dp
! epsi(111) = 144.1E0_dp
! !      epsi(111) = 105.0E0_dp
! mass(111) = 15.999E0_dp
! chname(111) = 'TIP4P-FQ O water '
! chemid(111)  = 'O  '

! ! --- TIP4P-FQ Hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(112) = 0.0E0_dp
! epsi(112) = 0.0E0_dp
! mass(112) = 1.0079E0_dp
! qelect(112) = 0.35E0_dp
! lij(112) = .false.
! lqchg(112) = .true.
! xiq(112) = 0.0E0_dp
! jayself(112) = (0.5E0_dp)*(503.2E0_dp)*(353.0E0_dp)
! chname(112) = 'TIP4P-FQ H water '
! chemid(112)  = 'H  '

! ! --- TIP4P-FQ Charge [Q] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(113) = 0.0E0_dp
! epsi(113) = 0.0E0_dp
! mass(113) = 0.0E0_dp
! qelect(113) = -0.70E0_dp
! lij(113) = .false.
! lqchg(113) = .true.
! xiq(113) = 34464.0E0_dp
! jayself(113) = (0.5E0_dp)*(503.2E0_dp)*(371.6E0_dp)
! chname(113) = 'TIP4P-FQ M water '
! chemid(113)  = 'M  '

! ! --- TraPPE carbon dioxide carbon in [C]O2-fq (jpotoff 2/15/00)
! sigi(131) = 2.80E0_dp
! epsi(131) = 28.5E0_dp
! mass(131) = 12.011E0_dp
! qelect(131) = 0.6512E0_dp
! lqchg(131) = .true.
! !      xiq(131) = (503.2E0_dp)*123.2E0_dp
! xiq(131) = 0.0E0_dp
! jayself(131) = (0.5E0_dp)*(503.2E0_dp)*(233.5E0_dp)
! chname(131) = 'Tr-FQ C in CO2   '
! chemid(131)  = 'C  '

! ! --- TraPPE carbon dioxide oxygen in C[O]2-fq (jpotoff 2/15/00)
! sigi(132) = 3.06E0_dp
! epsi(132) = 80.5E0_dp
! mass(132) = 15.999E0_dp
! qelect(132) = -0.3256E0_dp
! lqchg(132) = .true.
! !      xiq(132) = (503.2E0_dp)*201.56E0_dp
! xiq(132) = 39430.75E0_dp
! jayself(132) = (0.5E0_dp)*(503.2E0_dp)*(308.17E0_dp)
! chname(132) = 'Tr-FQ O in CO2   '
! chemid(132)  = 'O  '

! ! - CO2-FQ Carbon-Oxygen cross term (JCO)
! i = 131
! j = 132
! djay = (503.2E0_dp)*(133.905E0_dp)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ij) = djay
! jayq(ji) = djay

! ! - CO2-FQ Oxygen-Oxygen cross term (JOO)
! i = 132
! j = 132
! djay = (503.2E0_dp)*(1.09E0_dp)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- SPC-FQ water Oxygen-Hydrogen cross term
! i = 109
! j = 110
! djay = (503.2E0_dp)*(276.0E0_dp)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ij) = djay
! jayq(ji) = djay

! ! --- SPC-FQ water Hydrogen-Hydrogen cross term
! i = 110
! j = 110
! djay = (503.2E0_dp)*(196.0E0_dp)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- TIP4P water Charge-Hydrogen cross term
! i = 112
! j = 113
! djay = (503.2E0_dp)*(286.4E0_dp)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ji) = djay
! jayq(ij) = djay

! ! --- TIP4P water Hydrogen-Hydrogen cross term
! i = 112
! j = 112
! djay = (503.2E0_dp)*(203.6E0_dp)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- Methane C-H cross term
! i = 28
! j = 29
! ij = (i-1)*nntype + j
! ji = (j-i)*nntype + i
! jayq(ji) = 114855.0E0_dp
! jayq(ij) = 114855.0E0_dp

! ! --- Methane H-H cross term
! i = 29
! j = 29
! ij = (i-1)*nntype + j
! jayq(ij) = 112537.0E0_dp
  end subroutine init_ff

!> \brief Calculates nonbonding van der Waals potential using linear interpolation between two points.
  function lininter_vdW(r,typi,typj) result(tabulated_vdW)
    use util_math,only:polint
    use util_search,only:LOCATE
    real::tabulated_vdW
    real,intent(in)::r
    integer,intent(in)::typi,typj
    integer::low,high

    low=locate(rvdW(:,typi,typj),vdWsplits(typi,typj),r,2)
    high=low+1
    if (rvdW(low,typi,typj).gt.r.or.rvdW(high,typi,typj).lt.r) then
       write(io_output,*) 'problem in lininter_vdW!'
       write(io_output,*) 'r', r, ' typi', typi, ' typj ', typj
       write(io_output,*) 'low ', low, rvdW(low, typi, typj)
       write(io_output,*) 'high ', high, rvdW(high, typi, typj)
       write(io_output,*)
    end if
    call polint(rvdW(low:high,typi,typj),tabvdW(low:high,typi,typj),2,r,tabulated_vdW)
    return
  end function lininter_vdW

!> \brief Calculates nonbonding electrostatic potential using linear interpolation between two points.
  function lininter_elect(r,typi,typj) result(tabulated_elect)
    use util_math,only:polint
    use util_search,only:LOCATE
    real::tabulated_elect
    real,intent(in)::r
    integer,intent(in)::typi,typj
    integer::low,high

    low=locate(relect(:,typi,typj),electsplits(typi,typj),r,2)
    high=low+1
    if (relect(low,typi,typj).gt.r.or.relect(high,typi,typj).lt.r) then
       write(io_output,*) 'problem in lininter_elect!'
       write(io_output,*) 'r', r, ' typi', typi, ' typj ', typj
       write(io_output,*) 'low ', low, relect(low, typi, typj)
       write(io_output,*) 'high ', high, relect(high, typi, typj)
       write(io_output,*)
    end if
    call polint(relect(low:high,typi,typj),tabelect(low:high,typi,typj),2,r,tabulated_elect)
    return
  end function lininter_elect

!DEC$ ATTRIBUTES FORCEINLINE :: type_2body
  function type_2body(ntii,ntjj)
    use energy_garofalini,only:idx_garofalini
    integer::type_2body
    integer,intent(in)::ntii,ntjj

    if (lgaro) then
       type_2body = idx_garofalini(ntii,ntjj)
    else
       type_2body = (ntii-1)*nntype + ntjj
    end if

  end function type_2body

!> \brief Calculation of pair energy
!>
!> \note New potenial functional form needs to be added here,
!> and in coru (tail corrections to energy), prop_pressure::pressure
!> (force calculation), corp (tail corrections to pressure)
!DEC$ ATTRIBUTES FORCEINLINE :: U2
  function U2(rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
    real::U2
    real,intent(in)::rijsq
    integer,intent(in)::i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij

    real::sr,sr2,sr3,sr6,rs1,rs2,rs6,rs7,rs8,rs10,rs12,tmp,sigma2,epsilon2,qave,rij

    U2 = 0.0_dp
    if (lij(ntii).and.lij(ntjj)) then
       if (lgaro) then
          ! KEA Feuston-Garofalini potential; do NOT work with CBMC/boltz
          rij = sqrt(rijsq)
          U2=garofalini(rij,ntii,ntjj)
       else if (lsami) then
          U2=ljsami(rijsq,ntij)
       else if (lmuir) then
          U2=ljmuir(rijsq,ntij)
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

          if (lfepsi.and.i.ne.j) then ! NOT for intramolecular interactions
             sr6 = rijsq**3
             if ((.not.lqchg(ntii)).and.(.not.lqchg(ntjj))) then
                if (nunit(imolty).eq.4) then
                   !> \bug TIP-4P structure (temperary use?)
                   qave=(qqu(i,4)+qqu(j,4))/2.0_dp
                else
                   qave=(qqu(i,4)+qqu(i,5)+qqu(j,4)+qqu(j,5))*0.85_dp
                end if
             else
                qave=(qqu(i,ii)+qqu(j,jj))/2.0_dp
             end if
             U2=((aslope*(qave-a0)*(qave-a0)+ashift)/sr6-(bslope*(qave-b0)*(qave-b0)+bshift))/sr6*epsilon2
          else
             sr6=sr2**3
             U2=sr6*(sr6-1.0_dp)*epsilon2
          end if
       else if (nonbond_type(ntij).eq.2) then
          ! Buckingham exp-6
          rij = sqrt(rijsq)
          U2 = vvdW(1,ntij)*exp(vvdW(2,ntij)*rij) - vvdW(3,ntij)/(rijsq**3)
       else if (nonbond_type(ntij).eq.3) then
          ! Mie
          rij = sqrt(rijsq)
          sr = vvdW(2,ntij) / rij
          U2 = vvdW(1,ntij)*(sr**vvdW(3,ntij)-sr**vvdW(4,ntij))
       else if (nonbond_type(ntij).eq.4) then
          ! MMFF94
          rij = sqrt(rijsq)
          if (vvdW(2,ntij).ne.0) then
             rs1 = rij/vvdW(2,ntij)
             rs2 = rs1*rs1
             rs7 = rs1*rs2**3
             U2 = vvdW(1,ntij) * ((1.07_dp/(rs1+0.07_dp))**7) * (1.12_dp/(rs7+0.12_dp)-2.0_dp)
          end if
       else if (nonbond_type(ntij).eq.5) then
          ! LJ 9-6
          rij = sqrt(rijsq)
          sr = vvdW(2,ntij)/rij
          sr3 = sr**3
          sr6 = sr3*sr3
          U2 = vvdW(1,ntij)*sr6*(2.0_dp*sr3 - 3.0_dp)
       else if (nonbond_type(ntij).eq.6) then
          ! Generalized LJ
          rij = sqrt(rijsq)
          sr = vvdW(2,ntij) / rij
          if (rij.le.vvdW(2,ntij)) then
             tmp = sr**(vvdW(3,ntij)/2.0_dp)
             U2 = vvdW(1,ntij)*tmp*(tmp-2.0_dp)
          else
             tmp = sr**vvdW(4,ntij)
             U2 = vvdW(1,ntij)*tmp*(tmp-2.0_dp)
          end if
       else if (nonbond_type(ntij).eq.7) then
          ! LJ 12-6-8
          rs6=rijsq**3
          rs8=rs6*rijsq
          rs12=rs6*rs6
          U2=vvdW(1,ntij)/rs12-vvdW(2,ntij)/rs6-vvdW(3,ntij)/rs8
       else if (nonbond_type(ntij).eq.8) then
          ! DPD potential
          rij=sqrt(rijsq)
          if (rij<=vvDW(2,ntij)) U2 = vvdW(1,ntij)*(1.0_dp-rij/vvdW(2,ntij))**2
       else if (nonbond_type(ntij).eq.9) then
          ! Hard-core square-well
          if (rij.lt.vvdW(2,ntij)) then
             U2 = overlapValue
          else if (rij.lt.vvdw(3,ntij)) then
             U2 = -vvdW(1,ntij)
          end if
       else if (nonbond_type(ntij).eq.-1) then
          ! tabulated potential
          rij = sqrt(rijsq)
          U2=lininter_vdW(rij,ntii,ntjj)
       else if (nonbond_type(ntij).eq.10) then
          ! LJ 12-6-8-10
          rs6=rijsq**3
          rs8=rs6*rijsq
          rs10=rs8*rijsq
          rs12=rs6*rs6
          U2=vvdW(1,ntij)/rs12-vvdW(2,ntij)/rs6-vvdW(3,ntij)/rs8-vvdW(4,ntij)/rs10
       else if (nonbond_type(ntij).ne.0) then
          call err_exit(__FILE__,__LINE__,'U2: undefined nonbond type',myid+1)
       end if

       if (lshift) U2 = U2-ecut(ntij)
    end if

    if (i.eq.j) then
       ! intramolecular interactions
       U2=U2*ljscale(imolty,ii,jj)
       if (lainclu(imolty,ii,jj)) then
          ! OH 1-5 interaction
          U2=U2+a15(a15type(imolty,ii,jj))/(rijsq**6)
       end if
    end if

  end function U2

!> \brief Direct-space Coulomb interactions of point charges
!>
!> \note Boltz does NOT call this function. The reason is that
!> when using (neutral-)group-based cutoff, typically and as is
!> done here, charge leaderq is looked up to determine whether
!> two beads interact; however, they may not yet exist during
!> CBMC regrowth.
!DEC$ ATTRIBUTES FORCEINLINE :: Q2
  function Q2(rijsq,rcutsq,i,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo,ibox)
    real::Q2,rij
    real,intent(in)::rijsq,rcutsq,calpi
    integer,intent(in)::i,imolty,ii,ntii,j,jmolty,jj,ntjj,ibox
    logical,intent(in)::lqchgi
    logical,intent(inout)::lcoulo(numax,numax)

    integer::iii,jjj

    Q2=0.0E0_dp
    if(lqchgi.and.lqchg(ntjj)) then
       if (.not.lewald.or.lideal(ibox)) then
          if (.not.lchgall) then
             ! All-Atom charges (charge-group look-up table)
             iii = leaderq(imolty,ii)
             jjj = leaderq(jmolty,jj)
             if (iii.eq.ii .and. jjj.eq.jj)then
                ! set up the charge-interaction table
                if ( rijsq .lt. rcutsq ) then
                   lcoulo(iii,jjj) = .true.
                else
                   lcoulo(iii,jjj) = .false.
                end if
             end if
             ! set up table for neighboring groups- make sure they interact when
             ! leaderqs are only 2 bonds apart. For intramolecular charge interactions
             if (i.eq.j.and..not.lqinclu(imolty,iii,jjj)) then
                lcoulo(iii,jjj)  = .true.
             end if
          end if
          if (lchgall.or.lcoulo(iii,jjj) ) then
             rij = sqrt(rijsq)
             if (L_elect_table) then
                Q2=qqu(i,ii)*qqu(j,jj)*lininter_elect(rij,ntii,ntjj)
             else
                Q2=qqu(i,ii)*qqu(j,jj)/rij
             end if
          end if
       else if (lchgall.or.rijsq.lt.rcutsq) then
          rij = sqrt(rijsq)
          Q2=qqu(i,ii)*qqu(j,jj)*erfunc(calpi*rij)/rij
       end if
    end if
  end function Q2

!> \brief Read in tabulated potential for nonbonded pair interactions (vdW and elect) and set up linear interpolation
  subroutine read_table(file_tab,ntab,r,tab,splits,lists,lused)
    use util_search,only:indexOf
    use util_memory,only:reallocate
    use util_files,only:get_iounit
    use util_mp,only:mp_bcast
    character(len=*),intent(in)::file_tab
    integer,intent(out)::ntab
    integer,allocatable,intent(inout)::splits(:,:)
    real,allocatable,intent(inout)::r(:,:,:),tab(:,:,:)
    type(LookupTable),intent(inout)::lists
    logical,allocatable,intent(inout)::lused(:)

    integer::io_tab,mmm,ii,jj,i,jerr

    splits=0

    if (myid.eq.rootid) then
       io_tab=get_iounit()
       open(unit=io_tab,access='sequential',action='read',file=file_tab,form='formatted',iostat=jerr,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open tabulated potential file: '//file_tab,myid+1)

       read(io_tab,*) ntab
    end if

    call mp_bcast(ntab,1,rootid,groupid)

    do mmm=1,ntab
       ! ii and jj are bead types
       if (myid.eq.rootid) then
          read(io_tab,*) ii, jj
       end if

       call mp_bcast(ii,1,rootid,groupid)
       call mp_bcast(jj,1,rootid,groupid)

       ii=indexOf(lists,ii)
       jj=indexOf(lists,jj)
       if (ii.eq.0.or.jj.eq.0) call err_exit(__FILE__,__LINE__,'read_table: undefined bead',myid+1)
       lused(ii)=.true.
       lused(jj)=.true.
       i=1
       if (myid.eq.rootid) then
          do
             if (i.gt.size(r,1)) then
                call reallocate(r,1,2*size(r,1),1,size(r,2),1,size(r,3))
                call reallocate(tab,1,2*size(tab,1),1,size(tab,2),1,size(tab,3))
             end if
             read(io_tab,*,end=17) r(i,ii,jj),tab(i,ii,jj)
             if (r(i,ii,jj).eq.1000) exit
             ! write(io_tab+10,*) i,r(i,ii,jj),tab(i,ii,jj)
             i=i+1
          end do
17        splits(ii,jj)=i-1
       end if

       call mp_bcast(splits(ii,jj),1,rootid,groupid)

       if (splits(ii,jj).gt.size(r,1)) then
          call reallocate(r,1,splits(ii,jj),1,size(r,2),1,size(r,3))
          call reallocate(tab,1,splits(ii,jj),1,size(tab,2),1,size(tab,3))
       end if

       call mp_bcast(r(:,ii,jj),splits(ii,jj),rootid,groupid)
       call mp_bcast(tab(:,ii,jj),splits(ii,jj),rootid,groupid)
    end do
    if (myid.eq.rootid) close(io_tab)

    call reallocate(r,1,maxval(splits),1,size(r,2),1,size(r,3))
    call reallocate(tab,1,maxval(splits),1,size(tab,2),1,size(tab,3))
  end subroutine read_table

!> \brief Read in nonbonding van der Waals and electrostatic potential.
!> \since KM 12/03/08 vdW
!> \since KM 04/23/09 Q
  subroutine read_tabulated_ff_pair()
    use sim_system,only:L_elect_table
    integer,parameter::grid_size=1500
    integer::jerr

    if (nntype.le.0) then
       return
    end if

    if (ANY(nonbond_type(:).eq.-1)) then
       if (allocated(vdWsplits)) deallocate(vdWsplits,rvdW,tabvdW,stat=jerr)
       allocate(vdWsplits(1:nntype,1:nntype),rvdW(1:grid_size,1:nntype,1:nntype),tabvdW(1:grid_size,1:nntype,1:nntype),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'read_tabulated_potential_pair: allocation failed for vdW_table',myid+1)
       call read_table('fort.43',ntabvdW,rvdW,tabvdW,vdWsplits,atoms,lij)
    end if

    if (L_elect_table) then
       if (allocated(electsplits)) deallocate(electsplits,relect,tabelect,stat=jerr)
       allocate(electsplits(1:nntype,1:nntype),relect(1:grid_size,1:nntype,1:nntype),tabelect(1:grid_size,1:nntype,1:nntype)&
        ,stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'read_tabulated_potential_pair: allocation failed for elect_table',myid+1)
       call read_table('fort.44',ntabelect,relect,tabelect,electsplits,atoms,lqchg)
    end if
  end subroutine read_tabulated_ff_pair

!> \brief compute inter-molecular interaction using kd-tree in the sumup
  subroutine energy_inter_kd_tree_sumup(ibox, lvol, v, ovrlap)
     integer, intent(in) :: ibox
     logical, intent(in) :: lvol
     real, intent(in out) :: v(nEnergy)
     logical, intent(out) :: ovrlap

     integer :: i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij
     logical :: lqchgi, lcoulo(numax,numax)
     type(tree), pointer :: kd_tree
     real :: coord(3), rbcut, rijsq, calpi, rcutsq
     real, allocatable :: inter_list(:,:)
     integer :: iTree, iInter, inter_list_dim, dist_calc_num, dist_calc_num_temp

     ! if called from regular move, ibox = ibox
     ! if called from volume move, ibox = nbox + 1, i.e., to calculate the energy for the coordinates after the volume change
     if (.not. lvol) then
         kd_tree => mol_tree(ibox)%tree
     else
         iTree = 0
         do i = nbox+1, nbox+2
             if (associated(mol_tree(i)%tree)) then
                 if (mol_tree(i)%tree%box .eq. ibox) then
                     iTree = i
                     exit
                 end if
              end if
          end do
          if (iTree .eq. 0) call err_exit(__FILE__,__LINE__,'Error in update_box_kdtree: iTree not found',myid)
          kd_tree => mol_tree(iTree)%tree
     end if

     ! if not volume move, then update and output the height of the tree
     if (.not. lvol) then
         call update_tree_height(kd_tree)
         if (myid .eq. 0) write(io_output, *) "Height of the tree for box",ibox," is ", kd_tree%height
     end if

     rbcut = rcut(ibox)
     rcutsq = rbcut * rbcut
     calpi = calp(ibox)
     dist_calc_num = 0

     ! the kd-tree parallelize i, which differs from non kd-tree, which parallelizes j
     ! the load balancing is not as perfect as j, but considering usually i >> numprocs
     ! this parallelization is not bad
     do i = 1 + myid, nchain - 1, numprocs
         if (nboxi(i) .eq. ibox) then
             imolty = moltyp(i)
             do ii = 1, nunit(imolty)
                 coord(1) = rxu(i, ii)
                 coord(2) = ryu(i, ii)
                 coord(3) = rzu(i, ii)
                 ntii = ntype(imolty, ii)
                 lqchgi = lqchg(ntii)

                 call range_search(kd_tree, coord, nmax*numax, rmin, rbcut, i, ii, .true., ovrlap, inter_list_dim,&
                       inter_list, dist_calc_num_temp, .false.)

                 dist_calc_num = dist_calc_num + dist_calc_num_temp

                 if (ovrlap) then
                     if ((.not. lvol) .and. (myid .eq. rootid)) write(io_output, *) 'Overlap in sumup'
                     goto 399
                 else
                     do iInter = 1, inter_list_dim
                         rijsq = inter_list(1, iInter)
                         j = inter_list(2, iInter)
                         jj = inter_list(3, iInter)
                         jmolty = moltyp(j)
                         ntjj = ntype(jmolty, jj)
                         ntij = type_2body(ntii, ntjj)
                         v(ivInterLJ)=v(ivInterLJ)+U2(rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                         v(ivElect)=v(ivElect)+Q2(rijsq,rcutsq,i,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo,ibox)
                     end do
                 end if
             end do
         end if
     end do !< i

399  continue
     call mp_lor(ovrlap,1,groupid)
     if (ovrlap) return

     call mp_sum(v(ivInterLJ),1,groupid)
     call mp_sum(v(ivElect),1,groupid)
     return

  end subroutine energy_inter_kd_tree_sumup

!> \brief compute inter-molecular interaction using kd-tree in the energy
  subroutine energy_inter_kd_tree_energy(i, imolty, v, flagon, ibox, istart, iuend, ovrlap, lcoulo)
     integer, intent(in) :: ibox
     real, intent(in out) :: v(nEnergy)
     logical, intent(out) :: ovrlap

     integer :: i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij, flagon, istart, iuend, nchp2
     logical :: lqchgi, lcoulo(numax,numax), liji
     type(tree), pointer :: kd_tree
     real :: coord(3), rbcut, rijsq, calpi, rcutsq
     real, allocatable :: inter_list(:,:)
     integer :: iTree, iInter, inter_list_dim, dist_calc_num, dist_calc_num_temp

     kd_tree => mol_tree(ibox)%tree
     dist_calc_num = 0
     nchp2 = nchain + 2
     rbcut = rcut(ibox)
     calpi = calp(ibox)
     rcutsq = rbcut * rbcut

     ! loop over all units in chain i
     ! energy calculation using kd-tree does not loop over j
     ! so here it parallelizes the number of beads whose positions are changed
     ! this results in poor parallel efficiency if the number of beads changed in a move is
     ! smaller than the # of cores (e.g. > 2 cores for ethane, > 1 core for LJ)
     ! but it works well for most simulations, where # of cores > # of beads per molecule
     do ii = istart + myid, iuend, numprocs
         ntii = ntype(imolty,ii)
         liji = lij(ntii)
         lqchgi = lqchg(ntii)
         coord(1) = rxuion(ii, flagon)
         coord(2) = ryuion(ii, flagon)
         coord(3) = rzuion(ii, flagon)

         call range_search(kd_tree, coord, nmax*numax, rmin, rbcut, i, ii, .false., ovrlap, inter_list_dim,&
              inter_list, dist_calc_num_temp, .false.)
         dist_calc_num = dist_calc_num + dist_calc_num_temp
         if (ovrlap) goto 499

         ! After all the interaction sites within the cutoff are found
         ! Loop over all these sites and calculate the inter LJ as well as the electrostatics
         do iInter = 1, inter_list_dim
             rijsq = inter_list(1, iInter)
             j = inter_list(2, iInter)
             jj = inter_list(3, iInter)
             jmolty = moltyp(j)
             if (lexclu(imolty,ii,jmolty,jj)) cycle !< check exclusion table
             ntjj = ntype(jmolty, jj)
             ntij = type_2body(ntii, ntjj)

             v(ivInterLJ)=v(ivInterLJ)+U2(rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)

             ! In calculating electrostatics, nchp2 rather than i is used
             ! because the value here is only needed for determing the charge of this bead qqu(nchp2,ii)
             ! in the swatch move, qqu(nchp2,ii) may be different from qqu(i,ii), because i is the mol to be swatched
             ! and nchp2 is the mol to replace, so qqu(nchp2,ii) should be used here
             ! in all other moves, qqu(nchp2,ii) == qqu(i,ii)
             v(ivElect)=v(ivElect)+Q2(rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo,ibox)
         end do
     end do

499  continue

    call mp_lor(ovrlap,1,groupid)
    if (ovrlap) return

    call mp_sum(v(ivInterLJ),1,groupid)
    call mp_sum(v(ivElect),1,groupid)
    return

  end subroutine energy_inter_kd_tree_energy

!> \brief compute inter-molecular interaction using kd-tree in the boltz
  subroutine energy_inter_kd_tree_boltz(i, icharge, imolty, v, ibox, ntogrow, ovrlap, glist, itrial)
     integer, intent(in) :: ibox, glist(numax), ntogrow, itrial
     real, intent(in out) :: v(nEnergy)
     logical, intent(out) :: ovrlap

     integer :: i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij, icharge, count
     logical :: lqchgi, lcoulo(numax,numax), liji
     type(tree), pointer :: kd_tree
     real :: coord(3), rbcut, rijsq, calpi, rij
     real, allocatable :: inter_list(:,:)
     integer :: iTree, iInter, inter_list_dim, dist_calc_num, dist_calc_num_temp

     kd_tree => mol_tree(ibox)%tree
     dist_calc_num = 0

     if (ldual) then
         rbcut = rcutin
     else
         rbcut = rcut(ibox)
     end if

     ! loop over all beads of molecule i grown this step
     do count = 1, ntogrow
         ! assign bead type for ii
         ii = glist(count)
         ntii = ntype(imolty,ii)
         liji = lij(ntii)
         lqchgi = lqchg(ntii)

         ! assign positions to coord
         coord(1) = rxp(count,itrial)
         coord(2) = ryp(count,itrial)
         coord(3) = rzp(count,itrial)

         call range_search(kd_tree, coord, nmax*numax, rmin, rbcut, i, ii, .false., ovrlap, inter_list_dim,&
             inter_list, dist_calc_num_temp, .false.)

         dist_calc_num = dist_calc_num + dist_calc_num_temp

         if (ovrlap) return

         do iInter = 1, inter_list_dim
             rijsq = inter_list(1, iInter)
             j = inter_list(2, iInter)
             jj = inter_list(3, iInter)
             jmolty = moltyp(j)
             if (lexclu(imolty,ii,jmolty,jj)) cycle !< check exclusion table
             ntjj = ntype(jmolty, jj)
             ntij = type_2body(ntii, ntjj)
             v(ivInterLJ)=v(ivInterLJ)+U2(rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)

             if (L_Coul_CBMC.and.lqchg(ntii).and.lqchg(ntjj)) then
                 rij = sqrt(rijsq)
                 if (L_elect_table) then
                     v(ivElect) = v(ivElect) + qqu(icharge,ii)*qqu(j,jj)*lininter_elect(rij,ntii,ntjj)
                 else if (lewald) then
                     ! compute real space term of velect
                     v(ivElect) = v(ivElect) + qqu(icharge,ii)*qqu(j,jj)*erfunc(calp(ibox)*rij)/rij
                 else
                     ! compute all electrostatic interactions
                     v(ivElect) = v(ivElect) + qqu(icharge,ii)*qqu(j,jj)/ rij
                 end if
             end if
         end do
     end do

     return
  end subroutine energy_inter_kd_tree_boltz

!> \brief compute inter-molecular interaction using COM kd-tree in the boltz subroutine
!> \brief return the lij_list
!> ibox: the box of interest
!> i: the molecule of interest
!> r*ui: coordinates of interest
!> rbcut: rcut used (full rcut or rcutin)
!> rbmin: rmin used in the range_search
!> lfirst: whether the first bead in the swap move
!> lij_list: a quasi-1d list, j-th element indicates whether j mlcl has interaction with mlcl i
!> lsumup: whether it is the sumup call
!> ovrlap: whether there is overlap
  subroutine energy_inter_com_kd_tree(ibox, i, rxui, ryui, rzui, rbcut, rbmin, lfirst, lij_list, lsumup, ovrlap)
     integer :: ibox, i
     type(tree), pointer :: kd_tree
     real :: coord(3), rbcut_plus_buffer, rxui, ryui, rzui, rbcut, rbmin
     real, allocatable :: lij_list(:,:)
     integer::inter_list_dim, dist_calc_num_temp
     logical :: lfirst, ovrlap, lsumup

     kd_tree => mol_tree(ibox)%tree

     if (lfirst) then
         rbcut_plus_buffer = rbcut + kdtree_buffer_len(ibox)
     else
         rbcut_plus_buffer = rbcut + 2.0 * kdtree_buffer_len(ibox)
     end if

     ! do the range_search in kd-tree and get back the "logical" array
     ! indicating what molecules are within the COM cutoff (cutoff+rmsq)
     coord(1) = rxui
     coord(2) = ryui
     coord(3) = rzui
     call range_search(kd_tree, coord, nmax, rbmin, rbcut_plus_buffer, i, 1, lsumup, ovrlap, inter_list_dim,&
                       lij_list, dist_calc_num_temp, .false.)

     return
  end subroutine energy_inter_com_kd_tree
end MODULE energy_pairwise
