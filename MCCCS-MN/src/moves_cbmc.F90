MODULE moves_cbmc
  use var_type,only:dp
  use const_math,only:onepi,twopi
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use energy_kspace,only:recip
  use energy_pairwise,only:energy,boltz
  use energy_intramolecular,only:lininter_vib,lininter_bend,vtorso
  implicit none
  private
  save
  public::config,rosenbluth,schedule,explct,safeschedule,allocate_cbmc,init_cbmc,opt_safecbmc,output_cbmc_stats&
   ,output_safecbmc,read_checkpoint_cbmc,write_checkpoint_cbmc,align_lines,dihedral_rigrot

  ! CBMC.INC
  logical,public::llrig
  logical,allocatable,public::lexshed(:)& !< true if the bead exists at that time of the growth, false when a bead has not yet been grown this time
   ,llplace(:),lpnow(:)
  logical,allocatable::lsave(:)
  integer,allocatable,public::first_bead_to_swap(:)
  real,allocatable::bncb(:,:,:),bscb(:,:,:,:),fbncb(:,:,:),fbscb(:,:,:,:) !< temporary accumulators for conf.bias performance
  real::brvibmin(60),brvibmax(60)

  ! FIX.INC
  integer::endnum,wbefnum,nplace,nrigi,counttot,counthist
  real::hist(60,60,maxbin),probf(60,60,maxbin)
  integer,allocatable::iend(:),ipast(:,:),pastnum(:),fclose(:,:),fcount(:),iwbef(:),ibef(:,:),befnum(:),nextnum(:),inext(:,:)&
   ,rlist(:,:),rfrom(:),rprev(:),rnum(:),iplace(:,:),pfrom(:),pnum(:),pprev(:)
  logical::lcrank
  real,allocatable::xx(:),yy(:),zz(:),distij(:,:),vtvib(:),vtbend(:),vtgtr(:),bsum_tor(:)

  ! subroutine safecbmc
  real,allocatable::kforceb(:,:),equilb(:,:),flength(:,:),vequil(:,:),vkforce(:,:)

contains
!*****************************************************************
!> \brief Performs a length conserving configurational bias move
!> for linear, branched, anisotropic, and explicit atom molecules
!>
!> \b bncb(i,ibox,inb): number of trial attempts starting at unit inb for imolty i in box ibox \n
!> \b bscb(i,1,ibox,inb): number of successful generations of trial configuration \n
!> \b bscb(i,2,ibox,inb): number of accepted trial configurations
!>
!> \author rewritten from old config and branch subroutines by M.G. Martin 9-19-97
!*****************************************************************
  subroutine config()
    use sim_particle,only:update_neighbor_list_molecule,ctrmas,update_coord_in_tree
    use sim_cell,only:update_linked_cell

      logical::lterm,ovrlap,ltors,lfixnow

      integer::i,j,k,iii,ibox,iunit,igrow,icbu,islen,imolty,iutry
      integer::istt,iett,nchp1,ic,total,bin,count,findex,iw,grouptype

      real::v(nEnergy),vtornew,delen,deleo,vdum,wplace,wrig
      real::dchain,rchain,wnlog,wolog,wdlog,wratio,rcbmc
      real::vrecipn,vrecipo,cwtorfo,cwtorfn,x,y,z
      real::delta_vn,delta_vo
! ------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start CONFIG in ',myid
#endif
! select a chain at random ***
      vnew=0.
      vold=0.
      delta_vn = 0.0_dp
      delta_vo = 0.0_dp
      rchain  = random(-1)
      do icbu = 1,nmolty
         if ( rchain .lt. pmcbmt(icbu) ) then
            imolty = icbu
            exit
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

      if (temtyp(imolty).eq.0) return

! determine whether to use safe-cbmc or group-cbmc or not ***
      rcbmc = random(-1)
      if (rcbmc.lt.pmfix(imolty)) then
         lfixnow = .true.
         grouptype = 0
      else if (rcbmc.lt.pmgroup(imolty)) then
         lfixnow = .false.
         grouptype = 1
      else
         lfixnow = .false.
         grouptype = 0
      end if

      dchain = real(temtyp(imolty),dp)
      i = int( dchain*random(-1) + 1 )
      i = parall(imolty,i)
      ibox = nboxi(i)
      if ( moltyp(i) .ne. imolty ) call err_exit(__FILE__,__LINE__,'screwup config',myid+1)

! store number of units in iunit and # to be grown in igrow ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)

! store position of trial chain in r x/y/z cbu ***
      do icbu = 1, igrow
         rxnew(icbu) = rxu(i,icbu)
         rynew(icbu) = ryu(i,icbu)
         rznew(icbu) = rzu(i,icbu)
      end do

      if (lfixnow) then
         call safeschedule(igrow,imolty,islen,iutry,findex,1)
      else
         call schedule(igrow,imolty,islen,iutry,0,1,grouptype)
      end if

! determine how many beads are being regrown
      total = 0
      do icbu = 1,igrow
         if ( .not. lexshed(icbu) ) total = total + 1
      end do

      if (lfixnow) then
         fbncb(imolty,ibox,findex-1) = fbncb(imolty,ibox,findex-1) + 1.0E0_dp
      else
         bncb(imolty,ibox,total) = bncb(imolty,ibox,total) + 1.0E0_dp
      end if

! if ( lelect(imolty) ) then
! Call qqcheck to setup the group based qq cutoff
! call qqcheck(i,ibox,rxnew(1),rynew(1),rznew(1))
! end if

! grow new chain conformation
      if (grouptype .eq. 1) then
         call group_cbmc_grow(.true.,lterm,i,imolty,ibox)
      else
         call rosenbluth(.true.,lterm,i,i,imolty,islen,ibox,igrow,vdum,lfixnow,cwtorfn,1)
      end if

! termination of cbmc attempt due to walk termination ---
      if ( lterm ) then
! write(io_output,*) 'termination of growth',i
        return
      end if

      if (llrig) then
         call rigfix(.true.,i,ibox,imolty,lterm,wrig)
         if ( lterm ) return
         weight = weight * wrig
      end if

      if (llplace(imolty).and.lfixnow) then
         call place(.true.,lterm,i,imolty,ibox,islen,wplace)
         if ( lterm ) return
         weight = weight * wplace
      end if

! grow old chain conformation
      if (grouptype .eq. 1) then
         call group_cbmc_grow(.false.,lterm,i,imolty,ibox)
      else
         call rosenbluth(.false.,lterm,i,i,imolty,islen,ibox,igrow,vdum,lfixnow,cwtorfo,1)
      end if

! termination of old walk due to problems generating orientations
      if ( lterm ) then
         write(io_output,*) 'CONFIG:old growth rejected in box',ibox ,' for moltyp',imolty
         return
      end if

      if (llrig) then
         call rigfix(.false.,i,ibox,imolty,lterm,wrig)
         if ( lterm ) then
            write(io_output,*) 'CONFIG: old rigid fix rejected'
            return
         end if
         weiold = weiold * wrig
      end if

      if (llplace(imolty).and.lfixnow) then
         call place(.false.,lterm,i,imolty,ibox,islen,wplace)

         if ( lterm ) then
            write(io_output,*) 'CONFIG: old hydrogen placement rejected'
            return
         end if
         weiold = weiold * wplace
      end if

! -----------------------------------------------------------------------------
! Begin DC-CBMC, Explicit Atom and Ewald-sum Corrections

      if ( ldual .or. lewald .or. iunit .ne. igrow  .or. ((.not. lchgall) .and. lelect(imolty)) ) then
! Put on hydrogens for explicit AA model for calculation of COM
! and assign all of the grown new and old beads to rxuion
! with old = 1, new = 2
         do j=1,igrow
            rxuion(j,1)=rxu(i,j)
            ryuion(j,1)=ryu(i,j)
            rzuion(j,1)=rzu(i,j)
            qquion(j,1)=qqu(i,j)
         end do
         do j = 1,igrow
            rxuion(j,2) = rxnew(j)
            ryuion(j,2) = rynew(j)
            rzuion(j,2) = rznew(j)
            qquion(j,2) = qquion(j,1)
         end do
         nchp1=nchain+1
         nboxi(nchp1) = ibox
         moltyp(nchp1) = imolty
         moltion(1) = imolty
         moltion(2) = imolty

         if ( igrow .ne. iunit ) then
! iii = 1 old conformation
            do j = igrow+1, iunit
               rxuion(j,1) = rxu(i,j)
               ryuion(j,1) = ryu(i,j)
               rzuion(j,1) = rzu(i,j)
               qquion(j,1) = qqu(i,j)
            end do
! iii = 2 new conformation
            do j=1, igrow
               rxu(nchp1,j) = rxnew(j)
               ryu(nchp1,j) = rynew(j)
               rzu(nchp1,j) = rznew(j)
            end do
            call explct(nchp1,vtornew,.false.,.false.)
            do j=igrow+1, iunit
               rxuion(j,2) = rxu(nchp1,j)
               ryuion(j,2) = ryu(nchp1,j)
               rzuion(j,2) = rzu(nchp1,j)
               qquion(j,2) = qquion(j,1)
            end do
         end if
      end if

      if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) .or. (lchgall .and. lewald .and. (.not. ldual))) then
         istt = 1
         iett = igrow

! check new before old
         do iii = 2,1,-1
! calculate the Full rcut Lennard-Jones energy for the grown beads
! iii = 1 old conformation
! iii = 2 new conformation

            call energy(i,imolty,v,iii,ibox,istt,iett,.true.,ovrlap,.false.,.false.,.false.,.false.)

! write(98,*) '------------ ',iii
! write(98,*) v(ivTot),v(ivIntraLJ),v(ivInterLJ),v(ivExt),v(ivElect),v(ivEwald)
! write(98,*) '------------ new ',iii
! write(98,*) vnew(ivIntraLJ),vnew(ivInterLJ),vnew(ivExt),vnew(ivElect),vnew(ivEwald)

            if (ovrlap .and. (iii .eq. 1)) then
! if (ovrlap) then
               write(io_output,*) 'disaster: overlap in old conf config',i
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if

            if (iii .eq. 2) then
               delen = ( vnew(ivInterLJ) + vnew(ivExt) + vnew(ivElect) + vnew(ivEwald) + vnew(ivIntraLJ))
               if (lstagea) then
                  delen = (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*delen
               else if (lstageb) then
                  delen = etais*delen
               else if (lstagec) then
                  delen = (etais+(1.0E0_dp-etais)*lambdais)*delen
               end if
               delen = v(ivTot) - delen
! JLR 11-19-09 Commenting this out, it makes no sense and gives me an energy error!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM   this may not exactly right, but for
!!!!! numerical reasons I think we need it
! IF(delen*beta .LT. -2.3E0_dp*softcut) THEN
! delen=-2.3E0_dp*softcut/beta
! else if(delen*beta .GT. 2.3E0_dp*softcut) THEN
! delen=2.3E0_dp*softcut
! end if
! END JLR 11-19-09
               delta_vn= delen
               vnew(ivTot)     = vnew(ivTot) + delen
               vnew(ivInterLJ) = v(ivInterLJ)
               vnew(ivExt)   = v(ivExt)
               vnew(ivElect) = v(ivElect)
               vnew(ivIntraLJ) = v(ivIntraLJ)
               vnew(ivEwald) = v(ivEwald)
            else
               deleo = ( vold(ivInterLJ) + vold(ivExt) + vold(ivElect) + vold(ivEwald) + vold(ivIntraLJ))
               if (lstagea) then
                  deleo = (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*deleo
               else if (lstageb) then
                  deleo = etais*deleo
               else if (lstagec) then
                  deleo = (etais+(1.0E0_dp-etais)*lambdais)*deleo
               end if
               deleo = v(ivTot) - deleo
! JLR 11-19-09 Commenting this out, it gives me energy error!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM   this may not exactly right, but for
!!!!! numerical reasons I think we need it
! IF(deleo*beta .LT. -2.3E0_dp*softcut) THEN
! deleo=-2.3E0_dp*softcut/beta
! else if(deleo*beta .GT. 2.3E0_dp*softcut) THEN
! deleo=2.3E0_dp*softcut
! end if
! END JLR 11-19-09
               delta_vo = deleo
               vold(ivTot)     = vold(ivTot) + deleo
               vold(ivInterLJ) = v(ivInterLJ)
               vold(ivExt)   = v(ivExt)
               vold(ivElect) = v(ivElect)
               vold(ivIntraLJ) = v(ivIntraLJ)
               vold(ivEwald) = v(ivEwald)
            end if
         end do

      end if

      if ( iunit .ne. igrow ) then
         istt = igrow+1
         iett = iunit

! check new before old
         do iii = 2,1,-1
! calculate the true Lennard-Jones energy for the hydrogens
! iii=1 old conformation
! iii=2 new conformation
! hydrogens were placed and rxuion was assigned above

            if (iii .eq. 1) then
               ltors = .true.
            else
               ltors = .false.
            end if

! Calculate the energy of the non-backbone beads
            call energy(i,imolty,v,iii,ibox,istt,iett,.true.,ovrlap,ltors,.true.,.false.,.false.)

            if (iii .eq. 2) then
               if (ovrlap) return
               delen = v(ivTot) + vtornew
               if ( delen*beta .gt. (2.3E0_dp*softcut) ) then
! write(io_output,*) '##softcut in config caught explicit atoms'
                  return
               end if
               weight = weight*exp(-(beta*delen))
               vnew(ivTot)  = vnew(ivTot) + delen
               vnew(ivIntraLJ) = vnew(ivIntraLJ) + v(ivIntraLJ)
               vnew(ivInterLJ) = vnew(ivInterLJ) + v(ivInterLJ)
               vnew(ivExt)   = vnew(ivExt) + v(ivExt)
               vnew(ivTorsion) = vnew(ivTorsion) + vtornew
               vnew(ivElect) = vnew(ivElect) + v(ivElect)
               vnew(ivEwald) = vnew(ivEwald) + v(ivEwald)
            else
               if (ovrlap) then
                  write(io_output,*) 'ovrlap problem in old confomation', ' - CONFIG'
                  return
               end if
               deleo = v(ivTot) + v(ivTorsion)
               weiold = weiold*exp(-(beta*deleo))
               if ( weiold .lt. softlog ) then
                  write(io_output,*) '##old weight for explicit too low'
               end if
               vold(ivTot)     = vold(ivTot) + deleo
               vold(ivIntraLJ) = vold(ivIntraLJ) + v(ivIntraLJ)
               vold(ivInterLJ) = vold(ivInterLJ) + v(ivInterLJ)
               vold(ivExt)   = vold(ivExt) + v(ivExt)
               vold(ivTorsion)    = vold(ivTorsion) + v(ivTorsion)
               vold(ivElect) = vold(ivElect) + v(ivElect)
               vold(ivEwald) = vold(ivEwald) + v(ivEwald)
            end if
         end do
      end if

      vrecipn = 0.0E0_dp
      vrecipo = 0.0E0_dp
      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
! reciprocal space sum ---
! rxuion: 1= old configuration; 2= new configuration
         call recip(ibox,vrecipn,vrecipo,1)
         delen = vrecipn
         deleo = vrecipo
         vnew(ivElect) = vnew(ivElect) + vrecipn
         vold(ivElect) = vold(ivElect) + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
            vrecipn = (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipn
            vrecipo = (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipo
         else if (lstageb) then
            vrecipn = etais*vrecipn
            vrecipo = etais*vrecipo
         else if (lstagec) then
            vrecipn = (etais+(1.0E0_dp-etais)*lambdais)*vrecipn
            vrecipo = (etais+(1.0E0_dp-etais)*lambdais)*vrecipo
         end if
         vnew(ivTot) = vnew(ivTot) + vrecipn
         vold(ivTot) = vold(ivTot) + vrecipo
      end if

! End of DC-CBMC, Explicit Atom and Ewald-sum Corrections

! check for acceptance of trial configuration ***
      wnlog = log10 ( weight )
      wolog = log10 ( weiold )
! write(io_output,*) 'weight:',weight
! write(io_output,*) 'weiold:',weiold
      wdlog = wnlog - wolog - beta*(delta_vn + vrecipn - delta_vo - vrecipo)/log(10.0_dp)
      if ( wdlog .lt. -softcut ) then
! write(99,*) 'cbmc softcut',i
         return
      end if

      if (lfixnow) then
         wratio = weight * cwtorfo / ( weiold * cwtorfn)
         fbscb(imolty,1,ibox,findex-1) =  fbscb(imolty,1,ibox,findex-1) + 1.0E0_dp
      else
         wratio = weight / weiold
! write(99,*) weight,weiold
         bscb(imolty,1,ibox,total) = bscb(imolty,1,ibox,total) + 1.0E0_dp
      end if

! write(99,*) wratio
      wratio=wratio*exp(beta*(delta_vo+vrecipo-delta_vn-vrecipn))
! write(99,*) wratio,vold,vnew

      if ( random(-1) .le. wratio ) then
! write(io_output,*) 'CONFIG accepted',i,ibox
! we can now accept !!!!! ***
         if (lfixnow) then
            fbscb(imolty,2,ibox,findex-1) = fbscb(imolty,2,ibox,findex-1)  + 1.0E0_dp
         else
            bscb(imolty,2,ibox,total) = bscb(imolty,2,ibox,total) + 1.0E0_dp
         end if

         vbox(ivTot,ibox)    = vbox(ivTot,ibox)    + ( vnew(ivTot) - vold(ivTot) )
         vbox(ivInterLJ,ibox) = vbox(ivInterLJ,ibox) + (vnew(ivInterLJ) - vold(ivInterLJ))
         vbox(ivIntraLJ,ibox) = vbox(ivIntraLJ,ibox) + (vnew(ivIntraLJ)- vold(ivIntraLJ))
         vbox(ivStretching,ibox)   =  vbox(ivStretching,ibox)  + (vnew(ivStretching)- vold(ivStretching))
         vbox(ivTorsion,ibox)    = vbox(ivTorsion,ibox)    + (vnew(ivTorsion)- vold(ivTorsion))
         vbox(ivExt,ibox)   = vbox(ivExt,ibox)   + (vnew(ivExt) - vold(ivExt))
         vbox(ivBending,ibox)  = vbox(ivBending,ibox)  + (vnew(ivBending) - vold(ivBending))
         vbox(ivElect,ibox) = vbox(ivElect,ibox) + (vnew(ivElect) - vold(ivElect)) + (vnew(ivEwald) - vold(ivEwald))
         vbox(ivIpswb,ibox) = vbox(ivIpswb,ibox) + (vipswn-vipswo)
         vbox(ivWellIpswb,ibox) = vbox(ivWellIpswb,ibox) + (vwellipswn-vwellipswo)
         vipsw = vbox(ivIpswb,ibox)
         vwellipsw = vbox(ivWellIpswb,ibox)

         ! Update coordinates in kdtree
         if ((.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
             do ic = 1, iunit
                 rxu_update(ic) = rxnew(ic)
                 ryu_update(ic) = rynew(ic)
                 rzu_update(ic) = rznew(ic)
             end do
             call update_coord_in_tree(i, igrow, ibox, ibox, .true., .false.)
         end if

         ! Update coordinates in r*u arrays
         do ic = 1, igrow
            rxu(i,ic) = rxnew(ic)
            ryu(i,ic) = rynew(ic)
            rzu(i,ic) = rznew(ic)
         end do

         do ic = igrow+1, iunit
            rxu(i,ic)  = rxuion(ic,2)
            ryu(i,ic)  = ryuion(ic,2)
            rzu(i,ic)  = rzuion(ic,2)
         end do

         if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
! update reciprocal-space sum
            call recip(ibox,vdum,vdum,2)
         end if

         if (ldielect) then
            call dipole(ibox,1)
         end if

! update center of mass
         call ctrmas(.false.,ibox,i,7)
! update linkcell, if applicable
         if ( licell .and. (ibox.eq.boxlink)) then
            call update_linked_cell(i)
         end if

         ! update the neighbour map
         if (lneigh.or.lneighbor) then
            call update_neighbor_list_molecule(i)
         end if
       end if

       if (lpresim.or.lfixnow) then
! record bond distances for presimulation and reweighting
          counthist = counthist + 1
          do iw = 1, islen
             do count = 1, grownum(iw)
                k = growlist(iw,count)
                do j = 1, nunit(imolty)
                   if (k.eq.j) cycle
                   x = rxu(i,j) - rxu(i,k)
                   y = ryu(i,j) - ryu(i,k)
                   z = rzu(i,j) - rzu(i,k)

                   bin = anint(10.0E0_dp*sqrt(x**2+y**2+z**2))

                   if (bin.gt.maxbin) cycle

                   hist(j,k,bin) = hist(j,k,bin) + 1
                   hist(k,j,bin) = hist(k,j,bin) + 1
                end do
             end do
          end do
       end if

! -----------------------------------------------------------------
#ifdef __DEBUG__
       write(io_output,*) 'end CONFIG in ',myid,i
#endif
       return
  end subroutine config

!***************************************************************
!> \brief Performs a configurational bias move for branched molecules
!>
!> \param lnew true for new configurations
!> \param lterm true if early terminated
!> \param i perform rosenbluth growth for chain i
!> \param icharge usually same as i
!> \param imolty molecule type of chain i
!> \param ifrom number of grow-from points
!> \param ibox box number of chain i
!> \param igrow number of units to be grown
!> \param wadd rosenbluth weight for rigrot
!> \param lfixnow SAFE-CBMC
!> \param cwtorf rosenbluth weight of the crank-shaft move for the last torsion
!> \param movetype 1 = config moves;\n
!>           2 = swap/swatch moves for flexible molecule
!>           3 = swatch moves for rigid molecules that do not need rigrot but do need regrowth;\n
!>           4 = swatch moves for completely rigid molecule that regrow nothing
!>           5 = swatch moves for rigid molecules when nsampos=2
!>           0 = swatch moves for rigid molecules when nsampos=1
  subroutine rosenbluth(lnew,lterm,i,icharge,imolty,ifrom,ibox,igrow,wadd,lfixnow,cwtorf,movetype,first_bead,second_bead)
    use util_random,only:sphere
    use util_mp,only:mp_set_displs,mp_allgather

    ! variables passed to the subroutine
    logical::lnew,lterm,lwbef
    integer::i,j,ja,icharge,imolty,ifrom,ibox,igrow
    integer,optional::first_bead,second_bead

    ! local variables
    logical::ovrlap,ltorsion,lfixnow,lfixed,lreturn

    integer::glist(numax),iuprev,iufrom,ichoi,ntogrow,count
    integer::iu,iv,iw,ju,ip,ichtor,it,jut2,jut3,jut4,iwalk
    integer::angstart
    real::dum,xub,yub,zub,length,lengtha ,lengthb,wadd
    real::vdha,x,y,z,maxlen,vtorf,rbf,bsum,bs
    real::vbbtr,vvibtr,wei_vib,wbendv,dist
    real::bondlen(numax),bendang(numax),phi(numax),phidisp
    real::cwtorf,vphi

    ! new stuff
    integer::itor,bin,counta,movetype,ku
    real::bf_tor(nchtor_max),vtorsion(nchtor_max),phitors(nchtor_max),ctorf(nchtor_max),vfbbtr(nchtor_max),ctorf_acc(nchmax)&
     ,vfbbtr_acc(nchmax),ran_tor,wei_bend,jacobian

    ! MPI
    integer::rcounts(numprocs),displs(numprocs),my_start,my_end,blocksize,my_itrial,rid
    real::my_bf_tor(nchtor_max),my_vtorsion(nchtor_max),my_phitors(nchtor_max),my_ctorf(nchtor_max),my_vfbbtr(nchtor_max)
! ------------------------------------------------------------------

#ifdef __DEBUG__
    write(io_output,*) 'start ROSENBLUTH in ',myid
#endif

    lterm = .false.
    cwtorf = 1.0E0_dp
    wei_vib = 1.0E0_dp

! *******************************************
! Rosenbluth weight of trial conformation *
! *******************************************
    ! initialize conformation energies and weight
    if ( lnew ) then
       ! set the initial weight to unity ***
       weight  = 1.0E0_dp
       ! set total energy of trial configuration to zero ***
       vnew(ivTot) = 0.0E0_dp
       vnew(ivTorsion) = 0.0E0_dp
       vnew(ivBending) = 0.0E0_dp
       vnew(ivStretching) = 0.0E0_dp
       vnew(ivExt) = 0.0E0_dp
       vnew(ivIntraLJ) = 0.0E0_dp
       vnew(ivInterLJ) = 0.0E0_dp
       vnew(ivElect) = 0.0E0_dp
       vnew(ivEwald)= 0.0E0_dp
       vipswn  = 0.0E0_dp
       vwellipswn = 0.0E0_dp
    else
       ! old conformation
       ! set the initial weight of the old configuration to unity ***
       weiold  = 1.0E0_dp
       ! set total energy of trial configuration to zero ***
       vold(ivTot) = 0.0E0_dp
       vold(ivTorsion) = 0.0E0_dp
       vold(ivBending) = 0.0E0_dp
       vold(ivStretching) = 0.0E0_dp
       vold(ivExt) = 0.0E0_dp
       vold(ivIntraLJ) = 0.0E0_dp
       vold(ivInterLJ) = 0.0E0_dp
       vold(ivElect) = 0.0E0_dp
       vold(ivEwald)= 0.0E0_dp
       vipswo  = 0.0E0_dp
       vwellipswo = 0.0E0_dp
    end if

! for rigid molecules
! JLR 11-14-09 modifying for calls from swatch for rigid molecules
!   we don't want to do rigrot for rigid swatch when nsampos .ge. 3
    if (lrigid(imolty).and.movetype.ne.1) then
       wadd = 1.0E0_dp
       if (movetype.eq.2.or.movetype.eq.0.or.movetype.eq.5) then
          call rigrot(lnew,lterm,i,icharge,imolty,ibox,wadd,first_bead,second_bead,movetype)
       end if

       if (rindex(imolty).eq.0) then
          return
       end if

       if (movetype.eq.4) then
          return
       end if

       if (lterm) then
          return
       end if
    end if

    !Swatch and swap the same from here change imovetype to 2
    if (movetype.gt.2) movetype=2
! --- END JLR 11-24-09

    ! set lexist to lexshed
    do iu = 1,igrow
       lexist(iu) = lexshed(iu)
    end do

    ! calculate all bond vectors for lexist
    do iu = 1, igrow
       do iv = 1, invib(imolty,iu)
          ju = ijvib(imolty,iu,iv)
          if ( lexist(iu) .and. lexist(ju) ) then
             if ( lnew ) then
                ! use new coordinates
                xvec(iu,ju) = rxnew(ju) - rxnew(iu)
                yvec(iu,ju) = rynew(ju) - rynew(iu)
                zvec(iu,ju) = rznew(ju) - rznew(iu)
             else
                ! use old coordinates
                xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
                yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
                zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
             end if
             distij(iu,ju) = sqrt( xvec(iu,ju)**2 + yvec(iu,ju)**2 + zvec(iu,ju)**2 )
          end if
       end do
    end do

   ichoi = nchoi(imolty)
   ichtor = nchtor(imolty)

   ! MPI
   if (numprocs.gt.1) then
      rid=myid
   else
      rid=-1
   end if
   blocksize = ichtor/numprocs
   rcounts = blocksize
   blocksize = ichtor - blocksize * numprocs
   if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
   call mp_set_displs(rcounts,displs,blocksize,numprocs)
   my_start = displs(myid+1) + 1
   my_end = my_start + rcounts(myid+1) - 1

! *************************
! loop over trial units *
! *************************
    do iw = 1, ifrom
       ! set vibration and bending energies for this growth to 0.0
       if (llrig.and.lsave(iw)) cycle
       iufrom = growfrom(iw)
       ntogrow = grownum(iw)
       lfixed = .false.
       lwbef = .false.
       if (lfixnow) then
          do count = 1, ntogrow
             iu = growlist(iw,count)
             do j = 1, wbefnum
                do ja = 1, befnum(j)
                   if (iu.eq.ibef(j,ja)) then
                      ! time to do final crankshaft move
                      call safecbmc(3,lnew,i,iw,igrow,imolty,count,x,y,z,vphi,vtorf,wbendv ,lterm,movetype)
                      if (lterm) then
                         return
                      end if
                      lfixed = .true.
                      wei_bend = wbendv
                      if (lcrank) then
                         vvibtr = vphi
                      else
                         vvibtr = 0.0E0_dp
                      end if
                      do counta = 1, ntogrow
                         glist(counta) = growlist(iw,counta)
                      end do
                      ichoi = nchoi(imolty)
                      ! sometimes this loop makes the code skip geometry, which sets this
                      ! maxlen (the maximum bond length that CBMC will try to grow)
                      maxlen=2.0E0_dp
                      goto 250
                   end if
                end do
             end do
          end do
       end if

       ! perform the biased selection of bond angles and get lengths
       call geometry(lnew,iw,i,imolty,angstart,iuprev,glist,bondlen,bendang,phi,vvibtr,vbbtr,maxlen,wei_bend)

       ! for lfixnow check if there are two beads to go
       if (lfixnow) then
          fix_count:do count = 1, ntogrow
             iu = growlist(iw,count)
             do j = 1, wbefnum
                if (iu.eq.iwbef(j)) then
                   ! lets setup for two beads to go
                   call safecbmc(1,lnew,i,iw,igrow,imolty ,count,x,y,z,vphi,vtorf,wbendv ,lterm,movetype)
                   wei_vib = wei_vib * wbendv
                   vvibtr = vvibtr + vphi
                   lwbef = .true.
                   exit fix_count
                end if
             end do
          end do fix_count
       end if

       ! we now have the bond lengths and angles for the grown beads
       ! select nchoi trial positions based only on torsions
       do ip = 1,ichoi
          lreturn = .false.
205       continue
          ! set up the cone based on iuprev (could be grown if no prev)
          if (.not.lreturn.and.growprev(iw).eq.0) then
             ! calculate random vector on the unit sphere for the first bead
             count = 1
             if ( (.not. lnew) .and. ip .eq. 1 ) then
                ! use old molecule position
                iu = growlist(iw,count)
                length = bondlen(count)
                ! compute unit vector to be used in cone and torsion
                x = ( rxu(i,iu) - rxu(i,iufrom) )/length
                y = ( ryu(i,iu) - ryu(i,iufrom) )/length
                z = ( rzu(i,iu) - rzu(i,iufrom) )/length
                ! store this in xx yy zz
                xx(count) = x
                yy(count) = y
                zz(count) = z
             else
                ! choose randomly on the unit sphere
                call sphere(x,y,z,-1)
                xx(count) = x
                yy(count) = y
                zz(count) = z
             end if

             if ( ntogrow .gt. 1 ) then
                ! set up the cone
                xub = -x
                yub = -y
                zub = -z
                call cone(1,xub,yub,zub,dum,dum)
             end if

             if (lrigid(imolty)) then
                ! For a rigid molecule, the part that needs regrowth comes before all rigid beads.
                ! The first flexible bead have growprev = 0 but is connected to rigid part which may have torsion
                growprev(iw)=riutry(imolty,iw)+1
                lreturn = .true.
                goto 205
             end if
             ltorsion = .false.
          else
             ! set up the cone based on iuprev and iufrom
             length = distij(iuprev,iufrom)
             xub = xvec(iuprev,iufrom) / length
             yub = yvec(iuprev,iufrom) / length
             zub = zvec(iuprev,iufrom) / length
             call cone(1,xub,yub,zub,dum,dum)
             if (movetype.eq.2.and.lring(imolty).and.iw.eq.1) then
                ltorsion = .false.
             else
                ltorsion = .true.
             end if
          end if

          ! Begin loop to determine torsional angle
          if ( ltorsion ) then
             ! initialize bsum_tor
             bsum_tor(ip) = 0.0E0_dp

             my_itrial = 0
             do itor=my_start,my_end
                my_itrial = my_itrial + 1
                if ( (.not. lnew) .and. ip .eq. 1 .and. itor .eq. 1) then
                   ! old conformation - set phidisp to 0.0E0_dp
                   phidisp = 0.0E0_dp
                else
                   ! choose a random displacement angle from anglestart
                   ! assign the positions based on angles and lengths above
                   phidisp = twopi*random(rid)
                end if

                do count = angstart,ntogrow
                   call cone(2,x,y,z,bendang(count),phi(count) + phidisp)
                   ! store the unit vectors in xx, yy, zz
                   xx(count) = x
                   yy(count) = y
                   zz(count) = z
                end do

                ! set energies of trial position to zero ---
                vdha = 0.0E0_dp
                if (movetype.eq.2 .and.lring(imolty).and.iw.lt.3) then
                   goto 300
                end if

                ! compute torsion energy for given trial conformation
                do count = 1,ntogrow
                   iu = growlist(iw,count)

                   do it = 1, intor(imolty,iu)
                      jut2 = ijtor2(imolty,iu,it)
                      jut3 = ijtor3(imolty,iu,it)
                      if ( jut2 .eq. iufrom .and.  jut3 .eq. iuprev) then
                         jut4 = ijtor4(imolty,iu,it)

                         ! jut4 must already exist or we made a big mistake
                         if ( .not. lexist(jut4) )  then
                            ! allow regrowth where one torsion may already exist and one may not
                            cycle
                            ! write(io_output,*) 'jut4,jut3,jut2,iu',jut4,jut3,jut2,iu
                            ! call err_exit(__FILE__,__LINE__,'trouble jut4',myid+1)
                         end if
                         vdha = vdha + vtorso(xvec(jut4,jut3),yvec(jut4,jut3),zvec(jut4,jut3),xvec(jut3,jut2),yvec(jut3,jut2)&
                          ,zvec(jut3,jut2),xx(count),yy(count),zz(count),ittor(imolty,iu,it))
                      end if
                   end do
                end do
300             continue

                ! compute boltzmann factor and add it to bsum_tor
                my_bf_tor(my_itrial) = exp ( -vdha * beta )

                ! store vtorsion and phidisp for this trial
                my_vtorsion(my_itrial) = vdha
                my_phitors(my_itrial) = phidisp

                ! for safecbmc add extra weight to assure closure
                if (lfixnow) then
                   my_ctorf(my_itrial) = 1.0E0_dp
                   my_vfbbtr(my_itrial) = 0.0E0_dp

                   do count = 1, ntogrow
                      length = bondlen(count)

                      if (lnew) then
                         x = rxnew(iufrom) + xx(count)*length
                         y = rynew(iufrom) + yy(count)*length
                         z = rznew(iufrom) + zz(count)*length
                      else
                         x = rxu(i,iufrom) + xx(count)*length
                         y = ryu(i,iufrom) + yy(count)*length
                         z = rzu(i,iufrom) + zz(count)*length
                      end if

                      iu = growlist(iw,count)

                      if (movetype.eq.2.and.lnew) then
                         if (lwbef) then

                            ! determine special closing energies
                            call safecbmc(2,lnew,i,iw,igrow,imolty,count,x,y,z,vphi,vtorf,wbendv,lterm,movetype)

                            my_bf_tor(my_itrial) = my_bf_tor(my_itrial) * vtorf * exp( - beta * vphi )
                            my_ctorf(my_itrial) = my_ctorf(my_itrial) * vtorf

                            my_vfbbtr(my_itrial) =  my_vfbbtr(my_itrial) + vphi
                         else
                            do j = 1, fcount(iu)
                               ju = fclose(iu,j)
                               dist = sqrt((x-rxnew(ju))**2 + (y-rynew(ju))**2 + (z-rznew(ju))**2)
                               bin = anint(dist*10.0E0_dp)

                               my_bf_tor(my_itrial) = my_bf_tor(my_itrial) * probf(iu,ju,bin)
                               my_ctorf(my_itrial) = my_ctorf(my_itrial) * probf(iu,ju,bin)

                               if (iw.gt.2) then
                                  do counta = 1, pastnum(ju)
                                     ku = ipast(ju,counta)
                                     if (.not.lplace(imolty,ku)) then
                                        dist = sqrt((x-rxnew(ku))**2 + (y-rynew(ku))**2 + (z-rznew(ku))**2)
                                        bin = anint(dist*10.0E0_dp)

                                        my_bf_tor(my_itrial) = my_bf_tor(my_itrial) * probf(iu,ku,bin)
                                        my_ctorf(my_itrial) = my_ctorf(my_itrial) * probf(iu,ku,bin)
                                     end if
                                  end do
                               end if
                            end do
                         end if
                      else
                         if (lwbef) then
                            ! determine special closing energies
                            call safecbmc(2,lnew,i,iw,igrow,imolty,count,x,y,z,vphi,vtorf ,wbendv,lterm,movetype)

                            my_bf_tor(my_itrial) = my_bf_tor(my_itrial) * vtorf  * exp( - beta * vphi )
                            my_ctorf(my_itrial) = my_ctorf(my_itrial)  * vtorf

                            my_vfbbtr(my_itrial) = my_vfbbtr(my_itrial) + vphi
                         else
                            if (fcount(iu).gt.0) then
                               do j = 1, fcount(iu)
                                  ju = fclose(iu,j)
                                  dist = sqrt((x-rxu(i,ju))**2 + (y-ryu(i,ju))**2 + (z-rzu(i,ju))**2)
                                  bin = anint(dist*10.0E0_dp)

                                  my_bf_tor(my_itrial) = my_bf_tor(my_itrial) * probf(ju,iu,bin)
                                  my_ctorf(my_itrial) = my_ctorf(my_itrial) * probf(iu,ju,bin)

                                  if (pastnum(ju).ne.0) then
                                     do counta = 1, pastnum(ju)
                                        ku = ipast(ju,counta)
                                        if (.not.lplace(imolty,ku)) then
                                           dist = sqrt((x-rxu(i,ku))**2 + (y-ryu(i,ku))**2 + (z-rzu(i,ku))**2)
                                           bin = anint(dist*10.0E0_dp)

                                           my_bf_tor(my_itrial) = my_bf_tor(my_itrial) * probf(iu,ku,bin)
                                           my_ctorf(my_itrial) = my_ctorf(my_itrial) * probf(iu,ku,bin)
                                        end if

                                     end do
                                  end if
                               end do
                            end if
                         end if
                      end if
                   end do
                end if
                bsum_tor(ip) = bsum_tor(ip) + my_bf_tor(my_itrial)
             end do

             if ( lnew .or. ip .ne. 1 ) then
                ! choose one of the trial sites in a biased fashion
                ran_tor = random(rid)*bsum_tor(ip)
                bs = 0.0E0_dp
                do itor = 1,rcounts(myid+1)
                   bs = bs + my_bf_tor(itor)
                   if ( ran_tor .lt. bs ) then
                      ! save torsion energy of this trial position
                      vtgtr(ip) = my_vtorsion(itor)
                      ctorf_acc(ip) = my_ctorf(itor)
                      vfbbtr_acc(ip) = my_vfbbtr(itor)
                      ! assign the phidisp of this trial postion
                      phidisp = my_phitors(itor)
                      ! exit the loop
                      exit
                   end if
                end do
             else
                ! select the old conformation
                vtgtr(ip) = my_vtorsion(1)
                ctorf_acc(ip) = my_ctorf(1)
                vfbbtr_acc(ip) = my_vfbbtr(1)
                phidisp = my_phitors(1)
             end if

             if (numprocs.gt.1) then
                call mp_allgather(bsum_tor(ip),bf_tor,groupid)
                bsum_tor(ip)=sum(bf_tor(1:numprocs))
                call mp_allgather(vtgtr(ip),vtorsion,groupid)
                call mp_allgather(phidisp,phitors,groupid)
                call mp_allgather(ctorf_acc(ip),ctorf,groupid)
                call mp_allgather(vfbbtr_acc(ip),vfbbtr,groupid)

                if ( lnew .or. ip .ne. 1 ) then
                   ! choose one of the trial sites in a biased fashion
                   ran_tor = random(-1)*bsum_tor(ip)
                   bs = 0.0E0_dp
                   do itor = 1,numprocs
                      bs = bs + bf_tor(itor)
                      if ( ran_tor .lt. bs ) then
                         ! save torsion energy of this trial position
                         vtgtr(ip) = vtorsion(itor)
                         ctorf_acc(ip) = ctorf(itor)
                         vfbbtr_acc(ip) = vfbbtr(itor)
                         ! assign the phidisp of this trial postion
                         phidisp = phitors(itor)
                         ! exit the loop
                         exit
                      end if
                   end do
                else
                   ! select the old conformation
                   vtgtr(ip) = vtorsion(1)
                   ctorf_acc(ip) = ctorf(1)
                   vfbbtr_acc(ip) = vfbbtr(1)
                   phidisp = phitors(1)
                end if
             end if

             ! divide bsum by ichtor
             bsum_tor(ip) = bsum_tor(ip) / real(ichtor,dp)
          else
             ! no torsion energy, choose phidisp at random(-1except old)
             if ( (.not. lnew) .and. ip .eq. 1 ) then
                ! old conformation - set phidisp to 0.0E0_dp
                phidisp = 0.0E0_dp
             else
                ! choose a random displacement angle from anglestart
                ! assign the positions based on angles and lengths above
                phidisp = twopi*random(-1)
             end if

             ! set bsum_tor to 1.0E0_dp
             bsum_tor(ip) = 1.0E0_dp

             ! assign the torsional energy a value of 0.0
             vtgtr(ip) = 0.0E0_dp
             ctorf_acc(ip) = 1.0E0_dp
             vfbbtr_acc(ip) = 0.0E0_dp
          end if

          ! for accepted phidisp set up the vectors
          do count = angstart,ntogrow
             call cone(2,x,y,z,bendang(count),phi(count) + phidisp)
             ! store the unit vectors in xx, yy, zz
             xx(count) = x
             yy(count) = y
             zz(count) = z
          end do

          ! accepted coordinates, save them in r*p(trial)
          do count = 1,ntogrow
             length = bondlen(count)
             if ( lnew ) then
                ! use new positions
                rxp(count,ip) = rxnew(iufrom) + xx(count)*length
                ryp(count,ip) = rynew(iufrom) + yy(count)*length
                rzp(count,ip) = rznew(iufrom) + zz(count)*length
             else
                ! use old coordinates
                rxp(count,ip) = rxu(i,iufrom) + xx(count)*length
                ryp(count,ip) = ryu(i,iufrom) + yy(count)*length
                rzp(count,ip) = rzu(i,iufrom) + zz(count)*length
             end if
          end do
       end do

250    continue

       ! now that we have the trial site need to compute non-bonded energy
       call boltz(lnew,.false.,ovrlap,i,icharge,imolty,ibox,ichoi,iufrom,ntogrow,glist,maxlen)
       if (ovrlap) then
          lterm = .true.
          return
       end if
! ---------------------------------------------------------------------

       ! perform the walk according to the availibility of the choices ***
       ! and calculate the correct weight for the trial walk           ***
       bsum = 0.0E0_dp
       do ip = 1, ichoi
          ! include both the torsional and the LJ/qq
          bsum = bsum + bfac(ip)*bsum_tor(ip)
       end do

       if ( lnew ) then
          ! update new rosenbluth weight - include bending weight
          weight = weight * bsum * wei_bend * wei_vib

          if ( weight .lt. softlog ) then
             lterm=.true.
             return
          end if

          ! select one position at random ---
          rbf = bsum * random(-1)
          bs = 0.0E0_dp
          do ip = 1, ichoi
             if ( .not. lovr(ip) ) then
                bs = bs + bfac(ip)*bsum_tor(ip)
                if ( rbf .lt. bs ) then
                   ! select ip position ---
                   iwalk = ip
                   exit
                end if
             end if
          end do
       else
          ! old conformation, update weiold - include wei_bend
          weiold = weiold * bsum * wei_bend * wei_vib
          if (weiold .lt. softlog) write(io_output,*)  '###old weight too low'
       end if

       if (lfixed) then
          ! determine jacobian contribution for crankshaft
          jacobian = 1.0E0_dp
          do count = 1, ntogrow
             iu = growlist(iw,count)
             if (fcount(iu).gt.0) then
                counta = 1
                ju = fclose(iu,counta)
                if (lnew) then
                   if (movetype.eq.2) then
                      x = rxnew(ju) - rxnew(iufrom)
                      y = rynew(ju) - rynew(iufrom)
                      z = rznew(ju) - rznew(iufrom)
                      length = sqrt( x**2 + y**2 + z**2 )

                      x = rxnew(ju) - rxp(count,iwalk)
                      y = rynew(ju) - ryp(count,iwalk)
                      z = rznew(ju) - rzp(count,iwalk)
                      lengtha = sqrt( x**2 + y**2 + z**2 )
                   else
                      x = rxu(i,ju) - rxnew(iufrom)
                      y = ryu(i,ju) - rynew(iufrom)
                      z = rzu(i,ju) - rznew(iufrom)
                      length = sqrt( x**2 + y**2 + z**2 )

                      x = rxu(i,ju) - rxp(count,iwalk)
                      y = ryu(i,ju) - ryp(count,iwalk)
                      z = rzu(i,ju) - rzp(count,iwalk)
                      lengtha = sqrt( x**2 + y**2 + z**2 )
                   end if

                   x = rxp(count,iwalk) - rxnew(iufrom)
                   y = ryp(count,iwalk) - rynew(iufrom)
                   z = rzp(count,iwalk) - rznew(iufrom)
                   lengthb = sqrt( x**2 + y**2 + z**2 )
                else
                   x = rxu(i,ju) - rxu(i,iufrom)
                   y = ryu(i,ju) - ryu(i,iufrom)
                   z = rzu(i,ju) - rzu(i,iufrom)
                   length = sqrt( x**2 + y**2 + z**2 )

                   x = rxu(i,ju) - rxu(i,iu)
                   y = ryu(i,ju) - ryu(i,iu)
                   z = rzu(i,ju) - rzu(i,iu)
                   lengtha = sqrt( x**2 + y**2 + z**2 )

                   x = rxu(i,iu) - rxu(i,iufrom)
                   y = ryu(i,iu) - ryu(i,iufrom)
                   z = rzu(i,iu) - rzu(i,iufrom)
                   lengthb = sqrt( x**2 + y**2 + z**2 )
                end if
                jacobian = jacobian / (length*lengtha*lengthb)
             end if
          end do
          bsum = bsum * jacobian
       end if

       if ( lnew ) then
          if (lfixnow) then
             if (lwbef) then
                vbbtr = vbbtr + vfbbtr_acc(iwalk)
             end if
             if (lfixed) then
                vbbtr = vtbend(iwalk)
                vvibtr = vvibtr + vtvib(iwalk)
             else
                if (.not.(movetype.eq.2.and.lring(imolty).and.iw.eq.1)) then
                   cwtorf = cwtorf * ctorf_acc(iwalk)
                end if
             end if
          end if

! update new trial energies
          vnew(ivTot)     = vnew(ivTot)     + vtr(ivTot,iwalk)   + vtgtr(iwalk) + vvibtr + vbbtr
          vnew(ivStretching)  = vnew(ivStretching)  + vvibtr
          vnew(ivBending)    = vnew(ivBending)    + vbbtr
          vnew(ivTorsion)    = vnew(ivTorsion)    + vtgtr(iwalk)
          vnew(ivExt)   = vnew(ivExt)   + vtr(ivExt,iwalk)
          vnew(ivIntraLJ) = vnew(ivIntraLJ) + vtr(ivIntraLJ,iwalk)
          vnew(ivInterLJ) = vnew(ivInterLJ) + vtr(ivInterLJ,iwalk)
          vnew(ivElect) = vnew(ivElect) + vtr(ivElect,iwalk)
          vnew(ivEwald) = vnew(ivEwald) + vtr(ivEwald,iwalk)
          vipswn = vipswn+vipswnt(iwalk)
          vwellipswn = vwellipswn+vwellipswnt(iwalk)

          ! if (ldebug) then
          !    write(100+myid,*) 'iwalk: ',iwalk,'; vnewt: ',vnew(ivTot),'; vnewbb: ',vnew(ivBending),'; vnewtg: ',vnew(ivTorsion)
          ! end if
       else
          if (lfixnow) then
             if (lfixed) then
                vbbtr = vtbend(1)
                vvibtr = vvibtr + vtvib(1)
             else
                if (.not.(movetype.eq.2.and.lring(imolty).and.iw.eq.1)) then
                   cwtorf = cwtorf * ctorf_acc(1)
                end if
             end if
             if (lwbef) then
                vbbtr = vbbtr + vfbbtr_acc(1)
             end if
          end if

! update old trail energies
          vold(ivTot)     = vold(ivTot)     + vtr(ivTot,1)   + vtgtr(1) + vvibtr + vbbtr
          vold(ivStretching)  = vold(ivStretching)  + vvibtr
          vold(ivBending)    = vold(ivBending)    + vbbtr
          vold(ivTorsion)    = vold(ivTorsion)    + vtgtr(1)
          vold(ivExt)   = vold(ivExt)   + vtr(ivExt,1)
          vold(ivIntraLJ) = vold(ivIntraLJ) + vtr(ivIntraLJ,1)
          vold(ivInterLJ) = vold(ivInterLJ) + vtr(ivInterLJ,1)
          vold(ivElect) = vold(ivElect) + vtr(ivElect,1)
          vold(ivEwald) = vold(ivEwald) + vtr(ivEwald,1)
          vipswo = vipswo+vipswot(1)
          vwellipswo = vwellipswo+vwellipswot(1)

          ! if (ldebug) then
          !    write(100+myid,*) 'iwalk: ',iwalk,'; voldt: ',vold(ivTot),'; voldbb: ',vold(ivBending),'; voldtg: ',vold(ivTorsion)
          ! end if
       end if

       do count = 1,ntogrow
          iu = growlist(iw,count)
          if ( lnew ) then
             ! assign new positions to r*new
             rxnew(iu) = rxp(count,iwalk)
             rynew(iu) = ryp(count,iwalk)
             rznew(iu) = rzp(count,iwalk)
          end if

          ! set lexist(iu) to true so ewald sum computed properly
          lexist(iu) = .true.

          ! store new existing vectors between beads and iufrom
          if ( lnew ) then
             ! use r*new positions
             xvec(iu,iufrom) = rxnew(iufrom) - rxnew(iu)
             yvec(iu,iufrom) = rynew(iufrom) - rynew(iu)
             zvec(iu,iufrom) = rznew(iufrom) - rznew(iu)
          else
             ! use r*u positions
             xvec(iu,iufrom) = rxu(i,iufrom) - rxu(i,iu)
             yvec(iu,iufrom) = ryu(i,iufrom) - ryu(i,iu)
             zvec(iu,iufrom) = rzu(i,iufrom) - rzu(i,iu)
          end if
          distij(iu,iufrom) = sqrt( xvec(iu,iufrom)**2 + yvec(iu,iufrom)**2 + zvec(iu,iufrom)**2 )

          xvec(iufrom,iu) = - xvec(iu,iufrom)
          yvec(iufrom,iu) = - yvec(iu,iufrom)
          zvec(iufrom,iu) = - zvec(iu,iufrom)
          distij(iufrom,iu) = distij(iu,iufrom)

          if (lfixnow) then
             ! we must store new vectors with endpoints
             if (fcount(iu).gt.0) then
                do j = 1, fcount(iu)
                   ju = fclose(iu,j)

                   if (lnew) then
                      xvec(iu,ju) = rxnew(ju) - rxnew(iu)
                      yvec(iu,ju) = rynew(ju) - rynew(iu)
                      zvec(iu,ju) = rznew(ju) - rznew(iu)
                   else
                      xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
                      yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
                      zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
                   end if

                   xvec(ju,iu) = - xvec(iu,ju)
                   yvec(ju,iu) = - yvec(iu,ju)
                   zvec(ju,iu) = - zvec(iu,ju)

                   distij(iu,ju) = sqrt(xvec(iu,ju)**2 + yvec(iu,ju)**2 + zvec(iu,ju)**2)
                   distij(ju,iu) = distij(iu,ju)
                end do
             end if

          end if
       end do
! ********************************
! end of loop over trial units *
! ********************************
    end do

#ifdef __DEBUG__
    write(io_output,*) 'end ROSENBLUTH in ',myid
#endif

! ------------------------------------------------------------------

    return
  end subroutine rosenbluth

!*****************************************************************
!> \brief Computes the growth shedule for CBMC type moves
!>
!> \param movetype 1 = config moves;\n
!> 2 = swap moves for flexible molecules;\n
!> 3 = swatch moves for flexible molecules;\n
!> 4 = swap/swatch moves for partially rigid molecules;\n
!> 5 = swatch moves for molecules with more than 1 cuts
!> \param grouptype 0 = NOT use group-CBMC
!> 1 = group-CBMC scheduler for long chain molecules
!> 2 = group-CBMC scheduler for repeat unit (iutry is needed)
!*****************************************************************
  subroutine schedule(igrow,imolty,index,iutry,iprev,movetype,grouptype)
    logical::lfind(numax)
    integer::random_index
    integer::kickout,icbu,igrow,imolty,iutry,iut,invtry,iu,ju,gcbmc_prev,grouptype
    integer::ibead,count,ivib,idir,movetype,iprev,itry,i,ib2,isegment
    integer::temp_store(numax),temp_count,izz,outer_sites(numax),index,outer_num,outer_prev(numax),iufrom,outer_try
    integer :: gcbmc_imolty,unit_num_local
    real::dbgrow
    logical,parameter::lprint = .false.
! ------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start SCHEDULE in ',myid
#endif

    if (grouptype .eq. 1) then
       gcbmc_imolty = gcbmc_mol_list(imolty)
       unit_num_local = gcbmc_unit_num(gcbmc_imolty)
    end if

    kickout = 0

11  continue

    ! initialize temp_count
    temp_count = 0

    if (movetype.ne.5) then
       ! initialize lexshed so all beads currently exist
       lexshed(1:igrow) = .true.
       index = 0
    end if

    if (movetype.eq.1) then
       ! this part is just for config right now
       if (kickout.eq.500) call err_exit(__FILE__,__LINE__,'kickout is 500',myid+1)

       ! select the first bead to grow from
       if (lrigid(imolty)) then
          ! N.R. Randomly select one of the grow points
          random_index = int(real(rindex(imolty),dp)*random(-1)+1)
          iutry = riutry(imolty,random_index)
       else if (lrig(imolty).and.nrig(imolty).gt.0) then
          dbgrow = random(-1)*real(nrig(imolty),dp) + 1.0_dp
          iutry = irig(imolty,int(dbgrow))
       else if (grouptype .eq. 2) then
          ! group-CBMC of repeat unit, iutry is given
       else if (icbsta(imolty).gt.0) then
          iutry = icbsta(imolty)
       else if (icbsta(imolty).eq.0) then
          dbgrow = real(igrow,dp)
          iutry = int(dbgrow*random(-1)) + 1
       else
          dbgrow = real(igrow+icbsta(imolty)+1,dp)
          iutry = int(dbgrow*random(-1)) - icbsta(imolty)
       end if

       if (lrig(imolty).and.(nrig(imolty).eq.0)) then
          ju = int(random(-1)*real(nrigmax(imolty)-nrigmin(imolty)+1,dp)) + nrigmin(imolty)
       end if

       ! if group-CBMC is used, iutry is selected differently
       if (grouptype .eq. 1) then
          if (random(-1) .le. 0.5) then
              idir = 1
          else
              idir = -1
          end if

          ! choose which segment to grow
          ! for idir .eq. 1, isegment = (1,...,n-1)
          ! for idir .eq. -1, isegment = (2,...,n)
          isegment = int(random(-1) * (unit_num_local - 1)) + 1 !< 1,...,n - 1
          if (idir .eq. -1) then
             isegment = isegment + 1
          end if

          if (idir .eq. -1) then
             iutry = gcbmc_unit_list(gcbmc_imolty, isegment, 1)
             gcbmc_prev = gcbmc_unit_list(gcbmc_imolty, isegment, 2)
          else
             iutry = gcbmc_unit_list(gcbmc_imolty, isegment, 2)
             gcbmc_prev = gcbmc_unit_list(gcbmc_imolty, isegment, 1)
          end if
          growfrom(1) = iutry
       else if (grouptype .eq. 2) then
          idir = icbdir(imolty)
          growfrom(1) = iutry
       else
          idir = icbdir(imolty)
          growfrom(1) = iutry
       end if

       invtry = invib(imolty,iutry)
       index = 1

       if (invtry.eq.0) then
          ! problem, cannot do config move on a 1 bead molecule
          call err_exit(__FILE__,__LINE__,'cannot do CBMC on a one grow unit molecule',myid+1)
       else if (invtry.eq.1) then
          ! regrow entire molecule, check nmaxcbmc
          if (nmaxcbmc(imolty).lt.igrow-1) then
             kickout = kickout + 1
             goto 11
          end if

          ! group-CBMC does not allow that kind of move
          if (grouptype .eq. 1) call err_exit(__FILE__,__LINE__,'something wrong with group-CBMC scheduler',myid+1)

          growprev(1) = 0
          grownum(1) = invtry

          ivib = invtry
          iut = ijvib(imolty,iutry,ivib)
          if (idir.eq.1.and.iut.lt.iutry) then
             ! cannot grow from this bead, try again
             kickout = kickout + 1
             goto 11
          end if
          growlist(1,ivib) = iut
          lexshed(iut) = .false.
       else
          ! at a branch point, decide how many branches to regrow
          ! regrow all of legal branches (check idir )

          ! if group-CBMC of repeat unit, something is wrong
          if (grouptype.eq.2) call err_exit(__FILE__,__LINE__,'something wrong with group-CBMC scheduler',myid+1)

          count = 0
          do ivib = 1,invtry
             iut = ijvib(imolty,iutry,ivib)
             if (grouptype.eq.1 .and. iut.eq.gcbmc_prev) then
                growprev(1) = iut
             else if (idir.eq.1.and.iut.lt.iutry) then
                !> \bug this is the one and only previous nongrown bead. This is potentially a bug; certain numbering scheme will fail
                growprev(1) = iut
             else
                ! grow these branches
                temp_count = temp_count + 1
                temp_store(temp_count) = iut
             end if
          end do

          do izz = temp_count,1,-1
             ! choose grow bead randomly from temp_store
             itry = int(real(izz,dp)*random(-1)) + 1
             iut = temp_store(itry)
             count = count + 1
             growlist(1,count) = iut
             lexshed(iut) = .false.
             ! update temp_store for next iteration
             temp_store(itry) = temp_store(izz)
          end do
          grownum(1) = count

          if (count.eq.invtry) then
             if (random(-1).gt.pmall(imolty)) then
                ! not pmall so select a previous bead
20              ivib = int(random(-1)*real(invtry,dp)) + 1
                iu = growlist(1,ivib)
                if (lrig(imolty).and.nrig(imolty).gt.0) then
                   if (iu.ne.frig(imolty,int(dbgrow))) goto 20
                end if

                ! we should start to determine a random section to keep rigid
                if (lrigid(imolty)) then
                   if (iu.lt.riutry(imolty,1)) goto 20
                end if

                growprev(1) = iu
                ! replace this unit with the last unit in growlist
                growlist(1,ivib) = growlist(1,invtry)
                ! reduce numgrow by 1
                grownum(1) = grownum(1) - 1
                ! add this back into lexshed
                lexshed(iu) = .true.
             else
                ! we regrew all branches, no previous bead
                growprev(1) = 0
             end if
          else if (invtry-count.ne.1) then
             ! problem in logic, should only be one nongrown bead
             write(io_output,*) 'invtry,count',invtry,count
             write(io_output,*) 'igrow,imolty',igrow,imolty
             call err_exit(__FILE__,__LINE__,'logic problem in schedule',myid+1)
          end if
       end if
       ! end the part that is specific for config
    else if (movetype.eq.2.or.movetype.eq.3.or.movetype.eq.5) then
       ! begin the part that is specific for flexible molecules
       ! iutry is the first bead inserted - need to grow its neighbors
       if (iutry.eq.0) then
          ! no beads to be regrown via cbmc
          return
       end if

       growfrom(index+1) = iutry
       growprev(index+1) = iprev
       invtry = invib(imolty,iutry)

       if (invtry.ne.0) then
          ! grow all of the beads (except iprev) connected to bead iutry
          index = index + 1
          do ivib=1,invtry
             iut = ijvib(imolty,iutry,ivib)
             if (iut.ne.iprev) then
                ! grow these branches
                temp_count = temp_count + 1
                temp_store(temp_count) = iut
             end if
          end do

          count = 0
          do izz = temp_count,1,-1
             ! choose grow bead randomly from temp_store
             itry = int(real(izz,dp)*random(-1)) + 1
             iut = temp_store(itry)
             count = count + 1
             growlist(index,count) = iut
             lexshed(iut) = .false.
             ! update temp_store for next iteration
             temp_store(itry) = temp_store(izz)
          end do
          grownum(index) = count
       end if
       ! end the part that is specific for flexible molecules
    else if (movetype.eq.4) then
       ! begin part that is specific for rigid molecules

       !> \bug only the growth point listed last in riutry can have more than 1 beads to grow!
       index = rindex(imolty)
       if (index.ne.0) then
          ! molecule is partially rigid
          do i=1,index
             ! grow all non-rigid beads from the rigid part
             temp_count = 0
             ! for rigid molecules, have rigid vib last
             growfrom(i)= riutry(imolty,i)
             invtry = invib(imolty,growfrom(i)) - 1
             growprev(i)= ijvib(imolty,growfrom(i),invtry+1)
             grownum(i) = invtry

             do ivib = 1, invtry
                iut = ijvib(imolty,growfrom(i),ivib)
                ! grow these branches
                temp_count = temp_count + 1
                temp_store(temp_count) = iut
             end do

             count = 0
             do izz = temp_count,1,-1
                ! choose grow bead randomly from temp_store
                itry = int(real(izz,dp)*random(-1)) + 1
                iut = temp_store(itry)
                count = count + 1
                growlist(i,count) = iut
                lexshed(iut) = .false.

                ! update temp_store for next iteration
                temp_store(itry) = temp_store(izz)
             end do
          end do
       end if
       ! end part that is specific for rigid molecules
    else
       ! non-existent move type
       write(io_output,*) 'schedule movetype ',movetype
       call err_exit(__FILE__,__LINE__,'non-valid move type',myid+1)
    end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! from here on down config, swap, and swatch are the same

! OLD METHOD - not fully random - removed 6-13-98
!     ibead = 0
!     index = 0
! 50  if ( index .lt. islen ) then
!        index = index + 1
!        do icbu = 1,grownum(index)
!           ibead = ibead + 1
!           iu = growlist(index,icbu)
!           invtry = invib(imolty,iu)
!           if ( invtry .ne. 1) then
!              this will be the next growing bead
!              islen = islen + 1
!              growprev(islen) = growfrom(index)
!              growfrom(islen) = iu
!              count = 0
!              do ivib = 1,invtry
!                 iut = ijvib(imolty,iu,ivib)
!                 if ( iut .ne. growprev(islen) ) then
!                    count = count + 1
!                    growlist(islen,count) = iut
!                    lexshed(iut) = .false.
!                 end if
!              end do
!              grownum(islen) = count
!           end if
!        end do
!        goto 50
!     end if

    ! set up list of "outer" beads that have further growth sites
    ! this method implemented 6-13-98
    ibead=0
    if (index.gt.0) then
       outer_num = 0
       iufrom = growfrom(index)
       do icbu = 1,grownum(index)
          ! increment counter ibead for total number of beads grown
          ibead = ibead + 1

          ! determine whether this bead has any non-grown neighbors
          iu = growlist(index,icbu)
          invtry = invib(imolty,iu)

          if ( invtry .gt. 1 ) then
             ! add one to the stack of outer_sites
             outer_num = outer_num + 1
             outer_sites(outer_num) = iu
             outer_prev(outer_num) = iufrom
          end if
       end do

       ! begin while loop to grow all outer beads until done
70     if ( outer_num .gt. 0 ) then
          ! choose one site randomly from the stack
          outer_try = int(real(outer_num,dp)*random(-1)) + 1

          ! increment index to show this is the next growfrom
          index = index + 1

          iu = outer_sites(outer_try)
          iufrom = outer_prev(outer_try)
          ! assign growfrom and growprev for this index
          growfrom(index) = iu
          growprev(index) = iufrom

          invtry = invib(imolty,iu)
          ! assign the grow beads in random order
          temp_count = 0
          do ivib = 1,invtry
             iut = ijvib(imolty,iu,ivib)
             if ( iut .ne. iufrom ) then
                ! add to the list of beads to be grown from iu
                temp_count = temp_count + 1
                temp_store(temp_count) = iut
             end if
          end do

          count = 0
          do izz = temp_count, 1, -1
             itry = int(real(izz,dp)*random(-1)) + 1
             iut = temp_store(itry)
             count = count + 1
             ! assign growlist for current index and count
             growlist(index,count) = iut
             lexshed(iut) = .false.

             ! update temp_store for next iteration
             temp_store(itry) = temp_store(izz)
          end do
          ! assign grownum for this index
          grownum(index) = count

          ! update list of "outer" beads
          ! remove bead that was just grown from outer list
          outer_sites(outer_try) = outer_sites(outer_num)
          outer_prev(outer_try) = outer_prev(outer_num)
          outer_num = outer_num - 1

          ! add the new beads if they have more to be grown
          iufrom = iu
          do icbu = 1,grownum(index)
             ! increment counter ibead for total number of beads grown
             ibead = ibead + 1

             ! determine whether this bead has any non-grown neighbors
             iu = growlist(index,icbu)
             invtry = invib(imolty,iu)

             if ( invtry .gt. 1 ) then
                ! add one to the stack of outer_sites
                outer_num = outer_num + 1
                outer_sites(outer_num) = iu
                outer_prev(outer_num) = iufrom
             end if
          end do
          ! end of while loop 70
          goto 70
       end if
    end if

100 if ( (movetype .eq. 1) .and. (ibead .gt. nmaxcbmc(imolty)) ) then
       kickout = kickout + 1
       goto 11
    end if

    if ( lprint ) then
       write(io_output,*) 'movetype',movetype
       write(io_output,*) 'index',index
       do ibead = 1,index
          write(io_output,*) 'ibead',ibead
          write(io_output,*) 'growfrom(ibead)',growfrom(ibead)
          write(io_output,*) 'growprev(ibead)',growprev(ibead)
          write(io_output,*) 'grownum(ibead)',grownum(ibead)
          do count = 1,grownum(ibead)
             write(io_output,*) 'count,growlist(ibead,count)',count ,growlist(ibead,count)
          end do
       end do
       do iu = 1,igrow
          write(io_output,*) 'iu,lexshed(iu)',iu,lexshed(iu)
       end do
    end if

    if (lrig(imolty).and.nrig(imolty).eq.0) then
       ib2 = index - ju
       if (ib2.lt.1) then
          llrig = .false.
          return
       else
          do ibead = 1, index
             lsave(ibead) = .false.
          end do

          nrigi = 1
          llrig = .true.
          rfrom(1) = growfrom(ib2)
          rprev(1) = growprev(ib2)
          rnum(1) = grownum(ib2)
          lsave(ib2) = .true.

          do count = 1, grownum(ib2)
             iu = growlist(ib2,count)
             lfind(iu) = .true.
             rlist(1,count) = iu
             lexshed(iu) = .false.
          end do
       end if

       ! cycle through the rest of the rigid beads
       do ibead = ib2 + 1, index
          iufrom = growfrom(ibead)
          if (lfind(iufrom)) then
             lsave(ibead) = .true.
             do count = 1, grownum(ibead)
                iu = growlist(ibead,count)
                lfind(iu) = .true.
             end do
          end if
       end do
    else
       llrig = .false.
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end SCHEDULE in ',myid
#endif
    return
  end subroutine schedule

!*********************************************************************
!> \brief Determines the new geometry of the bond lengths and angles to be
!> rotated on the cone for rosenb.f
!>
!> For old computes the rosenbluth weight for bending and determines
!> the old bond lenghts, angles, and \a phi for growth
!> \param bondlen(count) is the bondlengths from the grow bead to count
!> \param  bendang(count) is the bond angle between iuprev,iufrom, and count
!> \param  phi(count) is the angle around the cone between count and 1
!> \par History
!> written by M.G. Martin 7-10-98 from geomnew and geomold \n
!> last modified by Neeraj Rai on 12/23/2008 for CG models
!*********************************************************************
  subroutine geometry(lnew,iw,i,imolty,angstart,iuprev,glist,bondlen,bendang,phi,vvibtr,vbbtr,maxlen,wei_bend)
    use util_math,only:cone_angle
    use util_mp,only:mp_set_displs,mp_allgather

    ! variables passed to/from the subroutine
    logical::lnew
    integer::iw,imolty,angstart,iuprev,glist(numax),i
    real::bondlen(numax),bendang(numax),phi(numax),vvibtr,vbbtr,maxlen
    real::wei_bend

    ! local variables
    integer::count,ntogrow,iugrow,iufrom,iv,juvib,jtvib,iu2back,ib,iulast,type,aaa,iuone
    real::equil,kforce,vvib,length,angle,phitwo,vangle,vphi

    ! new variables
    integer::ibend,nchben_a,nchben_b
    real::bsum_try,rsint,ang_trial(nchbn_max),vbend(nchbn_max),bfactor(nchbn_max),rbf,bs

    ! variables from geomold
    real::rxui,ryui,rzui,rxuij,ryuij,rzuij,xvecprev,yvecprev,zvecprev,distprev,xvecgrow,yvecgrow,zvecgrow,distgrow,anglec
    real::xub,yub,zub,dum,ux,uy,uz

    ! Neeraj: Adding for the lookup table for CG model
    real::distgrow2
    real::lengtha,lengthb,lengtha2,lengthb2,lengthc,lengthc2,lengthFP,lengthFP2

    ! Q. Paul C. -- for tabulated CBMC bending growth
    logical::use_bend_table=.false. ! whether to use this type of growth
    real::delta_prob
    real,dimension(2)::theta1_interval,theta2_interval,phi12_interval ! the angle interval being picked
    integer::ileft,iright,imedian,iprob,ilin,ibr  ! used for searching the angle interval when given a random number (sorting)
    integer::itheta1,itheta2,iphi12
    integer::branch_num !number of branch points
    integer::type2
    real::phi12

    ! MPI
    integer::rcounts(numprocs),displs(numprocs),my_start,my_end,blocksize,my_itrial,rid
    real::my_ang_trial(nchbn_max),my_vbend(nchbn_max),my_bfactor(nchbn_max)

#ifdef __DEBUG__
    write(io_output,*) 'START GEOMETRY in ',myid
#endif

    ! assign grownum and growfrom to local variables
    ntogrow = grownum(iw)
    iufrom = growfrom(iw)

    ! initialize trial energies
    vvibtr = 0.0E0_dp

    if ( .not. lnew ) then
       ! OLD store r*ui positions for unit iufrom
       rxui = rxu(i,iufrom)
       ryui = ryu(i,iufrom)
       rzui = rzu(i,iufrom)
    end if

    ! Begin Bond length selection based on Boltzmann rejection

    ! determine the bond lengths of the beads to be grown
    maxlen = 0.0E0_dp

    do count = 1,ntogrow
       iugrow = growlist(iw,count)
       glist(count) = iugrow

       ! determine the vibration (bond) type as jtvib
       do iv = 1, invib(imolty,iugrow)
          juvib = ijvib(imolty,iugrow,iv)
          if ( juvib .eq. iufrom ) then
             jtvib = itvib(imolty,iugrow,iv)
             exit
          end if
       end do

       if ( lnew ) then
          ! compute bond length
          equil = brvib(jtvib)
          kforce = brvibk(jtvib)
          call bondlength( jtvib,equil,kforce,beta,length,vvib)
       else
          ! compute bond length
          rxuij = rxu(i,iugrow) - rxui
          ryuij = ryu(i,iugrow) - ryui
          rzuij = rzu(i,iugrow) - rzui

          ! do not need mimage for intramolecular
          length = sqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)

          ! compute vibration energy
          if (L_vib_table) then
             vvib = lininter_vib(length,jtvib)
          else
             equil = brvib(jtvib)
             kforce = brvibk(jtvib)
             vvib = kforce * (length-equil )**2
          end if
       end if

       ! adjust maximum bond length of those being grown
       if ( length .gt. maxlen ) maxlen = length

       ! assign bondlength and add up vibrational energy
       bondlen(count) = length
       vvibtr = vvibtr + vvib
    end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         --- Begin Bond angle biased selection ---             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if ( growprev(iw) .eq. 0 ) then
       ! need to choose one bead to grow on unit sphere
       iuprev = growlist(iw,1)
       angstart = 2
    else
       iuprev = growprev(iw)
       angstart = 1
    end if

    if ( .not. lnew ) then
       ! compute the vector from iufrom to iuprev
       xvecprev = rxu(i,iuprev) - rxui
       yvecprev = ryu(i,iuprev) - ryui
       zvecprev = rzu(i,iuprev) - rzui
       distprev = sqrt( xvecprev*xvecprev + yvecprev*yvecprev  + zvecprev*zvecprev )
    end if

    if (L_bend_table) then
       if ( growprev(iw) .eq. 0 ) then
          lengthFP = bondlen(1)
          lengthFP2 = lengthFP*lengthFP
       else
          if(.not.lnew) then
             lengthFP = distprev
             lengthFP2 = lengthFP*lengthFP
          else
             rxuij = rxnew(iuprev) - rxnew(iufrom)
             ryuij = rynew(iuprev) - rynew(iufrom)
             rzuij = rznew(iuprev) - rznew(iufrom)
             lengthFP2 = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
             lengthFP = sqrt(lengthFP2)
          end if
       end if
    end if

    ! initialize wei_bend
    wei_bend = 1.0E0_dp
    vbbtr = 0.0E0_dp

    ! Q. Paul C. -- Determine if CBMC_bend_table is used
    if ( L_cbmc_bend ) then
        ! linear case
        if ( ntogrow .eq. 1 .and. angstart .eq. 1) then
            branch_num = 0
            iugrow = growlist(iw,ntogrow)

            do ib = 1, inben(imolty,iugrow)
                iulast = ijben2(imolty,iugrow,ib)
                if ( iulast .eq. iufrom ) then
                    iu2back = ijben3(imolty,iugrow,ib)
                    if ( iu2back .eq. iuprev ) then
                        type = itben(imolty,iugrow,ib)
                        equil = brben(type)
                        kforce = brbenk(type)
                        exit
                    end if
                end if
            end do

            ! Look for the tabulated values corresponding to the angle type
            use_bend_table = .false.
            do ilin=1,size(lin_bend_type)
                if ( lin_bend_type(ilin) .eq. type ) then
                    use_bend_table = .true.
                    exit
                end if
            end do

        ! one-branch case
        else if ( ntogrow .eq. 2 .and. angstart .eq. 1 ) then
            branch_num = 1

            ! find force field parameter for theta1 and theta2
            do count = angstart,ntogrow
                iugrow = growlist(iw,count)
                do ib = 1, inben(imolty,iugrow)
                    iulast = ijben2(imolty,iugrow,ib)
                    if ( iulast .eq. iufrom ) then
                        iu2back = ijben3(imolty,iugrow,ib)
                        if ( iu2back .eq. iuprev ) then
                            if (count .eq. angstart) then
                                type = itben(imolty,iugrow,ib)
                                equil = brben(type)
                                kforce = brbenk(type)
                            else
                                type2 = itben(imolty,iugrow,ib)
                            end if
                            exit
                        end if
                    end if
                end do
            end do

            ! Look for the tabulated values corresponding to the angle type
            use_bend_table = .false.
            do ibr=1,size(br_bend_type)
                if ( (br_bend_type(ibr) .eq. type) .and. (type .eq. type2) ) then
                    use_bend_table = .true.
                    exit
                end if
            end do

        ! double-branch case, not implemented for tabulated CBMC growth yet
        else
            use_bend_table = .false.
        end if
    end if

    ! if use tabulated CBMC_bend_table to grow a linear bead
    if (use_bend_table .and. branch_num .eq. 0) then
        bsum_try = 0.0E0_dp

        if ( lnew ) then
            ! new conformation
            ! first, perform the biased selection for (almost) equally-likely space
            rid=-1
            rbf = random(rid)

            ! log(n) selection of the tabulated value
            ileft=1
            iright=lin_bend_dim(ilin)
            do while ( iright-ileft .gt. 1)
                imedian=(ileft+iright)/2
                if ( rbf .gt. lin_bend_prob(ilin,imedian)) then
                    ileft=imedian
                else
                    iright=imedian
                end if
            end do

            ! Make sure that it is the leftest one (in case some values are identical)
            if (ileft .gt. 1) then
                do while ( lin_bend_prob(ilin,ileft-1) .eq. lin_bend_prob(ilin,ileft) )
                    ileft=ileft-1
                    if (ileft .eq. 1) then
                        exit
                    end if
                end do
            end if

            theta1_interval(1)=lin_bend_table(ilin,ileft)
            theta1_interval(2)=lin_bend_table(ilin,iright)

            ! then perform random selection among the interval
            rbf = random(rid)
            angle=rbf*(theta1_interval(2)-theta1_interval(1))+theta1_interval(1)
            vangle = kforce * (angle - equil)**2
            delta_prob=(lin_bend_prob(ilin,iright)-lin_bend_prob(ilin,ileft))
            bsum_try = bsum_try + sin(angle)*exp(-beta*vangle)*(theta1_interval(2)-theta1_interval(1))/delta_prob
        else
            ! old conformation
            xvecgrow = rxu(i,iugrow) - rxui
            yvecgrow = ryu(i,iugrow) - ryui
            zvecgrow = rzu(i,iugrow) - rzui
            distgrow = bondlen(ntogrow)
            distgrow2 = distgrow*distgrow
            ! dot product divided by lengths gives cos(angle)
            anglec = ( xvecprev*xvecgrow + yvecprev*yvecgrow  + zvecprev*zvecgrow ) / (distprev*distgrow)
            angle = acos(anglec)
            vangle = kforce * (angle - equil)**2

            ! log(n) selection of the angle interval
            ileft=1
            iright=lin_bend_dim(ilin)
            do while ( iright-ileft .gt. 1)
                imedian=(ileft+iright)/2
                if ( angle .gt. lin_bend_table(ilin,imedian)) then
                    ileft=imedian
                else
                    iright=imedian
                end if
            end do

            theta1_interval(1)=lin_bend_table(ilin,ileft)
            theta1_interval(2)=lin_bend_table(ilin,iright)
            delta_prob=(lin_bend_prob(ilin,iright)-lin_bend_prob(ilin,ileft))
            bsum_try = bsum_try + sin(angle)*exp(-beta*vangle)*(theta1_interval(2)-theta1_interval(1))/delta_prob
        end if

        ! propagate the rosenbluth weight
        wei_bend = wei_bend * bsum_try
        bendang(angstart) = angle
        vbbtr = vbbtr + vangle

        ! The following part is the same as the normal growth
        if ( lnew ) then
            ! assign phi(angstart) to 0.0
            phi(angstart) = 0.0E0_dp
        else if (angstart.le.ntogrow) then
            ! set up the cone using iuprev
            xub = -xvecprev/distprev
            yub = -yvecprev/distprev
            zub = -zvecprev/distprev
            call cone(1,xub,yub,zub,dum,dum)

            iugrow = growlist(iw,angstart)

            ! compute vector from iufrom to iugrow
            xvecgrow = rxu(i,iugrow) - rxui
            yvecgrow = ryu(i,iugrow) - ryui
            zvecgrow = rzu(i,iugrow) - rzui
            distgrow = bondlen(angstart)

            ! turn this into a unit vector
            ux = xvecgrow/distgrow
            uy = yvecgrow/distgrow
            uz = zvecgrow/distgrow

            call cone(3,ux,uy,uz,bendang(angstart),phi(angstart))
        end if

    ! if use tabulated CBMC_bend_table to grow single-branch two beads
    else if (use_bend_table .and. branch_num .eq. 1) then
        bsum_try = 0.0E0_dp

        if ( lnew ) then
            ! new conformation
            ! first, perform the biased selection for (almost) equally-likely space
            rid=-1
            rbf = random(rid)

            ! log(n) selection of the tabulated value, and back calculate the corresponding theta1,theta2,phi12
            ileft=1
            iright=br_bend_dim1(ibr)*br_bend_dim2(ibr)*br_bend_dim3(ibr)+1
            do while ( iright-ileft .gt. 1)
                imedian=(ileft+iright)/2
                if ( rbf .gt. br_bend_prob(ibr,imedian)) then
                    ileft=imedian
                else
                    iright=imedian
                end if
            end do

            ! Make sure that it is the leftest one (in case some values are identical)
            if (ileft .gt. 1) then
                do while ( br_bend_prob(ibr,ileft-1) .eq. br_bend_prob(ibr,ileft))
                    ileft = ileft - 1
                    if (ileft .eq. 1) then
                        exit
                    end if
                end do
            end if

            ! Converting the probability into a particular grid
            iphi12=mod(ileft,br_bend_dim3(ibr))
            if ( iphi12 .eq. 0) then
                iphi12 = br_bend_dim3(ibr)
            end if

            if ( mod(ileft,br_bend_dim3(ibr)*br_bend_dim2(ibr)) .eq. 0) then
                itheta2 = br_bend_dim2(ibr)
            else if (mod(ileft,br_bend_dim3(ibr)) .eq. 0) then
                itheta2=mod(floor(ileft*1.0E0_dp/br_bend_dim3(ibr)),br_bend_dim2(ibr))
            else
                itheta2=mod(floor(ileft*1.0E0_dp/br_bend_dim3(ibr)),br_bend_dim2(ibr))+1
            end if

            if ( mod(ileft,br_bend_dim3(ibr)*br_bend_dim2(ibr)) .eq. 0) then
                itheta1=floor(ileft*1.0E0_dp/(br_bend_dim3(ibr)*br_bend_dim2(ibr)))
            else
                itheta1=floor(ileft*1.0E0_dp/(br_bend_dim3(ibr)*br_bend_dim2(ibr)))+1
            end if

            theta1_interval(1)=br_bend_theta1(ibr,itheta1)
            theta1_interval(2)=br_bend_theta1(ibr,itheta1+1)
            theta2_interval(1)=br_bend_theta2(ibr,itheta1,itheta2)
            theta2_interval(2)=br_bend_theta2(ibr,itheta1,itheta2+1)
            phi12_interval(1)=br_bend_phi12(ibr,itheta1,itheta2,iphi12)
            phi12_interval(2)=br_bend_phi12(ibr,itheta1,itheta2,iphi12+1)
            delta_prob=br_bend_prob(ibr,iright)-br_bend_prob(ibr,ileft)

            ! then perform random selection for theta1
            rbf = random(rid)
            angle=rbf*(theta1_interval(2)-theta1_interval(1))+theta1_interval(1)
            vangle = kforce * (angle - equil)**2
            bendang(1) = angle
            phi(1) = 0.0E0_dp
            vbbtr = vbbtr + vangle

            ! theta2
            rbf = random(rid)
            angle=rbf*(theta2_interval(2)-theta2_interval(1))+theta2_interval(1)
            vangle = kforce * (angle - equil)**2
            bendang(2) = angle
            vbbtr = vbbtr + vangle

            ! phi12
            rbf = random(rid)
            phitwo=rbf*(phi12_interval(2)-phi12_interval(1))+phi12_interval(1)

            ! calculate theta12 from phi12, and calculate the bending potential based on theta12
            angle=cone_angle(bendang(1),phi(1),bendang(2),phitwo)
            vangle = kforce * (angle - equil)**2
            vbbtr = vbbtr + vangle

            if ( phitwo .gt. onepi ) then
                phitwo=phitwo-twopi
            end if

            phi(2) = phitwo
            delta_prob=br_bend_prob(ibr,iright)-br_bend_prob(ibr,ileft)
            bsum_try = sin(bendang(1))*sin(bendang(2))*exp(-beta*vbbtr)* &
                (theta1_interval(2)-theta1_interval(1))*(theta2_interval(2)-theta2_interval(1))* &
                (phi12_interval(2)-phi12_interval(1))/delta_prob

        else
            ! old conformation
            ! first calculate theta1 and theta2, and find out the corresponding probability
            do count = 1, ntogrow
                iugrow = growlist(iw,count)
                xvecgrow = rxu(i,iugrow) - rxui
                yvecgrow = ryu(i,iugrow) - ryui
                zvecgrow = rzu(i,iugrow) - rzui
                distgrow = bondlen(count)
                distgrow2 = distgrow*distgrow

                ! dot product divided by lengths gives cos(angle)
                anglec = ( xvecprev*xvecgrow + yvecprev*yvecgrow  + zvecprev*zvecgrow ) / (distprev*distgrow)
                angle = acos(anglec)
                bendang(count) = angle
                vangle = kforce * (angle - equil)**2
                vbbtr = vbbtr + vangle

                ! calculate phi1 and phi2 corresponding to theta1 and theta2
                if ( count .eq. 1) then
                    xub = -xvecprev/distprev
                    yub = -yvecprev/distprev
                    zub = -zvecprev/distprev
                    call cone(1,xub,yub,zub,dum,dum)
                end if

                ux = xvecgrow/distgrow
                uy = yvecgrow/distgrow
                uz = zvecgrow/distgrow
                call cone(3,ux,uy,uz,bendang(count),phi(count))

                ! find the unit element corresponding to the angle
                if ( count .eq. 1) then
                    ileft=1
                    iright=br_bend_dim1(ibr)+1
                    do while ( iright-ileft .gt. 1)
                         imedian=(ileft+iright)/2
                         if ( angle .gt. br_bend_theta1(ibr,imedian)) then
                            ileft=imedian
                         else
                            iright=imedian
                         end if
                     end do

                     itheta1=ileft
                     theta1_interval(1)=br_bend_theta1(ibr,itheta1)
                     theta1_interval(2)=br_bend_theta1(ibr,itheta1+1)
                     iprob = (itheta1-1)*br_bend_dim2(ibr)*br_bend_dim3(ibr)
                else
                     ! find itheta2
                     ileft=1
                     iright=br_bend_dim2(ibr)+1
                     do while ( iright-ileft .gt. 1)
                        imedian=(ileft+iright)/2
                        if ( angle .gt. br_bend_theta2(ibr,itheta1,imedian)) then
                            ileft = imedian
                        else
                            iright = imedian
                        end if
                     end do

                     itheta2=ileft
                     theta2_interval(1)=br_bend_theta2(ibr,itheta1,itheta2)
                     theta2_interval(2)=br_bend_theta2(ibr,itheta1,itheta2+1)
                     iprob = iprob+(itheta2-1)*br_bend_dim3(ibr)

                     ! find iphi12
                     phi12=phi(1)-phi(2)
                     if ( phi12 .lt. 0) then
                        phi12=phi12+twopi
                     end if

                     ileft=1
                     iright=br_bend_dim3(ibr)+1
                     do while ( iright-ileft .gt. 1)
                        imedian=(ileft+iright)/2
                        if ( phi12 .gt. br_bend_phi12(ibr,itheta1,itheta2,imedian)) then
                            ileft=imedian
                        else
                            iright=imedian
                        end if
                     end do

                     iphi12=ileft
                     phi12_interval(1)=br_bend_phi12(ibr,itheta1,itheta2,iphi12)
                     phi12_interval(2)=br_bend_phi12(ibr,itheta1,itheta2,iphi12+1)
                     iprob = iprob+iphi12
                end if
            end do

            ! calculate theta12
            angle=cone_angle(bendang(1),phi(1),bendang(2),phi(2))
            vangle = kforce * (angle - equil)**2
            vbbtr = vbbtr + vangle

            ! weight calculation
            delta_prob=(br_bend_prob(ibr,iprob+1)-br_bend_prob(ibr,iprob))
            bsum_try = sin(bendang(1))*sin(bendang(2))*exp(-beta*vbbtr)* &
              (theta1_interval(2)-theta1_interval(1))*(theta2_interval(2)-theta2_interval(1))* &
              (phi12_interval(2)-phi12_interval(1))/delta_prob
        end if

        ! propagate the rosenbluth weight
        wei_bend = wei_bend * bsum_try

    ! if not use tabulated CBMC_bend_table, regular CBMC growth
    else

        nchben_a = nchbna(imolty)
        ! MPI
        if (numprocs.gt.1) then
            rid=myid
        else
            rid=-1
        end if
        blocksize = nchben_a/numprocs
        rcounts = blocksize
        blocksize = nchben_a - blocksize * numprocs
        if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
            call mp_set_displs(rcounts,displs,blocksize,numprocs)
            my_start = displs(myid+1) + 1
            my_end = my_start + rcounts(myid+1) - 1

            ! determine the iugrow-iufrom-iuprev angles
            do count = angstart,ntogrow
                iugrow = growlist(iw,count)
                kforce = -1000.0E0_dp
                do ib = 1, inben(imolty,iugrow)
                    iulast = ijben2(imolty,iugrow,ib)
                    if ( iulast .eq. iufrom ) then
                        iu2back = ijben3(imolty,iugrow,ib)
                        if ( iu2back .eq. iuprev ) then
                            type = itben(imolty,iugrow,ib)
                            equil = brben(type)
                            kforce = brbenk(type)
                            exit
                        end if
                    end if
                end do

            if ( kforce .gt. 0.1E0_dp ) then
                ! flexible bond angle
                ! initialize bsum_try
                bsum_try = 0.0E0_dp

                if (L_bend_table) then
                    distgrow = bondlen(count)
                    distgrow2 = distgrow*distgrow
                end if

                ! compute trial angles and energies
                my_itrial = 0
                do ibend = my_start,my_end
                    my_itrial = my_itrial + 1
                    if (lnew.or.ibend.ne.1) then
                        ! new conformation or old conformation skipping 1st
                        ! choose the angle uniformly on sin(theta)
                        rsint = 2.0E0_dp*random(rid) - 1.0E0_dp
                        angle = acos(rsint)

                        ! calculate the bond angle energy
                        if (L_bend_table) then
                            lengthc2 = lengthFP2 + distgrow2 - 2.0E0_dp*lengthFP*distgrow*cos(angle)
                            lengthc = sqrt(lengthc2)
                            vangle = lininter_bend(lengthc,type)
                        else
                            vangle = kforce * (angle - equil)**2
                        end if

                        my_ang_trial(my_itrial) = angle
                        my_vbend(my_itrial)=vangle
                        my_bfactor(my_itrial) = exp(-beta*vangle)
                        bsum_try = bsum_try + my_bfactor(my_itrial)
                    else
                        ! first ibend is the old conformation
                        xvecgrow = rxu(i,iugrow) - rxui
                        yvecgrow = ryu(i,iugrow) - ryui
                        zvecgrow = rzu(i,iugrow) - rzui
                        distgrow = bondlen(count)
                        distgrow2 = distgrow*distgrow
                        ! dot product divided by lengths gives cos(angle)
                        anglec = ( xvecprev*xvecgrow + yvecprev*yvecgrow  + zvecprev*zvecgrow ) / (distprev*distgrow)
                        if(anglec.gt.1.0E0_dp) then
                           angle = 0.0E0_dp
                        else if(anglec.lt.-1.0E0_dp) then
                           angle = onepi
                        else
                           angle = acos(anglec)
                        end if
                        if (L_bend_table) then
                            lengthc2 = lengthFP2 + distgrow2 - 2.0E0_dp*lengthFP*distgrow*anglec
                            lengthc = sqrt(lengthc2)
                            vangle = lininter_bend(lengthc,type)
                        else
                            vangle = kforce * (angle - equil)**2
                        end if

                        my_ang_trial(1) = angle
                        my_vbend(1)=vangle
                        my_bfactor(1) = exp( -beta*vangle )
                        bsum_try = bsum_try + my_bfactor(1)
                    end if
                end do

                if ( lnew ) then
                    ! select one of the trial sites via bias
                    rbf = random(rid)*bsum_try
                    bs = 0.0E0_dp
                    do ibend = 1,rcounts(myid+1)
                        bs = bs + my_bfactor(ibend)
                        if ( rbf .lt. bs ) then
                            angle = my_ang_trial(ibend)
                            vangle=my_vbend(ibend)
                            exit
                        end if
                    end do
                else
                    ! select the old conformation
                    angle = my_ang_trial(1)
                    vangle=my_vbend(1)
                end if

                if (numprocs.gt.1) then
                    call mp_allgather(bsum_try,bfactor,groupid)
                    bsum_try=sum(bfactor(1:numprocs))
                    call mp_allgather(angle,ang_trial,groupid)
                    call mp_allgather(vangle,vbend,groupid)

                    if ( lnew ) then
                        ! select one of the trial sites via bias
                        rbf = random(-1)*bsum_try
                        bs = 0.0E0_dp
                        do ibend = 1,numprocs
                            bs = bs + bfactor(ibend)
                            if ( rbf .lt. bs ) then
                                angle = ang_trial(ibend)
                                vangle=vbend(ibend)
                                exit
                            end if
                        end do
                    else
                        ! select the old conformation
                        angle = ang_trial(1)
                        vangle=vbend(1)
                    end if
                end if

                ! propagate the rosenbluth weight
                wei_bend = wei_bend * bsum_try/dble(nchben_a)

            else if (kforce.lt.-0.1E0_dp) then
                ! freely-joint beads
                rsint = 2.0E0_dp*random(-1) - 1.0E0_dp
                angle = acos(rsint)
                vangle=0.0E0_dp
            else
                ! fixed bond angle
                angle = equil
                vangle = 0.0E0_dp
            end if
            bendang(count) = angle
            vbbtr = vbbtr + vangle
        end do

        ! Neeraj  iugrow-iufrom-iuprev bend angle has been selected!
        if ( lnew ) then
            ! assign phi(angstart) to 0.0
            phi(angstart) = 0.0E0_dp
        else if (angstart.le.ntogrow) then
            ! set up the cone using iuprev
            xub = -xvecprev/distprev
            yub = -yvecprev/distprev
            zub = -zvecprev/distprev
            call cone(1,xub,yub,zub,dum,dum)

            iugrow = growlist(iw,angstart)

            ! compute vector from iufrom to iugrow
            xvecgrow = rxu(i,iugrow) - rxui
            yvecgrow = ryu(i,iugrow) - ryui
            zvecgrow = rzu(i,iugrow) - rzui
            distgrow = bondlen(angstart)

            ! turn this into a unit vector
            ux = xvecgrow/distgrow
            uy = yvecgrow/distgrow
            uz = zvecgrow/distgrow

            call cone(3,ux,uy,uz,bendang(angstart),phi(angstart))
        end if

        nchben_b = nchbnb(imolty)
        ! MPI
        blocksize = nchben_b/numprocs
        rcounts = blocksize
        blocksize = nchben_b - blocksize * numprocs
        if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
        call mp_set_displs(rcounts,displs,blocksize,numprocs)
        my_start = displs(myid+1) + 1
        my_end = my_start + rcounts(myid+1) - 1

        ! determine the angles of the grown beads relative to anglestart
        ! skip angstart in the loop below
        do count = angstart+1,ntogrow
            iugrow = growlist(iw,count)

            ! initialize bsum_try
            bsum_try = 0.0E0_dp

            if (L_bend_table) then
                lengtha = bondlen(count)
                lengtha2 = lengtha * lengtha
            end if

            ! compute trial energies and weights
            my_itrial = 0
            do ibend = my_start,my_end
                my_itrial = my_itrial + 1
                vphi = 0.0E0_dp
                if (lnew.or.ibend.ne.1) then
                    ! perform all ibend for NEW and OLD (except for 1st in OLD)
                    ! determine a value of phitwo
                    phitwo = random(rid)*twopi

                    do aaa = angstart,count-1
                        iuone = growlist(iw,aaa)
                        angle=cone_angle(bendang(aaa),phi(aaa),bendang(count),phitwo)
                        do ib = 1, inben(imolty,iugrow)
                            iulast = ijben2(imolty,iugrow,ib)
                            if ( iulast .eq. iufrom ) then
                                iu2back = ijben3(imolty,iugrow,ib)
                                if ( iu2back .eq. iuone ) then
                                    type = itben(imolty,iugrow,ib)
                                    ! calculate the bond angle energy
                                    if (L_bend_table) then
                                        lengthb = bondlen(aaa)
                                        lengthb2 = lengthb*lengthb
                                        lengthc2 = lengtha2 + lengthb2 - 2.0E0_dp*lengtha*lengthb*cos(angle)
                                        lengthc = sqrt(lengthc2)
                                        vphi = vphi + lininter_bend (lengthc,type)
                                    else
                                        vphi = vphi + brbenk(type) * (angle - brben(type))**2
                                    end if
                                end if
                            end if
                        end do
                    end do

                    ! store the boltzmann factors and phi
                    my_ang_trial(my_itrial) = phitwo
                    my_vbend(my_itrial)=vphi
                    my_bfactor(my_itrial) = exp(-beta*vphi)
                    bsum_try = bsum_try + my_bfactor(my_itrial)
                else
                    ! compute vector from iufrom to iugrow
                    xvecgrow = rxu(i,iugrow) - rxui
                    yvecgrow = ryu(i,iugrow) - ryui
                    zvecgrow = rzu(i,iugrow) - rzui
                    distgrow = bondlen(count)

                    ! turn this into a unit vector
                    ux = xvecgrow/distgrow
                    uy = yvecgrow/distgrow
                    uz = zvecgrow/distgrow

                    call cone(3,ux,uy,uz,bendang(count),phitwo)

                    do aaa = angstart,count-1
                        iuone = growlist(iw,aaa)
                        angle=cone_angle(bendang(aaa),phi(aaa),bendang(count),phitwo)

                        do ib = 1, inben(imolty,iugrow)
                            iulast = ijben2(imolty,iugrow,ib)
                            if ( iulast .eq. iufrom ) then
                                iu2back = ijben3(imolty,iugrow,ib)
                                if ( iu2back .eq. iuone ) then
                                    type = itben(imolty,iugrow,ib)
                                    ! calculate the bond angle energy
                                    if (L_bend_table) then
                                        lengthb = bondlen(aaa)
                                        lengthb2 = lengthb * lengthb
                                        lengthc2 = lengtha2 + lengthb2 - 2.0E0_dp*lengtha*lengthb*cos(angle)
                                        lengthc = sqrt(lengthc2)
                                        vphi = vphi + lininter_bend(lengthc,type)
                                    else
                                        vphi = vphi + brbenk(type) * (angle - brben(type))**2
                                    end if
                                end if
                            end if
                        end do
                    end do

                    my_ang_trial(1) = phitwo
                    my_vbend(1)=vphi
                    my_bfactor(1) = exp( -beta * vphi )
                    bsum_try = bsum_try + my_bfactor(1)
                end if
            end do

            if ( lnew ) then
                ! select a value of phitwo in a biased fashion
                rbf = random(rid)*bsum_try
                bs = 0.0E0_dp
                do ibend = 1,rcounts(myid+1)
                    bs = bs + my_bfactor(ibend)
                    if ( rbf .lt. bs ) then
                        phitwo = my_ang_trial(ibend)
                        vphi=my_vbend(ibend)
                        exit
                    end if
                end do
            else
                ! select the OLD value of phitwo
                phitwo = my_ang_trial(1)
                vphi=my_vbend(1)
            end if

            if (numprocs.gt.1) then
                call mp_allgather(bsum_try,bfactor,groupid)
                bsum_try=sum(bfactor(1:numprocs))
                call mp_allgather(phitwo,ang_trial,groupid)
                call mp_allgather(vphi,vbend,groupid)

                if ( lnew ) then
                    ! select a value of phitwo in a biased fashion
                    rbf = random(-1)*bsum_try
                    bs = 0.0E0_dp
                    do ibend = 1,numprocs
                        bs = bs + bfactor(ibend)
                        if ( rbf .lt. bs ) then
                            phitwo = ang_trial(ibend)
                            vphi=vbend(ibend)
                            exit
                        end if
                    end do
                else
                    ! select the OLD value of phitwo
                    phitwo = ang_trial(1)
                    vphi=vbend(1)
                end if
            end if

            ! propagate angle weight
            wei_bend = wei_bend * (bsum_try/dble(nchben_b))

            ! store the angle for phitwo
            phi(count) = phitwo
            vbbtr = vbbtr + vphi
        end do
    end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End Bond angle biased selection ---           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef __DEBUG__
    write(io_output,*) 'FINISH GEOMETRY in ',myid
#endif
    return
  end subroutine geometry

!**********************************************************
!> \brief choose a bond angle
!>
!> \param equil is equilibrium bond angle
!> \param kforce is the force constant
!> \param betaT is 1/kT
!> \param angle is the angle returned by this subroutine
!> \param vangle is the energy of the angle
!> \author M.G. Martin
!**********************************************************
  subroutine bendangle(equil,kforce,betaT,angle,vangle)
    real::equil,kforce,betaT,angle,vangle,rr,v1,v2

#ifdef __DEBUG__
    write(io_output,*) 'start BENDANGLE in ',myid
#endif

    if ( kforce .gt. 0.1E0_dp ) then
       ! find a vector inside the unit sphere
80     v1 = 2.0E0_dp*random(-1)-1.0E0_dp
       v2 = 2.0E0_dp*random(-1)-1.0E0_dp
       rr = v1*v1 + v2*v2
       if (rr .ge. 1.0E0_dp ) goto 80

       ! select angle from a gaussian distribution
       angle = equil + v1*sqrt( (-log(rr)) /( kforce*betaT*rr) )

       if (angle .le. 0.0E0_dp .or. angle .ge. 3.1415 ) then
          write(io_output,*) 'chose angle outside of 0,Pi in bendangle'
          goto 80
       end if

       ! correct for the phase space of the angle
       if ( random(-1) .gt. sin(angle) ) goto 80
       ! if (L_bend_table) then
       ! call lininter_bend(angle, tabulated_bend, type)
       ! type isn't specified in this subroutine, but
       ! I don't think its called anywhere anymore
       vangle = kforce * (angle - equil)**2
    else
       ! fixed bond angle
       angle = equil
       vangle = 0.0E0_dp
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end BENDANGLE in ',myid
#endif
    return
  end subroutine bendangle

!*****************************************************************
!> note that the rotation matrix is saved after each call so it
!> needs to be reset when you wish to use another cone
!>
!> \param iinit if iinit = 1 then it sets up the rotation matrix for the cone
!> using x,y,z as a unit vector pointing in the +z direction \n
!> if iinit = 2 then it creates a unit vector that has an angle
!> of alpha from the +z direction (previous vector) and an angle
!> of gamma (0,2Pi) around the cone circle and returns this as
!> x,y,z \n
!> if iinit = 3 then this computes gamma (0,2Pi) for a unit
!> vector x,y,z with angle alpha from the -z direction
!> \par History
!> originally written prior to 1995 \n
!> last modified 02-12-2001 by M.G. Martin
!*****************************************************************
  subroutine cone(iinit,x,y,z,alpha,gamma)
    ! variables passed to/from the subroutine
    integer,intent(in)::iinit
    real,intent(in)::alpha
    real::x,y,z,gamma

    ! local variables
    real,save::a11,a12,a13,a21,a22,a31,a32,a33,cosalph,sinalph
    real::determ,sinthe,costhe,sinpsi,cospsi,singamma,cosgamma,uxtemp,uytemp,uztemp

#ifdef __DEBUG_VERBOSE__
    write(io_output,*) 'start CONE in ',myid
#endif

    if ( iinit .eq. 1 ) then
       ! setup the unit cone
       ! assume x,y,z is a rotation of vector (0,0,1)
       ! conversions from xyz to theta,psi
       ! x = -sin(theta)*cos(psi)
       ! y = sin(theta)*sin(psi)
       ! z = cos(theta)
       ! rotation matrix (note that sin(phi) = 0, cos(phi) = 1)
       ! a11 = cos(psi) cos(theta) cos(phi) - sin(psi) sin(phi)
       ! cos(psi) cos(theta)
       ! a12 = -sin(psi) cos(theta) cos(phi) - cos(psi) sin(phi)
       ! sin(psi) cos(theta)
       ! a13 = sin(theta) cos(phi)
       ! sin(theta)
       ! a21 = cos(psi) cos(theta) sin(phi) + sin (psi) cos(phi)
       ! sin(psi)
       ! a22 = -sin(psi) cos(theta) sin(phi) + cos (psi) cos(phi)
       ! cos(psi)
       ! a23 = sin(theta) sin(phi)
       ! 0.0
       ! a31 = -cos(psi) sin(theta)
       ! x
       ! a32 = sin(psi)  sin (theta)
       ! y
       ! a33 = cos(theta)
       ! z

       costhe = z
       sinthe = sqrt(1.0E0_dp - costhe*costhe)
       if (sinthe .ne. 0.0E0_dp) then
          sinpsi=(y)/sinthe
          cospsi=(-x)/sinthe
       else
          sinpsi = 0.0E0_dp
          cospsi = 1.0E0_dp
       end if

       ! rotation matrix
       a11 = cospsi*costhe
       a12 = -(sinpsi*costhe)
       a13 = sinthe
       a21 = sinpsi
       a22 = cospsi
       a31 = x
       a32 = y
       a33 = z
    else if ( iinit .eq. 2 ) then
       ! find the vector on a cone:
       ! alpha is the angle with the cone axis (theta)
       sinalph = sin(alpha)
       cosalph = cos(alpha)
       ! gamma is now passed to cone - (phi)
       uztemp = -cosalph
       uytemp = sinalph*cos(gamma)
       uxtemp = sinalph*sin(gamma)
       x = a11*uxtemp + a21*uytemp + a31*uztemp
       y = a12*uxtemp + a22*uytemp + a32*uztemp
       z = a13*uxtemp +              a33*uztemp
    else if ( iinit .eq. 3 ) then
       ! find gamma for a given unit vector on the cone
       ! use cramer's rule where a23 has already been replaced with 0
       ! we don't need the uztemp so it is not computed
       determ =  ( a11*a22*a33 + a21*(a32*a13 - a12*a33) - a31*a22*a13 )
       uxtemp =   ( x*a22*a33 + a21*(a32*z - y*a33) - a31*a22*z ) / determ
       uytemp =  ( a11*(y*a33 - a32*z) + x*(a32*a13 - a12*a33) + a31*(a12*z - y*a13) ) / determ

       sinalph = sin(alpha)
       singamma = uxtemp/sinalph
       cosgamma = uytemp/sinalph

       ! now need to find the gamma on [-Pi,Pi] that satisfies cos and sin
       if ( cosgamma .gt. 1.0E0_dp ) then
          gamma = 0.0E0_dp
       else if ( cosgamma .lt. -1.0E0_dp ) then
          gamma = onepi
       else
          gamma = acos(cosgamma)
       end if
       if ( singamma .lt. 0.0E0_dp ) gamma = -gamma
    else
       write(io_output,*) 'iinit ',iinit
       call err_exit(__FILE__,__LINE__,'non valid iinit in cone.f',myid+1)
    end if

#ifdef __DEBUG_VERBOSE__
    write(io_output,*) 'finish CONE in ',myid
#endif

    return
  end subroutine cone

!***********************************************************
!> \brief Computes the bond length for a given vibration type
!>
!> \param vibtype is the vibration type
!> \param requil is the equilibrium bond length
!> \param kvib is the force constant for the bond length
!> \param length is the bond length returned by this subroutine
!> \param betaT is 1/kT
!> \param vvib is the vibration energy for this bond
!> \author M.G. Martin  2-4-98
!***********************************************************
  subroutine bondlength(vibtype,requil,kvib,betaT,length,vvib)
    integer::vibtype
    real::length,bond,bf,vvib,betaT,kvib,requil,minRegrow,maxRegrow,regrowWidth,maxRegrowSq
    maxRegrow = maxRegrowVib(vibtype)
    maxRegrowSq = maxRegrow**2
    minRegrow = minRegrowVib(vibtype)
    vvib = 0.0E0_dp
    regrowWidth = maxRegrow-minRegrow
    if ( kvib .gt. 0.1E0_dp ) then
       ! random bond length from Boltzmann distribution ---
       ! selects a bondlength from the minimum to the maximum specified in
       ! fort.4
       ! for efficiency sake we first select a value from the distribution
       ! P(L) =  L^2 dL
107    bond = (regrowWidth*random(-1)+minRegrow)

       ! correct for jacobian by dividing by the max^2
       bf = bond*bond/(maxRegrowSq)
       if ( random(-1) .ge. bf ) goto 107
       bond = bond * requil

       ! correct for the bond energy by applying the boltzmann
       ! weight to the bond length selection.
       vvib = kvib * (bond-requil )**2
       bf = exp ( -(vvib * betaT) )
       if ( random(-1) .ge. bf ) goto 107
       length = bond

       ! tabulated bond stretching potential
       ! added 12/02/08 by KM
    else if (L_vib_table) then
       ! random bond length from Boltzmann distribution ---
       !     ---  +/- 25% of equilibrium bond length
108    bond = (regrowWidth*random(-1)+minRegrow)

    ! correct for jacobian by dividing by the max^2
       bf = bond*bond*bond/(maxRegrowSq)
       if ( random(-1) .ge. bf ) goto 108

       bond = bond * requil
       vvib=lininter_vib(bond,vibtype)

       bf = exp ( -(vvib * betaT) )

       if ( random(-1) .ge. bf ) goto 108
       length = bond

    else
       ! fixed bond length
       length = requil
    end if

    return
  end subroutine bondlength

  !> \brief Translate and rotate the rigid body such that the first and second beads match with given ones
  !> 1) Translate so that bead 1 matches with given bead 1
  !> 2) Rotate so that bead 2 matches with given bead 2
  !> If bead 3 is given as well, call subroutine align_rotate to do the rigid rotation around vector 12
  !> \param first_bead, second_bead: the two bead numbers that need to be aligned with two target beads
  !> \param nbead: number of repeat unit beads to move
  !> \param xi, yi, zi: input and output coordinates
  !> \param xtarget, ytarget, ztarget: size of 2, the coordinates of 2 given beads
  subroutine align_lines(first_bead, second_bead, nbead, xi, yi, zi, xtarget, ytarget, ztarget)
      use util_runtime, only : err_exit
      use util_math, only : cross_product,normalize_vector
      integer, intent(in) :: first_bead, second_bead, nbead
      real, intent(inout) :: xi(numax), yi(numax), zi(numax)
      real, intent(in) :: xtarget(2), ytarget(2), ztarget(2)
      real :: Qrot(3,3)
      real :: xtrans, ytrans, ztrans, costheta, sintheta, smallc, bigc, s, bead_err
      real :: vector_a(3), vector_b(3), vector_n(3), vector_w(3), vector_rot(3)
      real :: vector_diff(3), vector_cross(3)
      integer :: ibead

      !----- Translation part
      xtrans = xtarget(1) - xi(first_bead)
      ytrans = ytarget(1) - yi(first_bead)
      ztrans = ztarget(1) - zi(first_bead)

      do ibead = 1, nbead
          xi(ibead) = xi(ibead) + xtrans
          yi(ibead) = yi(ibead) + ytrans
          zi(ibead) = zi(ibead) + ztrans
      end do

      !----- Rotation part
      ! use cross_product to determine the normal vector to the plane
      ! the plane is determined by two vectors, xi(1)-->xi(2) and xtarget(1)-->xtarget(2)
      vector_a(1) = xtarget(2) - xtarget(1)
      vector_a(2) = ytarget(2) - ytarget(1)
      vector_a(3) = ztarget(2) - ztarget(1)
      vector_b(1) = xi(second_bead) - xi(first_bead)
      vector_b(2) = yi(second_bead) - yi(first_bead)
      vector_b(3) = zi(second_bead) - zi(first_bead)

      ! if (in the old growth) vector_a and vector_b are the same, skip the rotation part to avoid numerical issues
      vector_diff = vector_a - vector_b
      if (dot_product(vector_diff,vector_diff) .le. 1e-8) then
          return
      end if

      vector_cross = cross_product(vector_a, vector_b)
      vector_n = normalize_vector(vector_cross)
      vector_w = cross_product(vector_n, vector_a) !< vector_w and vector_a are the orthogonal vector in the planeof vector_a and vector_b

      ! calculate angle theta to rotate
      costheta = dot_product(vector_a, vector_b) / (sqrt(dot_product(vector_a,vector_a)*dot_product(vector_b,vector_b)))
      sintheta = sqrt(dot_product(vector_cross, vector_cross) / (dot_product(vector_a,vector_a)*dot_product(vector_b,vector_b)))
      s = atan2(dot_product(vector_w, vector_b), dot_product(vector_a, vector_b))
      if (atan2(dot_product(vector_w, vector_b), dot_product(vector_a, vector_b)) .gt. 0) then
          sintheta = -sintheta
      end if

      ! set up the rotation matrix (see wikipedia "rotation matrix, axis and angle" for theoretical details)
      ! for clarity, the same notation is used as in wikipedia (smallc, bigc, s)
      smallc = costheta
      s = sintheta
      bigc = 1 - smallc
      Qrot(1, 1) = vector_n(1) * vector_n(1) * bigc + smallc
      Qrot(1, 2) = vector_n(2) * vector_n(1) * bigc - vector_n(3) * s
      Qrot(1, 3) = vector_n(1) * vector_n(3) * bigc + vector_n(2) * s
      Qrot(2, 1) = vector_n(2) * vector_n(1) * bigc + vector_n(3) * s
      Qrot(2, 2) = vector_n(2) * vector_n(2) * bigc + smallc
      Qrot(2, 3) = vector_n(2) * vector_n(3) * bigc - vector_n(1) * s
      Qrot(3, 1) = vector_n(3) * vector_n(1) * bigc - vector_n(2) * s
      Qrot(3, 2) = vector_n(3) * vector_n(2) * bigc + vector_n(1) * s
      Qrot(3, 3) = vector_n(3) * vector_n(3) * bigc + smallc

      ! do the rotation, and determine the new coordinate for xi, yi and zi
      do ibead = 1, nbead
          vector_b(1) = xi(ibead) - xi(first_bead)
          vector_b(2) = yi(ibead) - yi(first_bead)
          vector_b(3) = zi(ibead) - zi(first_bead)
          vector_rot = matmul(Qrot, vector_b) ! the new vector after rotation
          xi(ibead) = xi(first_bead) + vector_rot(1)
          yi(ibead) = yi(first_bead) + vector_rot(2)
          zi(ibead) = zi(first_bead) + vector_rot(3)
      end do

      ! make sure it matches
      bead_err = abs((xi(second_bead) - xtarget(2))**2 + (yi(second_bead) - ytarget(2))**2 + (zi(second_bead) - ztarget(2))**2)
      if (bead_err .gt. 1e-4) then
         call err_exit(__FILE__,__LINE__,'error in align_lines subroutine, the second bead does not match with the target',-1)
      end if

  end subroutine align_lines

  !> \brief fix the first (or last) two beads, and rotate the third one by angle phi (0, 2pi)
  !> \param bead1_coord, bead2_coord: coordinates of first (or last) two beads that serve as axis, dim = 3
  !> \param nbead: number of repeat unit beads
  !> \param xi, yi, zi: input coordinates
  !> \param phi: the input angle phi (0, 2pi)
  !> \param xout, yout, zout: output coordinates
  subroutine dihedral_rigrot(bead1_coord, bead2_coord, nbead, xi, yi, zi, phi, xout, yout, zout)
      use util_runtime, only : err_exit
      use util_math, only : normalize_vector
      integer, intent(in) :: nbead
      real, intent(in) :: xi(numax), yi(numax), zi(numax), bead1_coord(3), bead2_coord(3)
      real, intent(out) :: xout(numax), yout(numax), zout(numax)
      real, intent(inout) :: phi

      integer :: ibead
      real :: vector_prev(3), vector_grow(3), vector_prev_normalized(3), W(3, 3), Qrot(3, 3)
      real :: I(3, 3)

      xout = xi
      yout = yi
      zout = zi

      ! compute the unit vector as the rotation axis
      vector_prev = bead2_coord - bead1_coord
      vector_prev_normalized = normalize_vector(vector_prev)

      ! set up the rotation matrix (according to Rodrigues rotation formula)
      W(1, 1) = 0.0
      W(1, 2) = - vector_prev_normalized(3)
      W(1, 3) = vector_prev_normalized(2)
      W(2, 1) = vector_prev_normalized(3)
      W(2, 2) = 0.0
      W(2, 3) = - vector_prev_normalized(1)
      W(3, 1) = - vector_prev_normalized(2)
      W(3, 2) = vector_prev_normalized(1)
      W(3, 3) = 0.0
      I = 0.0
      I(1, 1) = 1.0
      I(2, 2) = 1.0
      I(3, 3) = 1.0

      Qrot = I + sin(phi) * W + (2. * sin(.5 * phi)**2) * matmul(W, W)

      ! loop over all unit and grow
      do ibead = 1, nbead
          vector_grow(1) = xout(ibead) - bead2_coord(1)
          vector_grow(2) = yout(ibead) - bead2_coord(2)
          vector_grow(3) = zout(ibead) - bead2_coord(3)

          vector_grow = matmul(Qrot, vector_grow)

          xout(ibead) = vector_grow(1) + bead2_coord(1)
          yout(ibead) = vector_grow(2) + bead2_coord(2)
          zout(ibead) = vector_grow(3) + bead2_coord(3)
      end do

  end subroutine dihedral_rigrot


!************************************************************
!> \brief Performs a rotational configurational bias move
!>
!> \attention all rigid beads should come after riutry
!************************************************************
  subroutine rigrot(lnew,lterm,iskip,imol,imolty,ibox,wadd,first_bead,second_bead,mvtp)
    logical::lnew,ovrlap,lterm,ltors

    integer::ibox,igrow,j,ip,iwalk,iunit,imolty,iu,iskip,igrow2
    integer::ichoi,imol,glist(numax),ntogrow,count,mvtp,first_bead,second_bead

    real::rx,ry,rz,rxorig,ryorig,rzorig,rxnw,rynw,rznw,phi,ax1(3),ax2(3)
    real::xdgamma,ydgamma,zdgamma,xcosdg,xsindg,ycosdg,ysindg,zcosdg,zsindg,rbf,rxur(numax),ryur(numax),rzur(numax),length
    real::w,bsum,wadd,maxlen,bs
! ----------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start RIGROT in ',myid
#endif

! initialize conformation energies and weight
    wadd = 1.0E0_dp
    w = 0.0E0_dp
    ichoi = nchoir(imolty)
    ltors = .false.
    iunit = nunit(imolty)
    if(mvtp.eq.0) then  !BX
       igrow = first_bead
       igrow2 = 0
    else if(mvtp.eq.5) then
       igrow = second_bead
       igrow2 = first_bead
    else
       igrow = riutry(imolty,1)
       igrow2 = 0
    end if

    if (lnew) then
       rxorig = rxnew(igrow)
       ryorig = rynew(igrow)
       rzorig = rznew(igrow)
       do j = igrow+1, iunit
          rxur(j) = rxnew(j)
          ryur(j) = rynew(j)
          rzur(j) = rznew(j)
       end do
    else
       rxorig = rxu(imol,igrow)
       ryorig = ryu(imol,igrow)
       rzorig = rzu(imol,igrow)
       do j = igrow+1, iunit
          rxur(j) = rxu(imol,j)
          ryur(j) = ryu(imol,j)
          rzur(j) = rzu(imol,j)
       end do
    end if

    ! find maxlen
    maxlen = 0.0E0_dp
    do j = igrow+1, iunit
       length = (rxur(j)-rxorig)**2+(ryur(j)-ryorig)**2 +(rzur(j)-rzorig)**2
       if (length.gt.maxlen) maxlen = length
    end do
    maxlen = sqrt(maxlen)

    if(mvtp.eq.0.or.mvtp.eq.2) then
       do ip = 1, ichoi
          if (lnew.or.ip.ne.1) then
             xdgamma = twopi*random(-1)
             ydgamma = twopi*random(-1)
             zdgamma = twopi*random(-1)
             ! set up rotation matrix
             xcosdg = cos(xdgamma)
             xsindg = sin(xdgamma)
             ycosdg = cos(ydgamma)
             ysindg = sin(ydgamma)
             zcosdg = cos(zdgamma)
             zsindg = sin(zdgamma)

             ! set molecule to rotate around
             ! rotate around all axis
             count = 1
             do j = igrow+1, iunit
                ry = ryur(j) - ryorig
                rz = rzur(j) - rzorig
                rynw = xcosdg * ry + xsindg * rz
                rznw = xcosdg * rz - xsindg * ry
                ryp(count,ip) = ryorig + rynw
                rzp(count,ip) = rzorig + rznw

                rx = rxur(j) - rxorig
                rz = rzp(count,ip) - rzorig
                rxnw = ycosdg * rx - ysindg * rz
                rznw = ycosdg * rz + ysindg * rx
                rxp(count,ip) = rxorig + rxnw
                rzp(count,ip) = rzorig + rznw

                rx = rxp(count,ip) - rxorig
                ry = ryp(count,ip) - ryorig
                rxnw = zcosdg * rx + zsindg * ry
                rynw = zcosdg * ry - zsindg * rx
                rxp(count,ip) = rxorig + rxnw
                ryp(count,ip) = ryorig + rynw

                count = count + 1
             end do

          else
             count = 1
             do j = igrow+1, iunit
                rxp(count,ip) = rxu(imol,j)
                ryp(count,ip) = ryu(imol,j)
                rzp(count,ip) = rzu(imol,j)
                count = count + 1
             end do
          end if
       end do

    !BX: for nsampos.eq.2
    else

       do ip=1,ichoi

          if(lnew) then
             count = 1
             do j=1,iunit
                if(j.eq.first_bead.or.j.eq.second_bead) cycle
                rxp(count,ip) = rxnew(j)
                ryp(count,ip) = rynew(j)
                rzp(count,ip) = rznew(j)
                count = count + 1
             end do
             ax1(1) = rxnew(first_bead)
             ax1(2) = rynew(first_bead)
             ax1(3) = rznew(first_bead)
             ax2(1) = rxnew(second_bead)
             ax2(2) = rynew(second_bead)
             ax2(3) = rznew(second_bead)
          else
             count = 1
             do j=1,iunit
                if(j.eq.first_bead.or.j.eq.second_bead) cycle
                rxp(count,ip) = rxu(imol,j)
                ryp(count,ip) = ryu(imol,j)
                rzp(count,ip) = rzu(imol,j)
                count = count + 1
             end do
             ax1(1) = rxu(imol,first_bead)
             ax1(2) = ryu(imol,first_bead)
             ax1(3) = rzu(imol,first_bead)
             ax2(1) = rxu(imol,second_bead)
             ax2(2) = ryu(imol,second_bead)
             ax2(3) = rzu(imol,second_bead)
          end if

          if(lnew.or.ip.ne.1) then
             phi = twopi*random(-1)
             call dihedral_rigrot(ax1,ax2,iunit-2,rxp(:,ip),ryp(:,ip),rzp(:,ip),phi,rxp(:,ip),ryp(:,ip),rzp(:,ip))
          end if

       end do
    end if

    count = 1
    do j=1,iunit
       lexist(j) = .false.
       if(j.eq.igrow.or.(mvtp.eq.5.and.j.eq.igrow2)) then
          lexist(j) = .true.
          cycle
       end if
    end do
    do j=igrow+1, iunit
       glist(count) = j
       count = count + 1
    end do
    ntogrow = iunit-igrow

    call boltz( lnew,.false.,ovrlap,iskip,imol,imolty,ibox ,ichoi,igrow,ntogrow,glist,maxlen )

    if (ovrlap) then
       lterm = .true.
       return
    end if

    bsum = 0.0E0_dp
    do ip = 1, ichoi
       bsum = bsum + bfac(ip)
    end do

    wadd = wadd * bsum

    if (lnew) then
       if ( wadd .lt. softlog ) then
          lterm = .true.
          return
       end if

       ! select one position at random ---
       rbf = bsum * random(-1)
       bs = 0.0E0_dp
       do ip = 1, ichoi
          if (.not. lovr(ip) ) then
             bs = bs + bfac(ip)
             if ( rbf .lt. bs ) then
                ! select ip position
                iwalk = ip
                goto 15
             end if
          end if
       end do
       write(io_output,*) 'screwup in rigrot'
       call err_exit(__FILE__,__LINE__,'screwup in rigrot',myid+1)
15     continue
    else
       iwalk = 1
       if ( wadd .lt. softlog ) then
          write(io_output,*) '###old rigrot weight too low'
       end if
    end if

    if ( lnew ) then
       vnew(ivTot) = vnew(ivTot) + vtr(ivTot,iwalk)
       vnew(ivExt) = vnew(ivExt) + vtr(ivExt,iwalk)
       vnew(ivInterLJ) = vnew(ivInterLJ) + vtr(ivInterLJ,iwalk)
       vnew(ivIntraLJ) = vnew(ivIntraLJ) + vtr(ivIntraLJ,iwalk)
       vnew(ivElect) = vnew(ivElect) + vtr(ivElect,iwalk)
       vnew(ivEwald) = vnew(ivEwald) + vtr(ivEwald,iwalk)
       vipswn = vipswn+vipswnt(iwalk)
       vwellipswn = vwellipswn+vwellipswnt(iwalk)
    else
       vold(ivTot) = vold(ivTot) + vtr(ivTot,iwalk)
       vold(ivExt) = vold(ivExt) + vtr(ivExt,iwalk)
       vold(ivInterLJ) = vold(ivInterLJ) + vtr(ivInterLJ,iwalk)
       vold(ivIntraLJ) = vold(ivIntraLJ) + vtr(ivIntraLJ,iwalk)
       vold(ivElect) = vold(ivElect) + vtr(ivElect,iwalk)
       vold(ivEwald) = vold(ivEwald) + vtr(ivEwald,iwalk)
       vipswo = vipswo+vipswot(iwalk)
       vwellipswo = vwellipswo+vwellipswot(iwalk)
    end if

    !write(io_output,*) 'vtry', vtr(ivTot,iwalk),iwalk

    do j = 1, ntogrow
       iu = glist(j)
       if (lnew) then
          rxnew(iu) = rxp(j,iwalk)
          rynew(iu) = ryp(j,iwalk)
          rznew(iu) = rzp(j,iwalk)
       end if
    end do

#ifdef __DEBUG__
    write(io_output,*) 'end RIGROT in ',myid
#endif

    return
  end subroutine rigrot

!> \brief Adds H-atoms to a linear carbon chain
!>
!> Potential for methyl-group rotation
!> is: V = 0.5*E0*(1-cos(3*alpha)), where
!> alpha is Ryckaert torsion angle!
  subroutine explct(ichain,vmethyl,lcrysl,lswitch)
      integer::ichain,nngrow,negrow,i,iplus,imins,nn, imolty,iuend,ii,jj
      real::vmethyl,ch,cc,cch,ca,ah,hch2,hk,ck ,en0,aa1,b1,c1,dln1,a2,b2,c2,dln2,a3,b3,c3,x12,y12,z12,x32,y32,z32,xa,ya,za,r&
       ,rx,ry,rz,dr,ven,prob,a4,b4,c4,rn,hch ,ce,ratio
      real::oa,hoh,hoh2,oh,ok,om
      logical::lcrysl,londone,lswitch,lalkanol

! find intramolecular structure
      imolty = moltyp(ichain)
      nngrow = nugrow(imolty)
      negrow = nunit(imolty) - 3
      vmethyl = 0.0E0_dp
      if ((nunit(imolty)-6) .eq. 3*(nngrow-3) ) then
! alkanol cases
         lalkanol = .true.
         iuend = 2
      else
         lalkanol = .false.
         iuend = 1
      end if

      if ( nunit(imolty) .eq. 3 .or. nunit(imolty) .eq. 4) then

! Water case

         if ( nngrow .eq. 3 ) then

! TIP-4P geometry with 2H and O as the growing unit
! M site has been determined by their positions through geometry
! constraint

            om = brvib(itvib(imolty,4,1))
            aa1 = 0.5E0_dp*(rxu(ichain,2)+rxu(ichain,3))  - rxu(ichain,1)
            b1 = 0.5E0_dp*(ryu(ichain,2)+ryu(ichain,3))  - ryu(ichain,1)
            c1 = 0.5E0_dp*(rzu(ichain,2)+rzu(ichain,3))  - rzu(ichain,1)
            dr = sqrt(aa1*aa1+b1*b1+c1*c1)
            rxu(ichain,4) = rxu(ichain,1) + om*aa1/dr
            ryu(ichain,4) = ryu(ichain,1) + om*b1/dr
            rzu(ichain,4) = rzu(ichain,1) + om*c1/dr
            return

         else if ( nngrow .eq. 1 ) then

            if ( itvib(imolty,nngrow+1,1) .eq.  itvib(imolty,nngrow+2,1) ) then
               oh = brvib(itvib(imolty,nngrow+1,1))
               if ( nunit(imolty) .eq. 3 ) then
! SPC geometry

                  hoh = brben(itben(imolty,nngrow+1,1))
               else
! TIP-4P geometry

                  hoh = brben(itben(imolty,nngrow+1,1))
               end if
               if ( lcrysl ) then
                  hoh2 = hoh / 2.0E0_dp
                  hk = oh*sin(hoh2)
                  ok = oh*cos(hoh2)
                  rxu(ichain,2) = rxu(ichain,1)
                  ryu(ichain,2) = ryu(ichain,1) + hk
                  rzu(ichain,2) = rzu(ichain,1) + ok
                  rxu(ichain,3) = rxu(ichain,1)
                  ryu(ichain,3) = ryu(ichain,1) - hk
                  rzu(ichain,3) = rzu(ichain,1) + ok
               else
                  oa = oh*cos(onepi-hoh)
                  ah = oh*sin(onepi-hoh)
! generate a random vector on a sphere for the first H ---
 111              rx = 2.0E0_dp*random(-1) - 1.0E0_dp
                  ry = 2.0E0_dp*random(-1) - 1.0E0_dp
                  dr = rx*rx + ry*ry
                  if ( dr .gt. 1.0E0_dp ) goto 111
                  rz = 2.0E0_dp*sqrt(1-dr)
                  a2 = rx*rz
                  b2 = ry*rz
                  c2 = 1 - 2.0E0_dp*dr

                  rxu(ichain,2) = rxu(ichain,1) - oh*a2
                  ryu(ichain,2) = ryu(ichain,1) - oh*b2
                  rzu(ichain,2) = rzu(ichain,1) - oh*c2
! generate another random vector on a sphere for the second H ---
 222              rx = 2.0E0_dp*random(-1) - 1.0E0_dp
                  ry = 2.0E0_dp*random(-1) - 1.0E0_dp
                  dr = rx*rx + ry*ry
                  if ( dr .gt. 1.0E0_dp ) goto 222
                  rz = 2.0E0_dp*sqrt(1-dr)
                  rx = rx*rz
                  ry = ry*rz
                  rz = 1 - 2.0E0_dp*dr

! The two vectors above form a plane identified by n1 ---
                  aa1 = b2*rz - c2*ry
                  b1 = -(a2*rz - rx*c2)
                  c1 = a2*ry - rx*b2
                  dln1 = sqrt(aa1*aa1 + b1*b1 + c1*c1)
                  aa1 = aa1/dln1
                  b1 = b1/dln1
                  c1 = c1/dln1

! cross product n1 x n2 calc.-> n3 ---

                  a3 = b1*c2 - b2*c1
                  b3 = -(aa1*c2 - a2*c1)
                  c3 = aa1*b2 - a2*b1
                  rxu(ichain,3) = rxu(ichain,1) + a2*oa + ah*a3
                  ryu(ichain,3) = ryu(ichain,1) + b2*oa + ah*b3
                  rzu(ichain,3) = rzu(ichain,1) + c2*oa + ah*c3
               end if
               if ( nunit(imolty) .eq. 4 ) then
                  om = brvib(itvib(imolty,4,1))
                  aa1 = 0.5E0_dp*(rxu(ichain,2)+rxu(ichain,3))  - rxu(ichain,1)
                  b1 = 0.5E0_dp*(ryu(ichain,2)+ryu(ichain,3))  - ryu(ichain,1)
                  c1 = 0.5E0_dp*(rzu(ichain,2)+rzu(ichain,3))  - rzu(ichain,1)
                  dr = sqrt(aa1*aa1+b1*b1+c1*c1)
                  rxu(ichain,4) = rxu(ichain,1) + om*aa1/dr
                  ryu(ichain,4) = ryu(ichain,1) + om*b1/dr
                  rzu(ichain,4) = rzu(ichain,1) + om*c1/dr
               end if
               return
            else
! HF model
! generate a random vector on a sphere for the H and M sites ---
 333           rx = 2.0E0_dp*random(-1) - 1.0E0_dp
               ry = 2.0E0_dp*random(-1) - 1.0E0_dp
               dr = rx*rx + ry*ry
               if ( dr .gt. 1.0E0_dp ) goto 333
               rz = 2.0E0_dp*sqrt(1-dr)
               a2 = rx*rz
               b2 = ry*rz
               c2 = 1 - 2.0E0_dp*dr

               do i = 2,3
                  ch = brvib(itvib(imolty,i,1))
                  rxu(ichain,i) = rxu(ichain,1) + ch*a2
                  ryu(ichain,i) = ryu(ichain,1) + ch*b2
                  rzu(ichain,i) = rzu(ichain,1) + ch*c2
               end do
            end if
         end if

      else if ( nunit(imolty) .eq. 5 .or. nunit(imolty) .eq. 7 .or. nunit(imolty) .eq. 9 ) then

! Methane case or other rigid molecule with 5 units

         ch = brvib(itvib(imolty,nngrow+1,1))
         hch = brben(itben(imolty,nngrow+1,1))
         if ( nngrow .eq. 3 ) then

! for hydro-furan

! aa1 = rxu(ichain,3)-rxu(ichain,1)
! b1 = ryu(ichain,3)-ryu(ichain,1)
! c1 = rzu(ichain,3)-rzu(ichain,1)
! dln1 = sqrt(aa1*aa1 + b1*b1 + c1*c1)
! aa1 = aa1/dln1
! b1 = b1/dln1
! c1 = c1/dln1
! a2 = 0.5E0_dp*(rxu(ichain,1)+rxu(ichain,3))-
!     &           rxu(ichain,2)
! b2 = 0.5E0_dp*(ryu(ichain,1)+ryu(ichain,3))-
!     &           ryu(ichain,2)
! c2 = 0.5E0_dp*(rzu(ichain,1)+rzu(ichain,3))-
!     &           rzu(ichain,2)
! dln1 = sqrt(a2*a2 + b2*b2 + c2*c2)
! a2 = a2/dln1
! b2 = b2/dln1
! c2 = c2/dln1
! rxu(ichain,4) = rxu(ichain,2) + 2.2758059E0_dp*a2 +
!     &           0.765E0_dp*aa1
! ryu(ichain,4) = ryu(ichain,2) + 2.2758059E0_dp*b2 +
!     &           0.765E0_dp*b1
! rzu(ichain,4) = rzu(ichain,2) + 2.2758059E0_dp*c2 +
!     &           0.765E0_dp*c1
! rxu(ichain,5) = rxu(ichain,4) - 1.53E0_dp*aa1
! ryu(ichain,5) = ryu(ichain,4) - 1.53E0_dp*b1
! rzu(ichain,5) = rzu(ichain,4) - 1.53E0_dp*c1

            aa1 = rxu(ichain,2)-rxu(ichain,1)
            b1 = ryu(ichain,2)-ryu(ichain,1)
            c1 = rzu(ichain,2)-rzu(ichain,1)

            a2 = rxu(ichain,3)-rxu(ichain,2)
            b2 = ryu(ichain,3)-ryu(ichain,2)
            c2 = rzu(ichain,3)-rzu(ichain,2)

! cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(aa1*c2 - a2*c1)
            c3 = aa1*b2 - a2*b1
            dln1 = sqrt(a3*a3 + b3*b3 + c3*c3)
            a2 = a3/dln1
            b2 = b3/dln1
            c2 = c3/dln1

            aa1 = 0.5E0_dp*(rxu(ichain,1)+rxu(ichain,3))- rxu(ichain,2)
            b1 = 0.5E0_dp*(ryu(ichain,1)+ryu(ichain,3))- ryu(ichain,2)
            c1 = 0.5E0_dp*(rzu(ichain,1)+rzu(ichain,3))- rzu(ichain,2)
            dln1 = sqrt(aa1*aa1 + b1*b1 + c1*c1)
            aa1 = aa1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
            a3 = b1*c2 - b2*c1
            b3 = -(aa1*c2 - a2*c1)
            c3 = aa1*b2 - a2*b1

            rxu(ichain,4) = rxu(ichain,2) + 2.25E0_dp*aa1 +  0.41E0_dp*a2+0.77E0_dp*a3
            ryu(ichain,4) = ryu(ichain,2) + 2.25E0_dp*b1 +  0.41E0_dp*b2+0.77E0_dp*b3
            rzu(ichain,4) = rzu(ichain,2) + 2.25E0_dp*c1 +  0.41E0_dp*c2+0.77E0_dp*c3
            rxu(ichain,5) = rxu(ichain,4) - 1.54E0_dp*a3
            ryu(ichain,5) = ryu(ichain,4) - 1.54E0_dp*b3
            rzu(ichain,5) = rzu(ichain,4) - 1.54E0_dp*c3
         else if ( nugrow(imolty) .eq. 2) then

            ca = ch*cos(onepi-hch)
            ah = ch*sin(onepi-hch)

! only one H participate in the growing
! this subroutine put on the rest 3 hydrogens
! H-atoms for the first CH3-group---
! first define vector C-H ------
            x12 = rxu(ichain,2)-rxu(ichain,1)
            y12 = ryu(ichain,2)-ryu(ichain,1)
            z12 = rzu(ichain,2)-rzu(ichain,1)
! generate a random vector ------
 444        rx = 2.0E0_dp*random(-1) - 1.0E0_dp
            ry = 2.0E0_dp*random(-1) - 1.0E0_dp
            dr = rx*rx + ry*ry
            if ( dr .gt. 1.0E0_dp ) goto 444
            rz = 2.0E0_dp*sqrt(1-dr)
            rx = rx*rz
            ry = ry*rz
            rz = 1 - 2.0E0_dp*dr

! the two vectors above form a plane identified by n1 ------
            aa1 = y12*rz - z12*ry
            b1 = -(x12*rz - rx*z12)
            c1 = x12*ry - rx*y12
            dln1 = sqrt(aa1**2 + b1**2 + c1**2)
            aa1 = aa1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
! normalizing C-H vector -> n2 ----
            a2 = x12/ch
            b2 = y12/ch
            c2 = z12/ch
! cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(aa1*c2 - a2*c1)
            c3 = aa1*b2 - a2*b1
! point A ------
            xa = rxu(ichain,2)+a2*ca
            ya = ryu(ichain,2)+b2*ca
            za = rzu(ichain,2)+c2*ca
! H-atoms for the first methyl group------
            r = 0.0E0_dp
            do i = 1,3
               a4 = ah*a3*cos(r) + ah*aa1*sin(r)
               b4 = ah*b3*cos(r) + ah*b1*sin(r)
               c4 = ah*c3*cos(r) + ah*c1*sin(r)
               rxu(ichain,nngrow+i) = xa + a4
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0E0_dp*onepi/180.0E0_dp
            end do

         else
! only carbon has been grown

! METHANE ------
            if ( lcrysl ) then
               hch2 = brben(itben(imolty,2,1))/2
               hk = ch*sin(hch2)
               ck = ch*cos(hch2)
               rxu(ichain,2) = rxu(ichain,1) - ck
               ryu(ichain,2) = ryu(ichain,1)
               rzu(ichain,2) = rzu(ichain,1) + hk
               rxu(ichain,3) = rxu(ichain,1) - ck
               ryu(ichain,3) = ryu(ichain,1)
               rzu(ichain,3) = rzu(ichain,1) - hk
               hch2 = brben(itben(imolty,4,1))/2
               ch = brvib(itvib(imolty,4,1))
               hk = ch*sin(hch2)
               ck = ch*cos(hch2)
               rxu(ichain,4) = rxu(ichain,1) + ck
               ryu(ichain,4) = ryu(ichain,1) + hk
               rzu(ichain,4) = rzu(ichain,1)
               rxu(ichain,5) = rxu(ichain,1) + ck
               ryu(ichain,5) = ryu(ichain,1) - hk
               rzu(ichain,5) = rzu(ichain,1)

! KEEP THE OLD CONFIGURATION in SWITCH MOVE

            else if ( lswitch ) then
! if (imolty .eq. 1) then
! ratio = brvib(itvib(1,nngrow+1,1))
!     &                 /brvib(itvib(2,nngrow+1,1))
! else
! ratio = brvib(itvib(2,nngrow+1,1))
!     &                 /brvib(itvib(1,nngrow+1,1))
! end if
               ratio = 1.0E0_dp
               rxu(ichain,2) = rxu(ichain,1) + ratio*(rxu(ichain,2) -rxu(ichain,1))
               ryu(ichain,2) = ryu(ichain,1) + ratio*(ryu(ichain,2) -ryu(ichain,1))
               rzu(ichain,2) = rzu(ichain,1) + ratio*(rzu(ichain,2) -rzu(ichain,1))
               rxu(ichain,3) = rxu(ichain,1) + ratio*(rxu(ichain,3) -rxu(ichain,1))
               ryu(ichain,3) = ryu(ichain,1) + ratio*(ryu(ichain,3) -ryu(ichain,1))
               rzu(ichain,3) = rzu(ichain,1) + ratio*(rzu(ichain,3) -rzu(ichain,1))
               ratio = 0.7E0_dp/0.25E0_dp
               rxu(ichain,4) = rxu(ichain,1) + ratio*(rxu(ichain,4) -rxu(ichain,1))
               ryu(ichain,4) = ryu(ichain,1) + ratio*(ryu(ichain,4) -ryu(ichain,1))
               rzu(ichain,4) = rzu(ichain,1) + ratio*(rzu(ichain,4) -rzu(ichain,1))
               rxu(ichain,5) = rxu(ichain,1) + ratio*(rxu(ichain,5) -rxu(ichain,1))
               ryu(ichain,5) = ryu(ichain,1) + ratio*(ryu(ichain,5) -ryu(ichain,1))
               rzu(ichain,5) = rzu(ichain,1) + ratio*(rzu(ichain,5) -rzu(ichain,1))
            else
               ch = brvib(itvib(imolty,2,1))
               hch2 = brben(itben(imolty,2,1))/2
               ca = ch*cos(hch2)
               ah = ch*sin(hch2)
! generate a random vector for the first H ------
 555           rx = 2.0E0_dp*random(-1) - 1.0E0_dp
               ry = 2.0E0_dp*random(-1) - 1.0E0_dp
               dr = rx*rx + ry*ry
               if ( dr .gt. 1.0E0_dp ) goto 555
               rz = 2.0E0_dp*sqrt(1-dr)
               a2 = rx*rz
               b2 = ry*rz
               c2 = 1 - 2.0E0_dp*dr
! generate another random vector for the rotation ------
 666           rx = 2.0E0_dp*random(-1) - 1.0E0_dp
               ry = 2.0E0_dp*random(-1) - 1.0E0_dp
               dr = rx*rx + ry*ry
               if ( dr .gt. 1.0E0_dp ) goto 666
               rz = 2.0E0_dp*sqrt(1-dr)
               rx = rx*rz
               ry = ry*rz
               rz = 1 - 2.0E0_dp*dr
! the two vectors above form a plane identified by n1 ------
               aa1 = b2*rz - c2*ry
               b1 = -(a2*rz - rx*c2)
               c1 = a2*ry - rx*b2
               dln1 = sqrt(aa1**2 + b1**2 + c1**2)
               aa1 = aa1/dln1
               b1 = b1/dln1
               c1 = c1/dln1

               rxu(ichain,2) = rxu(ichain,1) + a2*ca + aa1*ah
               ryu(ichain,2) = ryu(ichain,1) + b2*ca + b1*ah
               rzu(ichain,2) = rzu(ichain,1) + c2*ca + c1*ah
               rxu(ichain,3) = rxu(ichain,1) + a2*ca - aa1*ah
               ryu(ichain,3) = ryu(ichain,1) + b2*ca - b1*ah
               rzu(ichain,3) = rzu(ichain,1) + c2*ca - c1*ah

               ch = brvib(itvib(imolty,4,1))
               hch2 = brben(itben(imolty,4,1))/2
               ca = ch*cos(hch2)
               ah = ch*sin(hch2)

! cross product n1 x n2 calc.-> n3 ---
               a3 = b1*c2 - b2*c1
               b3 = -(aa1*c2 - a2*c1)
               c3 = aa1*b2 - a2*b1

               rxu(ichain,4) = rxu(ichain,1) - a2*ca + a3*ah
               ryu(ichain,4) = ryu(ichain,1) - b2*ca + b3*ah
               rzu(ichain,4) = rzu(ichain,1) - c2*ca + c3*ah
               rxu(ichain,5) = rxu(ichain,1) - a2*ca - a3*ah
               ryu(ichain,5) = ryu(ichain,1) - b2*ca - b3*ah
               rzu(ichain,5) = rzu(ichain,1) - c2*ca - c3*ah

            end if
            if ( nunit(imolty) .eq. 9 ) then
               ce = brvib(itvib(imolty,6,1))
               ratio = ce / ch
               do i = 1,4
                  ii = i + 1
                  jj = i + 5
                  rxu(ichain,jj) = rxu(ichain,1) + ratio* (rxu(ichain,ii)-rxu(ichain,1))
                  ryu(ichain,jj) = ryu(ichain,1) + ratio* (ryu(ichain,ii)-ryu(ichain,1))
                  rzu(ichain,jj) = rzu(ichain,1) + ratio* (rzu(ichain,ii)-rzu(ichain,1))
               end do
            end if
         end if
         if ( nunit(imolty) .eq. 7 ) then
! for seven site water model
            rxu(ichain,6) =  4E0_dp*rxu(ichain,4)-3E0_dp*rxu(ichain,1)
            ryu(ichain,6) =  4E0_dp*ryu(ichain,4)-3E0_dp*ryu(ichain,1)
            rzu(ichain,6) =  4E0_dp*rzu(ichain,4)-3E0_dp*rzu(ichain,1)
            rxu(ichain,7) =  4E0_dp*rxu(ichain,5)-3E0_dp*rxu(ichain,1)
            ryu(ichain,7) =  4E0_dp*ryu(ichain,5)-3E0_dp*ryu(ichain,1)
            rzu(ichain,7) =  4E0_dp*rzu(ichain,5)-3E0_dp*rzu(ichain,1)
         end if
      else

! ALKANES (not for methane)
! find intramolecular structure for longer alkanes
! WARNING: work only for pure hydro- or perfluoro-carbons!!!

         cc = brvib(itvib(imolty,1,1))
         ch = brvib(itvib(imolty,nngrow+1,1))

         cch = brben(itben(imolty,nngrow+1,1))
         ca = ch*cos(onepi-cch)
         ah = ch*sin(onepi-cch)
         if ( nngrow .gt. 2) then
            hch = brben(itben(imolty,nngrow+4,1))
            hch2 = hch/2.0E0_dp
            hk = ch*sin(hch2)
            ck = ch*cos(hch2)
            en0 = 853.93E0_dp

! REGULAR N-ALKANE ------
! main loop
! calculates hydrogen positions for methylene groups
! WARNING only works for linear alkanes
            do i = 2,nngrow-iuend
               iplus = i+1
               imins = i-1
               aa1 = (ryu(ichain,imins)-ryu(ichain,i))*(rzu(ichain ,iplus)-rzu(ichain,i))&
                - (rzu(ichain,imins)-rzu(ichain ,i))*(ryu(ichain,iplus)-ryu(ichain,i))
               b1 = -(rxu(ichain,imins)-rxu(ichain,i))*(rzu(ichain ,iplus)-rzu(ichain,i))&
                + (rzu(ichain,imins)-rzu(ichain ,i))*(rxu(ichain,iplus)-rxu(ichain,i))
               c1 = (rxu(ichain,imins)-rxu(ichain,i))*(ryu(ichain,iplus) -ryu(ichain,i))&
                - (ryu(ichain,imins)-ryu(ichain,i))* (rxu(ichain,iplus)-rxu(ichain,i))
               dln1 = sqrt(aa1**2 + b1**2 + c1**2)
               aa1 = aa1/dln1
               b1 = b1/dln1
               c1 = c1/dln1
               a2 = (rxu(ichain,iplus)-rxu(ichain,imins))
               b2 = (ryu(ichain,iplus)-ryu(ichain,imins))
               c2 = (rzu(ichain,iplus)-rzu(ichain,imins))
               dln2 = sqrt(a2**2 + b2**2 + c2**2)
               a2 = a2/dln2
               b2 = b2/dln2
               c2 = c2/dln2
               a3 = b1*c2 - c1*b2
               b3 = -aa1*c2 + a2*c1
               c3 = aa1*b2 - b1*a2
               londone = .false.
               do nn = nngrow+4+2*(i-2), nngrow+5+2*(i-2)
                  if (.not.londone) then
                     rxu(ichain,nn) = rxu(ichain,i) + a3*ck + aa1*hk
                     ryu(ichain,nn) = ryu(ichain,i) + b3*ck + b1*hk
                     rzu(ichain,nn) = rzu(ichain,i) + c3*ck + c1*hk
                     londone = .true.
                  else
                     rxu(ichain,nn) = rxu(ichain,i) + a3*ck - aa1*hk
                     ryu(ichain,nn) = ryu(ichain,i) + b3*ck - b1*hk
                     rzu(ichain,nn) = rzu(ichain,i) + c3*ck - c1*hk
                  end if
               end do
            end do

! H-atoms for the first CH3-group---
! define c1c2c3 plane -> n1 ------
            x12 = rxu(ichain,1)-rxu(ichain,2)
            y12 = ryu(ichain,1)-ryu(ichain,2)
            z12 = rzu(ichain,1)-rzu(ichain,2)
            x32 = rxu(ichain,3)-rxu(ichain,2)
            y32 = ryu(ichain,3)-ryu(ichain,2)
            z32 = rzu(ichain,3)-rzu(ichain,2)
            aa1 = y12*z32 - z12*y32
            b1 = -(x12*z32 - x32*z12)
            c1 = x12*y32 - x32*y12
            dln1 = sqrt(aa1**2 + b1**2 + c1**2)
            aa1 = aa1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
! normalizing c2c1 vector -> n2 ----
            a2 = x12/cc
            b2 = y12/cc
            c2 = z12/cc
! cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(aa1*c2 - a2*c1)
            c3 = aa1*b2 - a2*b1
! point A ------
            xa = rxu(ichain,1)+a2*ca
            ya = ryu(ichain,1)+b2*ca
            za = rzu(ichain,1)+c2*ca
! H-atoms ------
            if ( lcrysl ) then
               r=0.0E0_dp
               ven = 0.0E0_dp
            else
 101           r = 2.0E0_dp*onepi*random(-1)
               ven = en0*(1.0E0_dp-cos(3.0E0_dp*r))
               prob = exp(-beta*ven)
               rn = random(-1)
               if (rn.gt.prob) goto 101
            end if
! rgr = r*180.0E0_dp/onepi
! write(io_output,*) 'final Ryckaert angle: ',rgr
            vmethyl = vmethyl + ven
            do i = 1,3
               a4 = ah*a3*cos(r+onepi) + ah*aa1*sin(r+onepi)
               b4 = ah*b3*cos(r+onepi) + ah*b1*sin(r+onepi)
               c4 = ah*c3*cos(r+onepi) + ah*c1*sin(r+onepi)

               rxu(ichain,nngrow+i) = xa + a4
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0E0_dp*onepi/180.0E0_dp
            end do

! if it is an alkanol molecule, it did not have the other ending CH3

            if ( lalkanol ) return

! H-atoms for the last CH3------
! define c1c2c3 plane -> n1 ------
            x12 = rxu(ichain,nngrow)-rxu(ichain,nngrow-1)
            y12 = ryu(ichain,nngrow)-ryu(ichain,nngrow-1)
            z12 = rzu(ichain,nngrow)-rzu(ichain,nngrow-1)
            x32 = rxu(ichain,nngrow-2)-rxu(ichain,nngrow-1)
            y32 = ryu(ichain,nngrow-2)-ryu(ichain,nngrow-1)
            z32 = rzu(ichain,nngrow-2)-rzu(ichain,nngrow-1)
            aa1 = y12*z32 - z12*y32
            b1 = -(x12*z32 - x32*z12)
            c1 = x12*y32 - x32*y12
            dln1 = sqrt(aa1**2 + b1**2 + c1**2)
            aa1 = aa1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
! normalizing c2c1 vector -> n2 ----
            a2 = x12/cc
            b2 = y12/cc
            c2 = z12/cc
! cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(aa1*c2 - a2*c1)
            c3 = aa1*b2 - a2*b1
! point A ------
            xa = rxu(ichain,nngrow)+a2*ca
            ya = ryu(ichain,nngrow)+b2*ca
            za = rzu(ichain,nngrow)+c2*ca
! H-atoms -------
            if (lcrysl) then
               r = 0.0E0_dp
               ven = 0.0E0_dp
            else
 200           r = 2.0E0_dp*onepi*random(-1)
               ven = en0*(1.0E0_dp-cos(3.0E0_dp*r))
               prob = exp(-beta*ven)
               rn = random(-1)
               if (rn.gt.prob) goto 200
            end if
            vmethyl = vmethyl + ven
            do i = 1,3
               a4 = ah*a3*cos(r+onepi) + ah*aa1*sin(r+onepi)
               b4 = ah*b3*cos(r+onepi) + ah*b1*sin(r+onepi)
               c4 = ah*c3*cos(r+onepi) + ah*c1*sin(r+onepi)
               rxu(ichain,negrow+i) = xa + a4
               ryu(ichain,negrow+i) = ya + b4
               rzu(ichain,negrow+i) = za + c4
               r = r + 120.0E0_dp*onepi/180.0E0_dp
            end do

         else if (nngrow .eq. 2 .and. nunit(imolty) .eq. 8) then
! ETHANE  -----
! define the new type en0 ------
! use en0 for torsion type 19
            en0 = 716.77E0_dp
! for ethane case ------
! H-atoms for the first CH3-group---
! first define vector C2C1 ------
            x12 = rxu(ichain,1)-rxu(ichain,2)
            y12 = ryu(ichain,1)-ryu(ichain,2)
            z12 = rzu(ichain,1)-rzu(ichain,2)
! generate a random vector ------
 777        rx = 2.0E0_dp*random(-1) - 1.0E0_dp
            ry = 2.0E0_dp*random(-1) - 1.0E0_dp
            dr = rx*rx + ry*ry
            if ( dr .gt. 1.0E0_dp ) goto 777
            rz = 2.0E0_dp*sqrt(1-dr)
            rx = rx*rz
            ry = ry*rz
            rz = 1 - 2.0E0_dp*dr

! the two vectors above form a plane identified by n1 ------
            aa1 = y12*rz - z12*ry
            b1 = -(x12*rz - rx*z12)
            c1 = x12*ry - rx*y12
            dln1 = sqrt(aa1**2 + b1**2 + c1**2)
            aa1 = aa1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
! normalizing c2c1 vector -> n2 ----
            a2 = x12/cc
            b2 = y12/cc
            c2 = z12/cc
! cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(aa1*c2 - a2*c1)
            c3 = aa1*b2 - a2*b1
! point A ------
            xa = rxu(ichain,1)+a2*ca
            ya = ryu(ichain,1)+b2*ca
            za = rzu(ichain,1)+c2*ca
! H-atoms for the first methyl group------
            r = 0.0E0_dp
            do i = 1,3
               a4 = ah*a3*cos(r) + ah*aa1*sin(r)
               b4 = ah*b3*cos(r) + ah*b1*sin(r)
               c4 = ah*c3*cos(r) + ah*c1*sin(r)
               rxu(ichain,nngrow+i) = xa + a4
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0E0_dp*onepi/180.0E0_dp
            end do
! H-atoms for the second methyl group ------
! point A' is opposite to point A ------
            xa = rxu(ichain,2)-a2*ca
            ya = ryu(ichain,2)-b2*ca
            za = rzu(ichain,2)-c2*ca
            if ( lcrysl ) then
               r=0.0E0_dp
               ven = 0.0E0_dp
            else
 100           r = 2.0E0_dp*onepi*random(-1)
               ven = en0*(1.0E0_dp-cos(3.0E0_dp*r))
               prob = exp(-beta*ven)
               rn = random(-1)
               if (rn.gt.prob) goto 100
            end if
! rgr = r*180.0E0_dp/onepi
! write(io_output,*) 'final Ryckaert angle: ',rgr
            vmethyl = vmethyl + ven

            do i = 4,6
               a4 = ah*a3*cos(r+onepi) + ah*aa1*sin(r+onepi)
               b4 = ah*b3*cos(r+onepi) + ah*b1*sin(r+onepi)
               c4 = ah*c3*cos(r+onepi) + ah*c1*sin(r+onepi)
               rxu(ichain,nngrow+i) = xa + a4
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0E0_dp*onepi/180.0E0_dp
            end do

         end if
      end if
      return
  end subroutine explct

!********************************************************
!> \brief Places hydrogens after the growth of the backbone of
!> a molecule for linear, branched or cylic molecules.
!> Uses CDCBMC to grow them
!********************************************************
  subroutine place(lnew,lterm,i,imolty,ibox,index,wplace)
    use util_math,only:cone_angle
    ! integer,parameter::max=numax
    logical::lnew,lterm,ovrlap

    integer::i,j,imolty,count,counta,iu,ju,ku,jtvib ,start,iv,index,ivib,nchvib,ibend,ib,type,ip,ichoi,niplace,iw,iufrom,it,jut2&
     ,jut3,jut4,ibox,glist,iwalk,iuprev,list,nchben_a,nchben_b,iu2back

    real::wplace,equil,kforce,bsum_try,mincb ,delcb,ux,uy,uz,r,vvib,bfactor,third,length,bs,rbf,vvibtr ,wei_vib,bendang,vangle&
     ,vbbtr,angle,vphi,thetac,rx,ry,rz ,rsint,dist,wei_bend,ang_trial,vdha,vbend,vtorsion ,bsum,alpha,gamma,dum,phi,thetatwo,phitwo

    dimension r(nchbn_max),bfactor(nchbn_max),bendang(numax,numax),ang_trial(nchbn_max),dist(numax),niplace(numax),vbend(nchmax)&
     ,vtorsion(nchmax),phi(numax),list(numax),glist(numax)

#ifdef __DEBUG__
    write(io_output,*) 'START PLACE in ',myid
#endif

      nchvib = nchbna(imolty)
      nchben_a = nchvib
      nchben_b = nchbnb(imolty)
      third = 1.0E0_dp / 3.0E0_dp
      wplace = 1.0E0_dp

      ichoi = nchoi(imolty)

      do j = 1, nunit(imolty)
         niplace(j) = 0
      end do

      do iw = 1, nplace
         do count = 1, pnum(iw)
            iu = iplace(iw,count)
            niplace(iu) = iw
         end do
      end do

      do iw = 1, nplace
         vvibtr = 0.0E0_dp
         wei_vib = 1.0E0_dp
         iufrom = pfrom(iw)
         do count = 1, pnum(iw)

            iu = iplace(iw,count)

            if (invib(imolty,iu).gt.1) then
               write(io_output,*) 'iu,invib',iu,invib(imolty,iu)
               call err_exit(__FILE__,__LINE__,'invib can no be larger than one for hydrogen',myid+1)
            end if

! determine bond lengths
            iv = 1

            ju = ijvib(imolty,iu,iv)

            if (iufrom.ne.ju) then
               write(io_output,*) 'iu,ju,iufrom',iu,ju,iufrom
               call err_exit(__FILE__,__LINE__,'ju not equal to iufrom',myid+1)
            end if

            jtvib = itvib(imolty,iu,iv)

            equil = brvib(jtvib)
            kforce = brvibk(jtvib)

            if (kforce.gt.1.0E-3_dp) then
! we will use flexible bond lengths
               bsum_try = 0.0E0_dp
               mincb = brvibmin(jtvib)**3
               delcb = brvibmax(jtvib)**3 - mincb

               if (.not. lnew) then
                  ux = rxu(i,ju) - rxu(i,iu)
                  uy = ryu(i,ju) - ryu(i,iu)
                  uz = rzu(i,ju) - rzu(i,iu)
                  r(1) = sqrt(ux**2 + uy**2 + uz**2)
                  if (L_vib_table) then
                     vvib = lininter_vib(r(1),jtvib)
                   else
                     vvib = kforce * (r(1) - equil)**2
                  end if
                  bfactor(1) = exp(-beta*vvib)
                  bsum_try = bsum_try + bfactor(1)
                  start = 2
               else
                  start = 1
               end if

               do ivib = start, nchvib
                  r(ivib) = (mincb + random(-1)*delcb)**third
                  if (L_vib_table) then
                     vvib = lininter_vib(r(ivib),jtvib)
                  else
                     vvib = kforce * ( r(ivib) - equil )**2
                  end if
                  bfactor(ivib) = exp(-beta*vvib)
                  bsum_try = bsum_try + bfactor(ivib)
               end do

               wei_vib = wei_vib * bsum_try

               if (lnew) then
! select one of the trial sites via bias
                  rbf = random(-1)*bsum_try
                  bs = 0.0E0_dp
                  do ivib = 1, nchvib
                     bs = bs + bfactor(ivib)
                     if (rbf .lt. bs ) then
                        length = r(ivib)
                        vvib = log(bfactor(ivib))/(-beta)
                        goto 5
                     end if
                  end do
 5                continue
               else
! select old conformation
                  length = r(1)
                  vvib = log(bfactor(1))/(-beta)
               end if
               vvibtr = vvibtr + vvib

            else

! our bond lengths are fixed
               length = equil

            end if

            distij(iu,ju) = length
            distij(ju,iu) = length
         end do
         if (lnew) then
            vnew(ivTot) = vnew(ivTot) + vvibtr
            vnew(ivStretching)  = vnew(ivStretching)  + vvibtr
         else
            vold(ivTot) = vold(ivTot) + vvibtr
            vold(ivStretching) = vold(ivStretching) + vvibtr
         end if

      end do

      do iw = 1, nplace

         iufrom = pfrom(iw)
         iuprev = pprev(iw)

! first, set up cone
         dist(2) = distij(iuprev,iufrom)

         rx = xvec(iuprev,iufrom) / dist(2)
         ry = yvec(iuprev,iufrom) / dist(2)
         rz = zvec(iuprev,iufrom) / dist(2)

         call cone(1,rx,ry,rz,dum,dum)

! now that we set up cone, we must determine other beads grown
! from iufrom
         counta = 2

         do iv = 1, invib(imolty,iufrom)
            ku = ijvib(imolty,iufrom,iv)

! make sure that ku is not equal to a site we are growing or iuprev
            if (ku.eq.iuprev) goto 95
            do count = 1, pnum(iw)
               iu = iplace(iw,count)
               if (iu.eq.ku) goto 95
            end do

! we must determine the angle associated with this one
            counta = counta + 1
            dist(counta) = distij(iufrom,ku)
            ux = xvec(iufrom,ku) / dist(counta)
            uy = yvec(iufrom,ku) / dist(counta)
            uz = zvec(iufrom,ku) / dist(counta)

! determine angle with iuprev
            thetac = -(ux*rx + uy*ry + uz*rz)
            bendang(ku,iuprev) = acos(thetac)

            alpha = bendang(ku,iuprev)

            call cone(3,ux,uy,uz,alpha,gamma)

            phi(counta) = gamma
            list(counta) = ku

 95         continue
         end do

         do ip = 1, ichoi

            wei_bend = 1.0E0_dp
            vbbtr = 0.0E0_dp
            vdha = 0.0E0_dp

            do count = 1, pnum(iw)

               iu = iplace(iw,count)

               ju = ijvib(imolty,iu,1)

               do ib = 1, inben(imolty,iu)

                  ku = ijben3(imolty,iu,ib)

                  if (ku.eq.pprev(iw)) then

                     if (ju.ne.ijben2(imolty,iu,ib)) then
                        write(io_output,*) 'ju,ijben2',ju,ijben2(imolty,iu,ib)
                        call err_exit(__FILE__,__LINE__,'ju not equal to ijben2 in place',myid+1)
                     end if

                     type = itben(imolty,iu,ib)
                     equil = brben(type)
                     kforce = brbenk(type)

! initialize bsum_try
                     bsum_try = 0

                     if (.not.lnew.and.ip.eq.1) then
! first ibend is the old conformation
! compute vector from iufrom to iugrow
                        ux = rxu(i,iu) - rxu(i,ju)
                        uy = ryu(i,iu) - ryu(i,ju)
                        uz = rzu(i,iu) - rzu(i,ju)
                        dist(1) = sqrt(ux**2 + uy**2 + uz**2)

! dot product divided by lengths gives cos(angle)
                        thetac = -( ux*rx + uy*ry  + uz*rz )  / (dist(1))
                        angle = acos(thetac)

! compute the energy of this angle
                        vangle = kforce * (angle - equil)**2
                        ang_trial(1) = angle
                        bfactor(1) = exp( -beta*vangle )
                        bsum_try = bsum_try + bfactor(1)

! skip first ibend in next loop
                        start = 2

                     else
! new conformation start at 1
                        start = 1
                     end if

! compute trial angles and energies
                     do ibend = start,nchben_a
! choose the angle uniformly on sin(theta)
                        rsint = 2.0E0_dp*random(-1) - 1.0E0_dp
                        angle = acos(rsint)
                        ang_trial(ibend) = angle

! calculate the bond angle energy
                        vangle = kforce * (angle - equil)**2
                        bfactor(ibend) = exp(-beta*vangle)
                        bsum_try = bsum_try + bfactor(ibend)
                     end do

                     if ( lnew.or.ip.ne.1 ) then
! select one of the trial sites via bias
                        rbf = random(-1)*bsum_try
                        bs = 0.0E0_dp
                        do ibend = 1,nchben_a
                           bs = bs + bfactor(ibend)
                           if ( rbf .lt. bs ) then
                              angle = ang_trial(ibend)
                              vangle = log(bfactor(ibend))/(-beta)
                              goto 10
                           end if
                        end do
 10                     continue
                     else
! select the old conformation
                        angle = ang_trial(1)
                        vangle = log(bfactor(1))/(-beta)
                     end if

! propagate the rosenbluth weight
                     wei_bend = wei_bend * bsum_try/dble(nchben_a)

                     bendang(iu,ku) = angle

                     vbbtr = vbbtr + vangle
                  end if

               end do

            end do

! now we must determine the second posible position for our sites

            do count = 1, pnum(iw)

               iu = iplace(iw,count)

! initialize bsum_try
               bsum_try = 0.0E0_dp

               if ( .not. lnew .and.ip.eq.1) then
                  ux = rxu(i,iu) - rxu(i,iufrom)
                  uy = ryu(i,iu) - ryu(i,iufrom)
                  uz = rzu(i,iu) - rzu(i,iufrom)

                  dist(1) = sqrt(ux**2 + uy**2 + uz**2)

                  ux = ux / dist(1)
                  uy = uy / dist(1)
                  uz = uz / dist(1)
                  thetatwo = bendang(iu,iuprev)

                  call cone(3,ux,uy,uz,thetatwo,phitwo)

                  vphi = 0.0E0_dp

                  do j = 3, counta + count - 1
                     ku = list(j)

                     angle=cone_angle(bendang(ku,iuprev),phi(j),thetatwo,phitwo)

                     do ib = 1, inben(imolty,iu)
                        iu2back = ijben3(imolty,iu,ib)

                        if (iu2back.eq.ku) then
                           type = itben(imolty,iu,ib)

                           vphi = vphi + brbenk(type) * (angle - brben(type))**2

                        end if
                     end do
                  end do

                  ang_trial(1) = phitwo
                  bfactor(1) = exp( -beta * vphi )
                  bsum_try = bsum_try + bfactor(1)

                  start = 2
               else
                  start = 1
                  thetatwo = bendang(iu,iuprev)
               end if

               do ibend = start, nchben_b
                  phitwo = random(-1) * twopi
                  vphi = 0.0E0_dp
                  do j = 3, counta + count - 1
                     ku = list(j)

                     angle=cone_angle(bendang(ku,iuprev),phi(j),thetatwo,phitwo)

                     do ib = 1, inben(imolty,iu)
                        iu2back = ijben3(imolty,iu,ib)
                        if (iu2back.eq.ku) then
                           type = itben(imolty,iu,ib)
                           vphi = vphi + brbenk(type)  * (angle - brben(type))**2
                        end if
                     end do
                  end do
! store the boltzmann factors and phi
                  bfactor(ibend) = exp(-beta*vphi)
                  ang_trial(ibend) = phitwo
                  bsum_try = bsum_try + bfactor(ibend)
               end do

               if ( lnew.or.ip.ne.1 ) then
                  rbf = random(-1) * bsum_try
                  bs = 0.0E0_dp
                  do ibend = 1, nchben_b
                     bs = bs + bfactor(ibend)
                     if (rbf .lt. bs) then
                        phitwo = ang_trial(ibend)
                        vphi = log( bfactor(ibend)) / (-beta)
                        goto 15
                     end if
                  end do
               else
                  phitwo = ang_trial(1)
                  vphi = log( bfactor(1) )/ (-beta)
               end if

 15            continue

               wei_bend = wei_bend * (bsum_try/dble(nchben_b))

               vbbtr = vbbtr + vphi

               phi(counta+count) = phitwo
               list(counta+count) = iu

! determine vectors associated with this

               call cone(2,ux,uy,uz,thetatwo,phitwo)

               xvec(iufrom,iu) = ux * distij(iufrom,iu)
               yvec(iufrom,iu) = uy * distij(iufrom,iu)
               zvec(iufrom,iu) = uz * distij(iufrom,iu)

               xvec(iu,iufrom) = -ux
               yvec(iu,iufrom) = -uy
               zvec(iu,iufrom) = -uz

! now to calculate all torsions

               do it = 1, intor(imolty,iu)

                  jut2 = ijtor2(imolty,iu,it)
                  jut3 = ijtor3(imolty,iu,it)
                  jut4 = ijtor4(imolty,iu,it)

                  if (niplace(jut4).lt.niplace(iu)) then

                     if (.not. lexist(jut4)) then
                        write(io_output,*) 'jut4,jut3,jut2,iu', jut4,jut3,jut2,iu
                        call err_exit(__FILE__,__LINE__,'trouble jut4 in place',myid+1)
                     end if

                     vdha = vdha + vtorso(xvec(iu,jut2),yvec(iu,jut2),zvec(iu,jut2),xvec(jut2,jut3),yvec(jut2,jut3)&
                      ,zvec(jut2,jut3),xvec(jut3,jut4),yvec(jut3,jut4),zvec(jut3,jut4),ittor(imolty,iu,it))
                  end if
               end do
            end do

            do count = 1, pnum(iw)
               iu = iplace(iw,count)
               if (lnew) then
                  rxp(count,ip) = xvec(iufrom,iu) + rxnew(iufrom)
                  ryp(count,ip) = yvec(iufrom,iu) + rynew(iufrom)
                  rzp(count,ip) = zvec(iufrom,iu) + rznew(iufrom)
               else
                  rxp(count,ip) = xvec(iufrom,iu) + rxu(i,iufrom)
                  ryp(count,ip) = yvec(iufrom,iu) + ryu(i,iufrom)
                  rzp(count,ip) = zvec(iufrom,iu) + rzu(i,iufrom)
               end if
               glist(count) = iu
            end do
            vtorsion(ip) = vdha
            vbend(ip) = vbbtr
            bsum_tor(ip) = exp(-beta * vdha) * wei_bend
         end do

! now we calculate the nonbonded interactions
         call boltz(lnew,.false.,ovrlap,i,i,imolty,ibox,ichoi,iufrom,pnum(iw),glist,0._dp)

         if (ovrlap) then
            lterm = .true.
            return
         end if

         bsum = 0.0E0_dp

         do ip = 1, ichoi
            if (.not. lovr(ip)) then
               bsum = bsum + bfac(ip) * bsum_tor(ip)
            end if
         end do

! update new rosenbluth weight + vibrations
         wplace = wplace * bsum * wei_vib

         if (wplace .lt. softlog) then
            lterm = .true.
            return
         end if

         if (lnew) then
            rbf = bsum * random(-1)
            bs = 0.0E0_dp

            do ip = 1, ichoi
               if (.not. lovr(ip)) then
                  bs = bs + bfac(ip) * bsum_tor(ip)
                  if (rbf .lt. bs) then
                     iwalk = ip
                     goto 20
                  end if
               end if
            end do

            call err_exit(__FILE__,__LINE__,'BIG TIME SCREWUP IN PLACE',myid+1)

         end if
 20      continue

         if (lnew) then

            vnew(ivTot) = vnew(ivTot) + vbend(iwalk) + vtorsion(iwalk) + vtr(ivIntraLJ,iwalk)
            vnew(ivBending)    = vnew(ivBending)    + vbend(iwalk)
            vnew(ivTorsion)    = vnew(ivTorsion)    + vtorsion(iwalk)
            vnew(ivExt)   = vnew(ivExt)   + vtr(ivExt,iwalk)
            vnew(ivIntraLJ) = vnew(ivIntraLJ) + vtr(ivIntraLJ,iwalk)
            vnew(ivInterLJ) = vnew(ivInterLJ) + vtr(ivInterLJ,iwalk)
            vnew(ivElect) = vnew(ivElect) + vtr(ivElect,iwalk)
            vnew(ivEwald) = vnew(ivEwald) + vtr(ivEwald,iwalk)
         else
            vold(ivTot) = vold(ivTot) + vbend(1) + vtorsion(1) + vtr(ivIntraLJ,1)
            vold(ivBending)    = vold(ivBending)    + vbend(1)
            vold(ivTorsion)    = vold(ivTorsion)    + vtorsion(1)
            vold(ivExt)   = vold(ivExt)   + vtr(ivExt,1)
            vold(ivIntraLJ) = vold(ivIntraLJ) + vtr(ivIntraLJ,1)
            vold(ivInterLJ) = vold(ivInterLJ) + vtr(ivInterLJ,1)
            vold(ivElect) = vold(ivElect) + vtr(ivElect,1)

            vold(ivEwald) = vold(ivEwald) + vtr(ivEwald,1)
         end if

         do count = 1, pnum(iw)
            iu = iplace(iw,count)

            if (lnew) then
               rxnew(iu) = rxp(count,iwalk)
               rynew(iu) = ryp(count,iwalk)
               rznew(iu) = rzp(count,iwalk)
            end if

            lexist(iu) = .true.

            if (lnew) then
               xvec(iu,iufrom) = rxnew(iufrom) - rxnew(iu)
               yvec(iu,iufrom) = rynew(iufrom) - rynew(iu)
               zvec(iu,iufrom) = rznew(iufrom) - rznew(iu)
            else
               xvec(iu,iufrom) = rxu(i,iufrom) - rxu(i,iu)
               yvec(iu,iufrom) = ryu(i,iufrom) - ryu(i,iu)
               zvec(iu,iufrom) = rzu(i,iufrom) - rzu(i,iu)
            end if
            distij(iu,iufrom) = sqrt( xvec(iu,iufrom)**2 + yvec(iu,iufrom)**2 + zvec(iu,iufrom)**2 )

            xvec(iufrom,iu) = - xvec(iu,iufrom)
            yvec(iufrom,iu) = - yvec(iu,iufrom)
            zvec(iufrom,iu) = - zvec(iu,iufrom)
            distij(iufrom,iu) = distij(iu,iufrom)
         end do
      end do

#ifdef __DEBUG__
      write(io_output,*) 'END PLACE in ',myid
#endif

      return
  end subroutine place

!*************************************************************
! Self-Adapting Fixed-Endpoint Configurational-Bias
! Monte Carlo SAFE-CBMC
!
! Most work is in safeschedule, safecbmc.f, and close.f.
!*************************************************************

!> \brief Finshes the last two steps for Fixed Endpoint CBMC
!> \param iinit iinit = 1  initial setup for two beads to go \n
!> iinit = 2  calculates closing probabilities
!> iinit = 3  does final crankshaft move
!> \see safeschedule.f FOR MORE INFORMATION
!> \author Originally completed by Collin Wick on 1-1-2000
  subroutine safecbmc(iinit,lnew,i,iw,igrow,imolty,count,ux,uy,uz,vphi,vtor,wei_bv,lterm,movetype)
      logical::lnew,lshit,lterm,ldo,lreturn ! lshit is used for diagnistics
      integer::igrow,imolty,count,counta,j,ja,ivib,iufrom,iuprev,iinit,iu,ju,ku,i,iv,juvib,jtvib,type,iu2,ib,iw,ntogrow,itor,ip&
       ,ichoi,ichtor,countb,bin,max,nu,iu1,dir,diracc,start,nchben_a,nchben_b,ibend,iopen,last,iclose,nchvib

      integer::it,jut2,jut3,jut4,movetype,lu,k,opencount

      real::vdha,phicrank,bf_tor,vtorsion,vbend,rbf,ran_tor,bs,ang_bend,bfactor,bsum_bend,wei_bv,bsum_try,third,vibtr

! to conserve memory, max is the maximum number of endpoints
! possible in place of numax
      parameter(max=10)

      real::x,y,z,equil,kforce,length,vvib,ux,uy,uz,hdist,lengtha,lengthb,vtor,vphi,thetac,angle,equila,kforcea,ovphi,alpha&
       ,phidisp,dum,rxt,ryt,rzt,phiacc,rxa,rya,rza,angles,bangles,r,mincb,delcb,vvibration,ovvib

      dimension alpha(max,numax),equila(max),rxa(max,max),rya(max,max),rza(max,max),rxt(max),ryt(max),rzt(max)&
       ,vbend(2*nchtor_max),phicrank(2*nchtor_max,max),phiacc(max),vtorsion(2*nchtor_max),bf_tor(2*nchtor_max)&
       ,dir(2*nchtor_max,max),diracc(max),kforcea(max),bfactor(nchbn_max),ang_bend(nchbn_max),angles(3),bangles(max,3)&
       ,iopen(2),r(nchbn_max),vvibration(2*nchtor_max)

!     ---------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'START SAFECMBC in ',myid
      write(io_output,*) 'iinit',iinit,'iw',iw,'igrow',igrow,'count',count
#endif

      vphi = 0.0E0_dp
      ovphi = 0.0E0_dp
      wei_bv = 1.0E0_dp

      lreturn = .false.

 400  continue

      if (iinit.eq.1.or.(.not.lreturn.and.lcrank)) then

         third = 1.0E0_dp / 3.0E0_dp
         nchvib = nchbna(imolty)
         ntogrow = grownum(iw)

         if (lcrank) then
            ntogrow = 1
         else
            ntogrow = grownum(iw)
         end if

! lets first determine our bond distances ***

! find vibrations for iu - ibef
         do 155 count = 1, ntogrow

            ja = 0

            if (lcrank) then
               wbefnum = 1
               iwbef(1) = growfrom(iw)
               iu = growfrom(iw)
            else
               iu = growlist(iw,count)
            end if
            do j = 1, wbefnum
               if (iu.eq.iwbef(j)) then
! ja defines iwbef in counter
                  ja = j
               end if
            end do

! we do not close with this bead
            if (ja.eq.0) goto 155

            do iv = 1, invib(imolty,iu)
               juvib = ijvib(imolty,iu,iv)
               do counta = 1, befnum(ja)
                  ju = ibef(ja,counta)
                  if (juvib.eq.ju) then
                     jtvib = itvib(imolty,iu,iv)
                     equil = brvib(jtvib)
                     kforce = brvibk(jtvib)
                     if (kforce.gt.0.1E0_dp) then
! we will use flexible bond lengths
                        bsum_try = 0.0E0_dp
                        mincb = brvibmin(jtvib)**3
                        delcb = brvibmax(jtvib)**3 - mincb

                        if (.not.lnew) then
                           x = rxu(i,iu) - rxu(i,ju)
                           y = ryu(i,iu) - ryu(i,ju)
                           z = rzu(i,iu) - rzu(i,ju)
                           r(1) = sqrt(x**2 + y**2 + z**2)
                           vvib = kforce * ( r(1) - equil )**2
                           bfactor(1) = exp(-beta*vvib)
                           bsum_try = bsum_try + bfactor(1)
                           start = 2
                        else
                           start = 1
                        end if

                        do ivib = start, nchvib

                           r(ivib) = (mincb + random(-1)*delcb)**third
                           vvib = kforce * ( r(ivib) - equil )**2
                           bfactor(ivib) = exp(-beta*vvib)
                           bsum_try = bsum_try + bfactor(ivib)
                        end do

                        if (lnew) then
! select one of the trial sites via vias
                           rbf = random(-1)*bsum_try
                           bs = 0.0E0_dp
                           do ivib = 1, nchvib
                              bs = bs + bfactor(ivib)
                              if (rbf .lt. bs ) then
                                 flength(iu,ju) = r(ivib)
                                 vvib = log(bfactor(ivib))/(-beta)
                                 goto 4
                              end if
                           end do
 4                         continue
                        else
! select old conformation
                           flength(iu,ju) = r(1)
                           vvib = log(bfactor(1))/(-beta)
                        end if

! add up vibrational energy
                        vphi = vphi + vvib

! propogate rosenbluth weight
                        wei_bv = wei_bv * bsum_try/dble(nchvib)

                     else
                        if (lnew) then
! compute new bond length
                           call bondlength(jtvib,equil,kforce,beta, length,vvib)
                           flength(iu,ju) = length

                        else
! compute old bond length
                           x = rxu(i,ju) - rxu(i,iu)
                           y = ryu(i,ju) - ryu(i,iu)
                           z = rzu(i,ju) - rzu(i,iu)
                           flength(iu,ju) = sqrt(x**2+y**2+z**2)
                        end if
                     end if
                  end if
               end do
            end do

! find vibrations for ibef - iend
            do counta = 1, befnum(ja)
               ju = ibef(ja,counta)
               do iv = 1, invib(imolty,ju)
                  juvib = ijvib(imolty,ju,iv)
                  do j = 1, fcount(ju)
                     ku = fclose(ju,j)
                     if (juvib.eq.ku) then
                        jtvib = itvib(imolty,ju,iv)
                        equil = brvib(jtvib)
                        kforce = brvibk(jtvib)
                        if (kforce.gt.0.1E0_dp) then
! we have flexible bond lengths

                           bsum_try = 0.0E0_dp
                           mincb = brvibmin(jtvib)**3
                           delcb = brvibmax(jtvib)**3 - mincb

                           if (j.gt.1) then
                              vequil(ju,ku) = equil
                              vkforce(ju,ku) = kforce
                              goto 112
                           end if

                           if (.not.lnew) then
                              x = rxu(i,ku) - rxu(i,ju)
                              y = ryu(i,ku) - ryu(i,ju)
                              z = rzu(i,ku) - rzu(i,ju)
                              r(1) = sqrt(x**2 + y**2 + z**2)
                              vvib = kforce * ( r(1) - equil )**2
                              bfactor(1) = exp(-beta*vvib)
                              bsum_try = bsum_try + bfactor(1)
                              start = 2
                           else
                              start = 1
                           end if

                           do ivib = start, nchvib

                              r(ivib) = (mincb + random(-1)*delcb)**third
                              vvib = kforce * ( r(ivib) - equil )**2
                              bfactor(ivib) = exp(-beta*vvib)
                              bsum_try = bsum_try + bfactor(ivib)
                           end do

                           if (lnew) then
! select one of the trial sites via vias
                              rbf = random(-1)*bsum_try
                              bs = 0.0E0_dp
                              do ivib = 1, nchvib
                                 bs = bs + bfactor(ivib)
                                 if (rbf .lt. bs ) then
                                    flength(ju,ku) = r(ivib)
                                    vvib = log(bfactor(ivib)) /(-beta)
                                    goto 6
                                 end if
                              end do
 6                            continue
                           else
! select old conformation
                              flength(ju,ku) = r(1)
                              vvib = log(bfactor(1))/(-beta)

                           end if

! add up vibrational energy
                           vphi = vphi + vvib

! propogate rosenbluth weight
                           wei_bv = wei_bv * bsum_try/dble(nchvib)

 112                       continue

                        else

                           if (lnew) then
! compute new bond length
                              call bondlength(jtvib,equil,kforce,beta ,length,vvib)
                              flength(ju,ku) = length
                           else
! compute old bond length
                              x = rxu(i,ku) - rxu(i,ju)
                              y = ryu(i,ku) - ryu(i,ju)
                              z = rzu(i,ku) - rzu(i,ju)

                              flength(ju,ku) = sqrt(x**2 + y**2 + z**2)

                           end if
                        end if
                     end if
                  end do
               end do
            end do

! determine angles for iwbef-ibef-iend

            if (lcrank) then
               do ib = 1, inben(imolty,iu)
                  iu2 = ijben3(imolty,iu,ib)
                  type = itben(imolty,iu,ib)
                  do counta = 1, grownum(iw)
                     ju = growlist(iw,counta)
                     if (fcount(ju).ne.0) then
                        do j = 1, fcount(ju)
                           ku = fclose(ju,j)
                           if (ku.eq.iu2) then
                              equilb(iu,ku) = brben(type)
                              kforceb(iu,ku) = brbenk(type)
                           end if
                        end do
                     end if
                  end do
               end do
            else
               do ib = 1, inben(imolty,iu)
                  iu2 = ijben3(imolty,iu,ib)
                  type = itben(imolty,iu,ib)
                  do j = 1, fcount(iu)
                     ju = fclose(iu,j)
                     if (ju.eq.iu2) then
                        equilb(iu,ju) = brben(type)
                        kforceb(iu,ju) = brbenk(type)
                     end if
                  end do
               end do
            end if
 155     continue

         if (lcrank) then
! we need to calculate the new angle
            iufrom = growfrom(iw)
            do count = 1, grownum(iw)
               iu = growlist(iw,count)
               if (fcount(iu).ne.0) then
                  j = 1
                  ju = fclose(iu,j)
                  lengtha = flength(iufrom,iu)
                  lengthb = flength(iu,ju)

                  x = rxu(i,ju) - rxnew(iufrom)
                  y = ryu(i,ju) - rynew(iufrom)
                  z = rzu(i,ju) - rznew(iufrom)

                  hdist = sqrt( x**2 + y**2 + z**2 )
! use law of cosines to calculate bond angle
                  thetac = (lengtha**2 + lengthb**2 - hdist**2) / (2.0E0_dp * lengtha * lengthb)

! check to make sure this will give a number

                  if (abs(thetac).gt.1.0E0_dp) then
                     vtor = 0
                     vphi = 0
                     lterm = .true.
                     return
                  end if
                  angle = acos(thetac)

                  ovphi = ovphi + kforceb(iufrom,ju) * (angle - equilb(iufrom,ju))**2

! write(io_output,*) iufrom,iu,ju,kforceb(iufrom,ju)*(angle
!     &                 - equilb(iufrom,ju))**2

                  wei_bv = wei_bv * exp( - beta * ovphi )
               end if
            end do
            vibtr= vphi

            lreturn = .true.
            goto 400
         end if

!     ********************************************************************
      else if (iinit.eq.2.or.iinit.eq.4) then
! lets determine closing energy for this bead alone
         lshit = .false.

! if (iinit.eq.4) then
! lshit = .true.
! end if

         vtor = 1.0E0_dp

! this case iw and count is sent in from rosenbluth

         iu = growlist(iw,count)

! determine iwbef count
         do j = 1, wbefnum
            if (iu.eq.iwbef(j)) then
               ja = j
               goto 100
            end if
         end do
 100     continue

         do j = 1, fcount(iu)
            ku = fclose(iu,j)

            if (movetype.eq.2.and.lnew) then
               x = rxnew(ku) - ux
               y = rynew(ku) - uy
               z = rznew(ku) - uz
            else
               x = rxu(i,ku) - ux
               y = ryu(i,ku) - uy
               z = rzu(i,ku) - uz
            end if

            hdist = sqrt( x**2 + y**2 + z**2 )

            if (j.gt.1) then
! we will use a phoney probability since we don't know our bond
! distances yet
               bin = anint( hdist * 10.0E0_dp )
               vtor = vtor * probf(iu,ku,bin)
            else
! we can calculate an angle here
               do counta = 1, befnum(ja)
                  ju = ibef(ja,counta)

                  do countb = 1, fcount(ju)
                     nu = fclose(ju,countb)
                     if (nu.eq.ku) goto 150
                  end do
                  goto 175

 150              continue

                  lengtha = flength(iu,ju)
                  lengthb = flength(ju,ku)

! use law of cosines to calculate bond angle

                  thetac = (lengtha**2 + lengthb**2 - hdist**2) / (2.0E0_dp * lengtha * lengthb)

! check to make sure this will give a number

                  if (abs(thetac).gt.1.0E0_dp) then
                     vtor = 0
                     vphi = 0
                     return
                  end if
                  angle = acos( thetac )

                  vphi = vphi + kforceb(iu,ku) * ( angle - equilb(iu,ku) )**2

 175              continue

               end do
            end if

! determine torsion interaction with growpast if it exists

            if (pastnum(ku).ne.0) then

               do counta = 1, pastnum(ku)
                  nu = ipast(ku,counta)
                  if (.not.lplace(imolty,nu)) then
                     if (movetype.eq.2.and.lnew) then
                        x = rxnew(nu) - ux
                        y = rynew(nu) - uy
                        z = rznew(nu) - uz
                     else
                        x = rxu(i,nu) - ux
                        y = ryu(i,nu) - uy
                        z = rzu(i,nu) - uz
                     end if
                     hdist = sqrt( x**2 + y**2 + z**2 )
                     bin = anint( hdist * 10.0E0_dp )
                     vtor = vtor * probf(iu,nu,bin)

                     if (nextnum(nu).ne.0) then
                        do k = 1, nextnum(nu)
                           lu = inext(nu,k)
                           if (.not.lplace(imolty,lu)) then
                              if (movetype.eq.2.and.lnew) then
                                 x = rxnew(lu) - ux
                                 y = rynew(lu) - uy
                                 z = rznew(lu) - uz
                              else
                                 x = rxu(i,lu) - ux
                                 y = ryu(i,lu) - uy
                                 z = rzu(i,lu) - uz
                              end if
                              hdist = sqrt( x**2 + y**2 + z**2 )
                              bin = anint( hdist * 10.0E0_dp )
                              vtor = vtor * probf(iu,lu,bin)
                           end if
                        end do
                     end if
                  end if
               end do
            end if
         end do

! *********************************************************************
      else

! CRANKSHAFT MOVE
         ntogrow = grownum(iw)
         iufrom = growfrom(iw)
         iuprev = growprev(iw)
         iopen(1) = 0
         iopen(2) = 0
         iclose = 0
         ovvib = 0
! determine bond distances for iuprev - iufrom

         opencount = 0

         do count = 1, ntogrow
            iu = growlist(iw,count)
            if (fcount(iu).ne.0) then

               do j = 1, fcount(iu)
                  ju = fclose(iu,j)
                  if (pastnum(ju).ne.0) then
                     do counta = 1, pastnum(ju)
                        ku = ipast(ju,counta)

! determine angles for iu - iend - ipast if they exist

                        do ib = 1, inben(imolty,iu)
                           iu2 = ijben3(imolty,iu,ib)
                           type = itben(imolty,iu,ib)
                           if (iu2.eq.ku) then
                              equilb(iu,ku) = brben(type)
                              kforceb(iu,ku) = brbenk(type)
                              goto 125
                           end if
                        end do
                        write(io_output,*) 'iu,ju,ku',iu,ju,ku
                        call err_exit(__FILE__,__LINE__,'no bond angle for these',myid+1)
 125                    continue

                        if (pastnum(ju).gt.1) then
                           do countb = counta+1, pastnum(ju)
                              lu = ipast(ju,countb)

                              do ib = 1, inben(imolty,ku)
                                 iu2 = ijben3(imolty,ku,ib)
                                 type = itben(imolty,ku,ib)
                                 if (iu2.eq.lu) then
                                    equilb(ku,lu) = brben(type)
                                    kforceb(ku,lu) = brbenk(type)
                                 end if
                              end do
                           end do
                        end if

                     end do
                  end if

                  if (j.gt.1) goto 126

! calculate distances from iufrom to iend
! count is iu count, and ju is iend bead
                  if (lnew) then
                     if (movetype.eq.2) then
                        xvec(iufrom,ju) = rxnew(ju) - rxnew(iufrom)
                        yvec(iufrom,ju) = rynew(ju) - rynew(iufrom)
                        zvec(iufrom,ju) = rznew(ju) - rznew(iufrom)
                     else
                        xvec(iufrom,ju) = rxu(i,ju) - rxnew(iufrom)
                        yvec(iufrom,ju) = ryu(i,ju) - rynew(iufrom)
                        zvec(iufrom,ju) = rzu(i,ju) - rznew(iufrom)
                     end if
                  else
                     xvec(iufrom,ju) = rxu(i,ju) - rxu(i,iufrom)
                     yvec(iufrom,ju) = ryu(i,ju) - ryu(i,iufrom)
                     zvec(iufrom,ju) = rzu(i,ju) - rzu(i,iufrom)
                  end if
                  hdist = sqrt( xvec(iufrom,ju)**2 + yvec(iufrom,ju)**2 + zvec(iufrom,ju)**2 )

! normalize these distances to one for cone
                  xvec(iufrom,ju) = xvec(iufrom,ju) / hdist
                  yvec(iufrom,ju) = yvec(iufrom,ju) / hdist
                  zvec(iufrom,ju) = zvec(iufrom,ju) / hdist

! calculate alpha
! count is iu count, and ju is iend bead
                  lengtha = flength(iufrom,iu)
                  lengthb = flength(iu,ju)

                  thetac = (lengtha**2 + hdist**2 - lengthb**2) / (2.0E0_dp * lengtha * hdist)

                  if (abs(thetac).gt.1.0E0_dp) then
                     vphi = 0
                     vtor = 0
                     lterm = .true.
                     return
                  end if
                  alpha(count,ju) = acos(-1.0E0_dp) - acos(thetac)
 126              continue
               end do
            end if

! determine angle for    - iu - iufrom -
            do ib = 1, inben(imolty,iu)
               iu1 = ijben2(imolty,iu,ib)
               type = itben(imolty,iu,ib)

               if (iu1.eq.iufrom) then
                  iu2 = ijben3(imolty,iu,ib)
                  if (iu2.eq.iuprev) then
                     equila(count) = brben(type)
                     kforcea(count) = brbenk(type)
                  else
                     equilb(iu,iu2) = brben(type)
                     kforceb(iu,iu2) = brbenk(type)
                  end if
               end if
            end do

            if (fcount(iu).gt.1) then
               do j = 1, fcount(iu) - 1
                  ju = fclose(iu,j)
                  do counta = j+1, fcount(iu)
                     ku = fclose(iu,counta)
                     do ib = 1, inben(imolty,ku)
                        iu2 = ijben3(imolty,ku,ib)
                        type = itben(imolty,ku,ib)

                        if (iu2.eq.ju) then
                           equilb(ju,ku) = brben(type)
                           kforceb(ju,ku) = brbenk(type)
                        end if
                     end do
                  end do
               end do
            end if

! ---------------------------------------------------------------
! begin part that is only for closures with an open bead

            if (fcount(iu).eq.0) then

               opencount = opencount + 1

! we first have to find the bond length here
               do iv = 1, invib(imolty,iufrom)
                  juvib = ijvib(imolty,iufrom,iv)
                  if (juvib.eq.iu)  then
                     jtvib = itvib(imolty,iufrom,iv)
                     equil = brvib(jtvib)
                     kforce = brvibk(jtvib)
                     if (kforce.gt.0.1E0_dp) then
! we will use flexible bond lengths
                        bsum_try = 0.0E0_dp
                        mincb = brvibmin(jtvib)**3
                        delcb = brvibmax(jtvib)**3 - mincb

                        if (.not.lnew) then
                           x = rxu(i,iufrom) - rxu(i,iu)
                           y = ryu(i,iufrom) - ryu(i,iu)
                           z = rzu(i,iufrom) - rzu(i,iu)
                           r(1) = sqrt(x**2 + y**2 + z**2)
                           vvib = kforce * ( r(1) - equil )**2
                           bfactor(1) = exp(-beta*vvib)
                           bsum_try = bsum_try + bfactor(1)
                           start = 2
                        else
                           start = 1
                        end if

                        do ivib = start, nchvib

                           r(ivib) = (mincb  + random(-1)*delcb)**third
                           vvib = kforce *  ( r(ivib) - equil )**2
                           bfactor(ivib) = exp(-beta*vvib)
                           bsum_try = bsum_try + bfactor(ivib)
                        end do

                        if (lnew) then
! select one of the trial sites via vias
                           rbf = random(-1)*bsum_try
                           bs = 0.0E0_dp
                           do ivib = 1, nchvib
                              bs = bs + bfactor(ivib)
                              if (rbf .lt. bs ) then
                                 flength(iufrom,iu) = r(ivib)
                                 vvib = log(bfactor(ivib))/(-beta)
                                 goto 61
                              end if
                           end do
 61                        continue
                        else
! select old conformation
                           flength(iufrom,iu) = r(1)
                           vvib = log(bfactor(1))/(-beta)
                        end if

! add up vibrational energy
                        ovvib = ovvib + vvib

! propogate rosenbluth weight
                        wei_bv = wei_bv * bsum_try/dble(nchvib)

!     *************************************
                     else

                        if (lnew) then
! compute new bond length
                           call bondlength(jtvib,equil,kforce,beta, length,vvib)
                           flength(iufrom,iu) = length
                        else
! compute old bond length
                           x = rxu(i,iufrom) - rxu(i,iu)
                           y = ryu(i,iufrom) - ryu(i,iu)
                           z = rzu(i,iufrom) - rzu(i,iu)
                           flength(iufrom,iu) = sqrt(x**2+y**2+z**2)

                        end if
                     end if
                  end if
               end do

! find bond angles

               bsum_bend = 0
               if (.not. lnew) then

                  lengtha = flength(iufrom,iu)
                  lengthb = distij(iufrom,iuprev)

                  thetac = -( (rxu(i,iu) - rxu(i,iufrom)) * xvec(iuprev,iufrom)  + (ryu(i,iu) - ryu(i,iufrom))&
                   * yvec(iuprev,iufrom) + (rzu(i,iu) - rzu(i,iufrom)) * zvec(iuprev,iufrom))  / (lengtha*lengthb)
                  angle = acos(thetac)
                  vphi =  kforcea(count) * (angle-equila(count))**2

! write(io_output,*) 'b',iu,iufrom,iuprev,vphi

                  ang_bend(1) = angle
                  bfactor(1) = exp( -beta*vphi)
                  bsum_bend = bsum_bend + bfactor(1)
                  start = 2
               else
                  start = 1
               end if

               nchben_a = nchbna(imolty)

               do ibend = start, nchben_a
! choose angle uniformly on sin(angle)
                  thetac = 2.0E0_dp*random(-1) - 1.0E0_dp
                  angle = acos(thetac)
                  ang_bend(ibend) = angle

! find bend energy
                  vphi = kforcea(count) * (angle-equila(count))**2
                  bfactor(ibend) = exp(-beta*vphi)
                  bsum_bend = bsum_bend + bfactor(ibend)
               end do

               if (lnew) then
! select one of the trial sites at random
                  rbf = random(-1)*bsum_bend
                  bs = 0.0E0_dp
                  do ibend = 1, nchben_a
                     bs = bs + bfactor(ibend)
                     if (rbf.lt.bs) then
                        bangles(count,1) = ang_bend(ibend)
                        ovphi = ovphi + log(bfactor(ibend))/(-beta)

! write(io_output,*) 'c',iu,iufrom,iuprev,log(bfactor(ibend))/(-beta)

                        goto 10
                     end if
                  end do
 10               continue
               else
                  bangles(count,1) = ang_bend(1)
                  ovphi = ovphi + log(bfactor(1))/(-beta)
               end if

               wei_bv = wei_bv * bsum_bend/dble(nchben_a)

! now we have to determine the angle with a crankshaft bead
! for the opencount = 1, or the other free bead for opencount > 1

               bsum_bend = 0

               do counta = 1, ntogrow

                  if (counta.ne.count) then
                     ju = growlist(iw,counta)

! check to see if the conditions just stated are true
                     if (opencount.gt.1) then
! we only want the angle with a free bead
                        if (fcount(ju).ne.0) then
                           goto 25
                        end if
                     else
! we want the angle with the closing bead
                        if (fcount(ju).eq.0) then
                           goto 25
                        else if (iclose.ne.0) then
                           goto 25
                        else
                           iclose = counta
                        end if
                     end if

                     if (.not. lnew) then

                        lengthb = flength(iufrom,ju)

                        thetac = ((rxu(i,ju)-rxu(i,iufrom)) * (rxu(i,iu)-rxu(i,iufrom)) + (ryu(i,ju)-ryu(i,iufrom))&
                         * (ryu(i,iu)-ryu(i,iufrom)) + (rzu(i,ju)-rzu(i,iufrom)) * (rzu(i,iu)-rzu(i,iufrom)))/(lengtha*lengthb)
                        angle = acos(thetac)

                        vphi = kforceb(iu,ju) * (angle-equilb(iu,ju))**2

! write(io_output,*) 'd',iu,iufrom,ju,vphi

                        ang_bend(1) = angle
                        bfactor(1) = exp( -beta*vphi )
                        bsum_bend = bsum_bend + bfactor(1)

                        start = 2
                     else
                        start = 1
                     end if

                     nchben_b = nchbnb(imolty)
                     do ibend = start, nchben_b
! choose angle uniformly on sin(angle)
                        thetac = 2.0E0_dp*random(-1) - 1.0E0_dp
                        angle = acos(thetac)
                        ang_bend(ibend) = angle

! find bend energy
                        vphi = kforceb(iu,ju) * (angle-equilb(iu,ju))**2

                        bfactor(ibend) = exp(-beta*vphi)
                        bsum_bend = bsum_bend + bfactor(ibend)
                     end do

                     if (lnew) then
! select one of the trial sites at random
                        rbf = random(-1)*bsum_bend
                        bs = 0.0E0_dp
                        do ibend = 1, nchben_b
                           bs = bs + bfactor(ibend)
                           if (rbf.lt.bs) then
                              bangles(count,2) = ang_bend(ibend)
                              ovphi = ovphi  + log(bfactor(ibend))/(-beta)

! write(io_output,*) 'd',iu,iufrom,ju,log(bfactor(ibend))/(-beta)

                              goto 20
                           end if
                        end do
 20                     continue
                     else
                        bangles(count,2) = ang_bend(1)
                        ovphi = ovphi + log(bfactor(1))/(-beta)

                     end if

                     wei_bv = wei_bv * bsum_bend/dble(nchben_b)
                   end if
 25               continue
               end do
            end if
         end do

! end part for the case with an open closing bead

! -----------------------------------------------------------------
! loop over all choices
         ichoi = nchoi(imolty)
! double ichtor to give extra help for this move
         ichtor = nchtor(imolty) * 2

         do ip = 1, ichoi

            bsum_tor(ip) = 0
            do itor = 1, ichtor

               vvib = 0
               countb = 0
               lshit = .false.

               if (lnew.and.ip.eq.19 .and.itor.eq.122) then
                  lshit = .true.
               end if

               if (.not.lnew .and.ip.eq.1.and.itor.eq.1) then
                  lshit = .true.
               end if

               vdha = 0
               vphi = 0

               ldo = .false.
               start = 1
               last = ntogrow
 30            continue
               do count = start, last

                  dir(itor,count) = 0

                  iu = growlist(iw,count)

                  if (fcount(iu).gt.0) then
! set up cone

                     ju = fclose(iu,1)
                     call cone(1,xvec(iufrom,ju),yvec(iufrom,ju),zvec(iufrom,ju),dum,dum)

! determine phidisp
                     if (.not.lnew.and.ip.eq.1.and.itor.eq.1) then
! give old unit vector for connection
                        xx(count) = rxu(i,iu) - rxu(i,iufrom)
                        yy(count) = ryu(i,iu) - ryu(i,iufrom)
                        zz(count) = rzu(i,iu) - rzu(i,iufrom)

                        ux = rxu(i,iu)
                        uy = ryu(i,iu)
                        uz = rzu(i,iu)

                     else
                        phidisp = twopi*random(-1)

                        call cone(2,x,y,z,alpha(count,ju),phidisp)

                        phicrank(itor,count) = phidisp

                        xx(count) = x * flength(iufrom,iu)
                        yy(count) = y * flength(iufrom,iu)
                        zz(count) = z * flength(iufrom,iu)

                        if (lnew) then
                           ux = xx(count) + rxnew(iufrom)
                           uy = yy(count) + rynew(iufrom)
                           uz = zz(count) + rznew(iufrom)
                        else
                           ux = xx(count) + rxu(i,iufrom)
                           uy = yy(count) + ryu(i,iufrom)
                           uz = zz(count) + rzu(i,iufrom)
                        end if
                     end if

                     if (fcount(iu).gt.1) then
! calculate distances to the other endpoints

                        do j = 2, fcount(iu)

                           ju = fclose(iu,j)

                           equil = vequil(iu,ju)
                           kforce = vkforce(iu,ju)

                           if (movetype.eq.2.and.lnew) then
                              x = rxnew(ju) - ux
                              y = rynew(ju) - uy
                              z = rznew(ju) - uz
                           else
                              x = rxu(i,ju) - ux
                              y = ryu(i,ju) - uy
                              z = rzu(i,ju) - uz
                           end if

                           xvec(iu,ju) = x
                           yvec(iu,ju) = y
                           zvec(iu,ju) = z

                           xvec(ju,iu) = -x
                           yvec(ju,iu) = -y
                           zvec(ju,iu) = -z

                           length = sqrt(x**2 + y**2 + z**2)

                           flength(iu,ju) = length

                           vvib = vvib + kforce * (length - equil)**2

! we need to calculate the new angle

                           lengtha = flength(iufrom,iu)
                           lengthb = length

                           if (lnew) then
                              if (movetype.eq.2) then
                                 x = rxnew(ju) - rxnew(iufrom)
                                 y = rynew(ju) - rynew(iufrom)
                                 z = rznew(ju) - rznew(iufrom)
                              else
                                 x = rxu(i,ju) - rxnew(iufrom)
                                 y = ryu(i,ju) - rynew(iufrom)
                                 z = rzu(i,ju) - rznew(iufrom)
                              end if
                           else
                              x = rxu(i,ju) - rxu(i,iufrom)
                              y = ryu(i,ju) - ryu(i,iufrom)
                              z = rzu(i,ju) - rzu(i,iufrom)
                           end if

                           hdist = sqrt( x**2 + y**2 + z**2 )
! use law of cosines to calculate bond angle
                           thetac = (lengtha**2 + lengthb**2  - hdist**2) / (2.0E0_dp * lengtha * lengthb)

! check to make sure this will give a number

                           if (abs(thetac).gt.1.0E0_dp) then
                              vtor = 0
                              vphi = 0
                              vvib = 0
                              bf_tor(itor) = 0
                              goto 190
                           end if
                           angle = acos(thetac)

                           vphi = vphi + kforceb(iufrom,ju) * (angle - equilb(iufrom,ju))**2

! if (lshit) then
! write(io_output,*) iufrom,iu,ju,kforceb(iufrom,ju)*(angle-equilb(iufrom,ju))**2
! end if

                        end do
                     end if
                  else if (fcount(iu).eq.0) then
                     ! we use the angle with the closed bond
                     if (.not.ldo) then
                        ldo = .true.
                        countb = countb + 1
                        iopen(countb) = count
                        goto 40
                     else if (countb.eq.1.and.iopen(1).ne.count) then
                        countb = countb + 1
                        iopen(countb) = count
                        goto 40
                     else if (countb.eq.2) then
                        countb = countb + 1
                     else
                        ldo = .false.
                        countb = countb + 1
                     end if

                     if (.not.lnew.and.ip.eq.1.and.itor.eq.1) then
! give old unit vector for connection
                        xx(count) = rxu(i,iu) - rxu(i,iufrom)
                        yy(count) = ryu(i,iu) - ryu(i,iufrom)
                        zz(count) = rzu(i,iu) - rzu(i,iufrom)
                     else
                        do counta = 1, ntogrow
                           if (counta.ne.count) then

                              if (opencount.gt.1) then
                                 if (countb.eq.3) then
! we only want the angle with the closing bead
                                    if (iopen(2).eq.counta) then
                                       goto 35
                                    end if
                                 else if (countb.eq.4) then
! we want the angle with the open bead
                                    if (iopen(1).ne.counta) then
                                       goto 35
                                    end if
                                 end if
                              else if(iclose.ne.counta) then
                                 goto 35
                              end if
                              ju = growlist(iw,counta)

! determine position by angle with this bond

                              lengtha = distij(iufrom,iuprev)
                              lengthb = flength(iufrom,ju)

                              rxt(1) = xvec(iufrom,iuprev)/lengtha
                              ryt(1) = yvec(iufrom,iuprev)/lengtha
                              rzt(1) = zvec(iufrom,iuprev)/lengtha

                              rxt(2) = xx(counta)/lengthb
                              ryt(2) = yy(counta)/lengthb
                              rzt(2) = zz(counta)/lengthb

                              angles(1) = bangles(count,1)
                              angles(2) = bangles(count,2)

                              call close(3,rxt,ryt,rzt,dum ,angles,lterm)

                              if (lterm) then
                                 lterm = .false.
                                 vphi = 0
                                 vdha = 0
                                 vvib = 0
                                 bf_tor(itor) = 0
                                 goto 190
                              end if

! since there are two possibilities, choose one at random
                              lengtha = flength(iufrom,iu)
                              if (random(-1).lt.0.5E0_dp) then
                                 xx(count) = rxt(3) * lengtha
                                 yy(count) = ryt(3) * lengtha
                                 zz(count) = rzt(3) * lengtha
                                 dir(itor,count) = 1
                              else
                                 xx(count) = rxt(4) * lengtha
                                 yy(count) = ryt(4) * lengtha
                                 zz(count) = rzt(4) * lengtha
                                 dir(itor,count) = 2
                              end if

!     ------------------------------------------------------------------
 35                           continue
                           end if
                        end do
                     end if

                  else
! we have to close differently

                     goto 132

!     **************************************************************

                     if (ip.eq.1.and.itor.eq.1) then
                        if (lnew) then
                           rxt(1) = rxnew(iufrom)
                           ryt(1) = rynew(iufrom)
                           rzt(1) = rznew(iufrom)
                        else
                           rxt(1) = rxu(i,iufrom)
                           ryt(1) = ryu(i,iufrom)
                           rzt(1) = rzu(i,iufrom)
                        end if

                        do j = 1, fcount(iu)
                           ju = fclose(iu,j)
                           rxt(j+1) = rxu(i,ju)
                           ryt(j+1) = ryu(i,ju)
                           rzt(j+1) = rzu(i,ju)
                        end do

                        call close(1,rxt,ryt,rzt,flength(iufrom,iu) ,angles,lterm)

                        if (lterm) then
                           return
                        end if

                        if (lnew) then
                           x = rxt(1) - rxnew(iufrom)
                           y = ryt(1) - rynew(iufrom)
                           z = rzt(1) - rznew(iufrom)
                        else
                           x = rxt(1) - rxu(i,iufrom)
                           y = ryt(1) - ryu(i,iufrom)
                           z = rzt(1) - rzu(i,iufrom)
                        end if

                        length = sqrt(x**2+y**2+z**2)

                        if (lnew) then
                           x = rxt(2) - rxnew(iufrom)
                           y = ryt(2) - rynew(iufrom)
                           z = rzt(2) - rznew(iufrom)
                        else
                           x = rxt(2) - rxu(i,iufrom)
                           y = ryt(2) - ryu(i,iufrom)
                           z = rzt(2) - rzu(i,iufrom)
                        end if

                        length = sqrt(x**2+y**2+z**2)

                        rxa(count,1) = rxt(1)
                        rya(count,1) = ryt(1)
                        rza(count,1) = rzt(1)

                        rxa(count,2) = rxt(2)
                        rya(count,2) = ryt(2)
                        rza(count,2) = rzt(2)
                     end if

! determine position at random

                     j = int(2.0E0_dp * random(-1)) + 1

                     if (.not.lnew.and.ip.eq.1.and.itor.eq.1) then
! give old unit vector for connection
                        xx(count) = rxu(i,iu) - rxu(i,iufrom)
                        yy(count) = ryu(i,iu) - ryu(i,iufrom)
                        zz(count) = rzu(i,iu) - rzu(i,iufrom)
                     else
                        dir(itor,count) = j
                        if (lnew) then
                           xx(count) = rxa(count,j) - rxnew(iufrom)
                           yy(count) = rya(count,j) - rynew(iufrom)
                           zz(count) = rza(count,j) - rznew(iufrom)
                        else
                           xx(count) = rxa(count,j) - rxu(i,iufrom)
                           yy(count) = rya(count,j) - ryu(i,iufrom)
                           zz(count) = rza(count,j) - rzu(i,iufrom)
                        end if
                     end if
!     *****************************************************************

 132                 continue

                  end if

                  xvec(iufrom,iu) = xx(count)
                  xvec(iu,iufrom) = -xx(count)
                  yvec(iufrom,iu) = yy(count)
                  yvec(iu,iufrom) = -yy(count)
                  zvec(iufrom,iu) = zz(count)
                  zvec(iu,iufrom) = -zz(count)

! determine position of trial spot
                  if (lnew) then
                     rxt(count) = rxnew(iufrom) + xx(count)
                     ryt(count) = rynew(iufrom) + yy(count)
                     rzt(count) = rznew(iufrom) + zz(count)
                  else
                     rxt(count) = rxu(i,iufrom) + xx(count)
                     ryt(count) = ryu(i,iufrom) + yy(count)
                     rzt(count) = rzu(i,iufrom) + zz(count)
                  end if

                  lengtha = flength(iufrom,iu)

! now that we have the trial spot, determine the bending
! energies that are applicable

                  if (iopen(1).eq.count.or.iopen(2).eq.count) then
! we already calculated all necessary bending energies for these
                     goto 40
                  end if

! first determine bending energy for iuprev-iufrom-iu
                  if (iuprev.ne.0) then
                     length = distij(iufrom,iuprev)
                     thetac = -(xx(count)*xvec(iuprev,iufrom) + yy(count)*yvec(iuprev,iufrom) + zz(count)*zvec(iuprev,iufrom))&
                      /(lengtha*length)
                     angle = acos(thetac)

                     vphi = vphi + kforcea(count) * (angle-equila(count) )**2

! if (lshit) then
! write(io_output,*) iu,iufrom,iuprev,kforcea(count) *(angle-equila(count))**2
! end if
                  end if

                  do j = 1, fcount(iu)
! determine vectors from trial spot to iend

                     ju = fclose(iu,j)
                     if (movetype.eq.2.and.lnew) then
                        xvec(iu,ju) = rxnew(ju) - rxt(count)
                        yvec(iu,ju) = rynew(ju) - ryt(count)
                        zvec(iu,ju) = rznew(ju) - rzt(count)
                     else
                        xvec(iu,ju) = rxu(i,ju) - rxt(count)
                        yvec(iu,ju) = ryu(i,ju) - ryt(count)
                        zvec(iu,ju) = rzu(i,ju) - rzt(count)
                     end if
                     xvec(ju,iu) = - xvec(iu,ju)
                     yvec(ju,iu) = - yvec(iu,ju)
                     zvec(ju,iu) = - zvec(iu,ju)

                     lengthb = flength(iu,ju)

! now determine bending energy for iu-iend-ipast if it exists
                     if (pastnum(ju).ne.0) then
                        do counta = 1, pastnum(ju)
                           ku = ipast(ju,counta)
                           length = distij(ju,ku)
                           thetac = -(xvec(iu,ju)*xvec(ju,ku) + yvec(iu,ju)*yvec(ju,ku) + zvec(iu,ju)*zvec(ju,ku))/(length*lengthb)
                           angle = acos( thetac )

                           vphi = vphi + kforceb(iu,ku) * (angle -equilb(iu,ku))**2

! if (lshit) then
! write(io_output,*) iu,ju,ku,kforceb(iu,ku)*(angle-equilb(iu,ku))**2
! end if
                        end do
                     end if
                  end do
 40               continue
               end do
               if (ldo) then

                  if (countb.eq.3) then
                     start = iopen(2)
                     last = iopen(2)
                  else
                     start = iopen(1)
                     last = iopen(1)
                  end if
                  goto 30
               end if

! calculate four sets of torsion energies
! iu to iend to ipast to inext

! ------------------------------------------------------------------------
! first calculate torsion for iu-iufrom-
               do count = 1, ntogrow
                  iu = growlist(iw,count)
                  do it = 1, intor(imolty,iu)
                     jut2 = ijtor2(imolty,iu,it)

                     if (jut2.eq.iufrom) then
                        jut3 = ijtor3(imolty,iu,it)
                        jut4 = ijtor4(imolty,iu,it)

                        if (lpnow(jut4)) goto 41

! jut4 must already exist or we made a big mistake
                        if (.not. lexist(jut4)) then
                           if (jut4.gt.iring(imolty)) then
                              goto 41
                           end if
                           write(io_output,*) 'jut4,jut3,jut2,iu', jut4,jut3,jut2,iu
                           call err_exit(__FILE__,__LINE__,'trouble jut4 in crankshaft',myid+1)
                        end if

                        vdha = vdha + vtorso(xvec(iu,jut2),yvec(iu,jut2),zvec(iu,jut2),xvec(jut2,jut3),yvec(jut2,jut3)&
                         ,zvec(jut2,jut3),xvec(jut3,jut4),yvec(jut3,jut4),zvec(jut3,jut4),ittor(imolty,iu,it))

 41                     continue

                     end if
                  end do

! now lets calculate torsions for iend-iu-iufrom-iuprev

                  do j = 1, fcount(iu)
                     ju = fclose(iu,j)
                     do it = 1, intor(imolty,ju)
                        jut2 = ijtor2(imolty,ju,it)
                        jut3 = ijtor3(imolty,ju,it)
                        jut4 = ijtor4(imolty,ju,it)

                        if (jut2.eq.iu.and.jut4.eq.iuprev) then
! add torsion energy to vdha
                           vdha = vdha + vtorso(xvec(ju,jut2),yvec(ju,jut2),zvec(ju,jut2),xvec(jut2,jut3),yvec(jut2,jut3)&
                            ,zvec(jut2,jut3),xvec(jut3,jut4),yvec(jut3,jut4),zvec(jut3,jut4),ittor(imolty,ju,it))
                        end if
                     end do

! calculate torsions for ipast-iend-iu-

                  if (pastnum(ju).ne.0) then
                  do counta = 1, pastnum(ju)
                     ku = ipast(ju,counta)
                     do it = 1, intor(imolty,ku)
                        jut2 = ijtor2(imolty,ku,it)
                        jut3 = ijtor3(imolty,ku,it)

                        if (jut2.eq.ju.and.jut3.eq.iu) then
                           jut4 = ijtor4(imolty,ku,it)

                           if (lpnow(jut4)) goto 42

! jut4 must already exist or we made a big mistake
                           if (.not. lexist(jut4)) then
                              write(io_output,*) 'jut4,jut3,jut2,iu', jut4,jut3,jut2,iu
                              call err_exit(__FILE__,__LINE__,'trouble jut4 in crankshaft',myid+1)
                           end if

! add torsion energy to vdha
                           vdha = vdha + vtorso(xvec(ku,jut2),yvec(ku,jut2),zvec(ku,jut2),xvec(jut2,jut3),yvec(jut2,jut3)&
                            ,zvec(jut2,jut3),xvec(jut3,jut4),yvec(jut3,jut4),zvec(jut3,jut4),ittor(imolty,ku,it))
 42                        continue
                        end if
                     end do

! calculate torsions for inext-ipast-iend-iu

                  if (nextnum(ku).ne.0) then
                     do ja = 1, nextnum(ku)
                     nu = inext(ku,ja)

                     do it = 1, intor(imolty,nu)
                        jut2 = ijtor2(imolty,nu,it)
                        jut3 = ijtor3(imolty,nu,it)
                        jut4 = ijtor4(imolty,nu,it)

                        if (jut2.eq.ku.and.jut3.eq.ju.and.jut4.eq.iu) then
! add torsion energy to vdha
                           vdha = vdha + vtorso(xvec(nu,jut2),yvec(nu,jut2),zvec(nu,jut2),xvec(jut2,jut3),yvec(jut2,jut3)&
                            ,zvec(jut2,jut3),xvec(jut3,jut4),yvec(jut3,jut4),zvec(jut3,jut4),ittor(imolty,nu,it))
                        end if
                     end do
                  end do
                  end if
! leave these unidented to save space
                  end do
                  end if
                  end do
               end do
! done determining torsions
! ------------------------------------------------------------------------

! determine angles for iend - iu - iend
               do count = 1, ntogrow
                  iu = growlist(iw,count)

                  if (fcount(iu).gt.1) then

                     do j = 1, fcount(iu) - 1
                        ju = fclose(iu,j)
                        lengtha = flength(iu,ju)

                        do counta = j + 1, fcount(iu)

                           ku = fclose(iu,counta)

                           lengthb = flength(iu,ku)

                           thetac = (xvec(iu,ju)*xvec(iu,ku) + yvec(iu,ju)*yvec(iu,ku) + zvec(iu,ju)*zvec(iu,ku))/(lengtha*lengthb)

                           if (abs(thetac).gt.1) then
                              write(io_output,*) '*********************' ,'****************************'
                              write(io_output,*) iu,ku,xvec(iu,ku) ,yvec(iu,ku) ,zvec(iu,ku),lengthb
                              call err_exit(__FILE__,__LINE__,'thetac outsie of range',myid+1)
                           end if

                           angle = acos(thetac)

                           vphi = vphi + kforceb(ju,ku) * (angle-equilb(ju,ku))**2

! if (lshit) then
! write(io_output,*) ju,iu,ku,kforceb(ju,ku)*(angle-equilb(ju,ku))**2
! end if

                        end do
                     end do
                  end if
               end do

! ------------------------------------------------------------------------
! now lets figure out the rest of the bending for ntogrow > 1

               if (ntogrow.gt.1) then
                  do count = 1, ntogrow - 1
                     iu = growlist(iw,count)
                     do counta = count+1, ntogrow
                        ju = growlist(iw,counta)
                        if (.not.((iopen(1).eq.count .or.iopen(1).eq.counta) .and.((iclose.eq.count.or.iclose.eq.counta)&
                         .or.iopen(2).eq.count .or.iopen(2).eq.counta))) then
! we already calculated these

                           lengtha = flength(iufrom,iu)
                           lengthb = flength(iufrom,ju)

                           thetac = (xx(count)*xx(counta) + yy(count)*yy(counta) + zz(count)*zz(counta)) / (lengtha*lengthb)
                           angle = acos(thetac)

                           vphi = vphi + kforceb(iu,ju) * (angle -equilb(iu,ju))**2

! if (lshit) then
! write(io_output,*) iu,iufrom,ju,kforceb(iu,ju)*(angle-equilb(iu,ju))**2.0E0_dp
! end if
                        end if
                     end do
                  end do
               end if

               bf_tor(itor) = exp(-(vvib + vphi + vdha)*beta)
 190           continue

               vvibration(itor) = vvib
               vtorsion(itor) = vdha
               vbend(itor) = vphi
               bsum_tor(ip) = bsum_tor(ip) + bf_tor(itor)
            end do

            if (lnew.or.ip.ne.1) then
! choose one of the trial sites in a biased fashion
               ran_tor = random(-1) * bsum_tor(ip)
               bs = 0
               do itor = 1, ichtor
                  bs = bs + bf_tor(itor)
                  if (ran_tor .lt. bs ) then

! save torsion energy of this trial position
                     vtgtr(ip) = vtorsion(itor)
                     vtbend(ip) = vbend(itor) + ovphi
                     vtvib(ip) = vvibration(itor) + ovvib
! assign the phidisp of this trial position
                     do count = 1, ntogrow
                        phiacc(count) = phicrank(itor,count)
                        diracc(count) = dir(itor,count)
                     end do
! exit the loop
                     goto 200
                  end if
               end do

 200           continue
            else
! select old conformation
               vtgtr(ip) = vtorsion(1)
               vtbend(ip) = vbend(1) + ovphi
               vtvib(ip) = vvibration(1) + ovvib
            end if

! divide bsum by ichtor
            bsum_tor(ip) = bsum_tor(ip) / dble(ichtor)

            start = 1
            last = ntogrow
            ldo = .false.
            countb = 0
 60         continue

! for accepted phidisp set up the vectors
            do count = 1,ntogrow

               iu = growlist(iw,count)

               if (.not.lnew.and.ip.eq.1) then
                  length = flength(iufrom,iu)
                  xx(count) = (rxu(i,iu) - rxu(i,iufrom)) / length
                  yy(count) = (ryu(i,iu) - ryu(i,iufrom)) / length
                  zz(count) = (rzu(i,iu) - rzu(i,iufrom)) / length
               else
                  if (fcount(iu).gt.0) then
                     ju = fclose(iu,1)
                     call cone(1,xvec(iufrom,ju),yvec(iufrom,ju),zvec(iufrom,ju),dum,dum)
                     phidisp = phiacc(count)

                     call cone(2,x,y,z,alpha(count,ju),phidisp)

! store the unit vectors in xx, yy, zz
                     xx(count) = x
                     yy(count) = y
                     zz(count) = z

! else if (fcount(iu).gt.1) then

! j = diracc(count)

! if (ip.gt.1.and..not.lnew) call err_exit(__FILE__,__LINE__,'',myid+1)

! if (lnew) then
! xx(count) = rxa(count,j) - rxnew(iufrom)
! yy(count) = rya(count,j) - rynew(iufrom)
! zz(count) = rza(count,j) - rznew(iufrom)
! else
! xx(count) = rxa(count,j) - rxu(i,iufrom)
! yy(count) = rya(count,j) - ryu(i,iufrom)
! zz(count) = rza(count,j) - rzu(i,iufrom)
! end if

! xx(count) = xx(count) / length
! yy(count) = yy(count) / length
! zz(count) = zz(count) / length

                  else

                     if (.not.ldo) then
                        ldo = .true.
                        countb = countb + 1
                        iopen(countb) = count
                        goto 70
                     else if (countb.eq.1.and.iopen(1).ne.count) then
                        countb = countb + 1
                        iopen(countb) = count
                        goto 70
                     else if (countb.eq.2) then
                        countb = countb + 1
                     else
                        ldo = .false.
                        countb = countb + 1
                     end if

                     do counta = 1, ntogrow

                        if (counta.ne.count) then
                           if (opencount.gt.1) then
                              if (countb.eq.3) then
! we only want the angle with the closing bead

                                 if (iopen(2).eq.counta) then
                                    goto 65
                                 end if
                              else if (countb.eq.4) then
! we want the angle with the open bead

                                 if (iopen(1).ne.counta) then
                                    goto 65
                                 end if
                              end if
                           else if(iclose.ne.counta) then
                              goto 65
                           end if

                           lengtha = distij(iufrom,iuprev)

                           rxt(1) = xvec(iufrom,iuprev)/lengtha
                           ryt(1) = yvec(iufrom,iuprev)/lengtha
                           rzt(1) = zvec(iufrom,iuprev)/lengtha

                           rxt(2) = xx(counta)
                           ryt(2) = yy(counta)
                           rzt(2) = zz(counta)

                           angles(1) = bangles(count,1)
                           angles(2) = bangles(count,2)

                           call close(3,rxt,ryt,rzt,dum ,angles,lterm)

                           if (lterm) then
                              return
                           end if

! choose the given possibility
                           if (diracc(count).eq.1) then
                              lengtha = flength(iufrom,iu)
                              xx(count) = rxt(3)
                              yy(count) = ryt(3)
                              zz(count) = rzt(3)
                           else
                              lengtha = flength(iufrom,iu)
                              xx(count) = rxt(4)
                              yy(count) = ryt(4)
                              zz(count) = rzt(4)
                           end if
                        end if

 65                     continue

                     end do
                  end if
               end if
 70            continue
            end do

            if (ldo) then
               if (countb.eq.3) then
                  start = iopen(2)
                  last = iopen(2)
               else
                  start = iopen(1)
                  last = iopen(1)
               end if
               goto 60
            end if

! accepted coordinates, save them in r*p(trial)
            do count = 1,ntogrow
               iu = growlist(iw,count)
               length = flength(iufrom,iu)
               if ( lnew ) then
! use new positions
                  rxp(count,ip) = rxnew(iufrom) + xx(count)*length
                  ryp(count,ip) = rynew(iufrom) + yy(count)*length
                  rzp(count,ip) = rznew(iufrom) + zz(count)*length
               else
! use old coordinates
                  rxp(count,ip) = rxu(i,iufrom) + xx(count)*length
                  ryp(count,ip) = ryu(i,iufrom) + yy(count)*length
                  rzp(count,ip) = rzu(i,iufrom) + zz(count)*length
               end if
            end do
         end do

         if (lcrank) then
            vphi = vibtr
         end if
      end if

#ifdef __DEBUG__
      write(io_output,*) 'END SAFECMBC in ',myid,'. init: ',iinit
#endif
      return
  end subroutine safecbmc

!> \brief Determines logic for a CBMC Move Between Fixed End Points
!> for Linear, Branched, and Cyclic Molecules
!>
!> Presently, works for linear molecules with rigid bond
!> lengths, and can (with little program changes) work for
!> branched molecules with rigid bonds, as long as it closes
!> at a binary or tertiary segment. \n
!> Does work for any branched molecule with flexible bond
!> lengths.
!> \par NEW LOGIC ONLY FOR SAFE-CMBC
!> \b iend = beads to grow to \n
!> \b ipast = one bead past iend \n
!> \b inext = two beads past iend \n
!> \b ibef = one bead before iend \n
!> \b iwbef = two beads before iend \n
!> \b fclose(iu) = beads to calculate interaction with from iu \n
!> \b fcount(iu) = number of fcloses for iu \n
!> \b movetype = 2 swap move
!> \b          = 5 SAFE-swatch move
!> \b findex = how many backbone beads to regrow
!> \b COLLIN = MASTER of the known universe
!> \author Originally completed by Collin Wick around 1-1-2000
!******************************************************************
  subroutine safeschedule(igrow,imolty,islen,iutry,findex,movetype,iprev)
      logical::lcount,lpick,lterm,lfixed,lfix,lfind
      integer::igrow,imolty,count,counta,iw,ivib,iv,iu ,ju,iutry
      integer::j,ja,kickout,invtry,index,fintnum,fint,k ,islen
      integer::ffrom,fprev,flist,fnum,fnuma,findex ,countb,iv1
      integer::movetype,fmaxgrow,kickouta,iufrom ,iuprev
      integer::num,inum,inuma,max

      integer,optional::iprev

      parameter(max=10)
      dimension fint(numax),ffrom(numax,max),fprev(numax,max),flist(numax,max,max),fnum(numax),fnuma(numax,max),lpick(numax)&
       ,lfix(numax),inum(max),inuma(max),lfind(numax)

!     --------------------------------------------------------------------

! safecbmc scheduler ****

! now determine findex excluding 1, it won't do anything

! set begining conditions
      if (movetype.eq.2) then
         if (.not.lring(imolty)) then
            call err_exit(__FILE__,__LINE__,'you can not use safecbmc for swap unless it is a ring',myid+1)
         end if
         fmaxgrow = nunit(imolty)
      else if (movetype.eq.5) then
         if (.not. present(iprev)) then
            call err_exit(__FILE__,__LINE__,'you can not use safe-swatch without indicating iprev',myid+1)
         end if
      else
         fmaxgrow = maxgrow(imolty) + 1
      end if
      lterm = .false.
      kickout = 0

      kickouta = 0

 100  continue

      do j = 1, nunit(imolty)
         lfix(j) = .false.
      end do

      if (kickouta.eq.5) call err_exit(__FILE__,__LINE__,'',myid+1)

      if (movetype .eq. 2) then
         findex = fmaxgrow
      else if (movetype .eq. 5) then
         ! don't do anything here because for SAFE-SWATCH, findex is an input variable
      else
         findex = int( random(-1) * dble(fmaxgrow - 1)) + 2
      end if

!     ******************************
! findex = 7
!     ******************************

      lfixed = .false.

      if (findex.eq.2) then
         lcrank = .true.
      else
         lcrank = .false.
      end if

      do j = 1, igrow
         lexshed(j) = .true.
      end do

      if (kickout.gt.250) then
         call err_exit(__FILE__,__LINE__,'SAFESCHEDULE KICKED YOU OUT',myid+1)
      end if

! find iutry
      if (movetype.eq.2) then

         iutry = 1
         invtry = invib(imolty,iutry) - 1

         if (invtry.lt.1) then
            call err_exit(__FILE__,__LINE__,'You need a ring to do this',myid+1)
         end if

         ffrom(1,1) = iutry
         fprev(1,1) = 0
         do iv = 1, invtry
            flist(1,1,iv) = ijvib(imolty,iutry,iv)
            lexshed(flist(1,1,iv)) = .false.
         end do
         fnum(1) = 1
         fnuma(1,1) = invtry

      else if (movetype.eq.5) then
      ! Paul -- SAFE-SWATCH part
         ! ffrom(1,1) and fprev(1,1) are read from input ifrom and iprev
         ffrom(1,1) = iutry
         fprev(1,1) = iprev
         invtry = invib(imolty,iutry)

         if (invtry .lt. 2) then
            call err_exit(__FILE__,__LINE__,'cannot do safe-swatch starting from an end bead',myid+1)
         end if

         count = 0
         do iv = 1, invtry
            if (ijvib(imolty,iutry,iv).ne.iprev .and. .not.lplace(imolty,ijvib(imolty,iutry,iv))) then
                count = count + 1
                ! take down the bead number that is connected to iutry and is about to be regrown next
                flist(1,1,count) = ijvib(imolty,iutry,iv)
                lexshed(flist(1,1,count)) = .false.
            end if
         end do

         fnum(1) = 1
         fnuma(1,1) = count
      else

         iutry = int( random(-1) * dble(iring(imolty)+icbsta(imolty)))+1 -icbsta(imolty)! factoring in icbsta for safecbmc

!     *************************
! iutry = 1
!     *************************

         if (lplace(imolty,iutry)) then
            kickout = kickout + 1
            goto 100
         end if

         invtry = invib(imolty,iutry) ! determining the number of vibrations that bead no. iutry has

         if (invtry.eq.0) then
            call err_exit(__FILE__,__LINE__,'cant do safecbmc on single bead',myid+1)
         else if(invtry.eq.1) then  ! At the end point of a molecule
            kickout = kickout + 1  ! we will let regular cbmc handle the end points
            goto 100

            fprev(1,1) = 0
            ivib = 0
         else if (invtry.eq.3.and.maxgrow(imolty).eq.1) then
            kickout = kickout + 1 ! do not start 1 bead safecbmc from a branch point
            goto 100
         else

 13         ivib = int(random(-1) * dble(invtry)) + 1  ! at a branch point, decide which way not to grow

!     ********************************
! ivib = 2
!     *******************************

            fprev(1,1) = ijvib(imolty,iutry,ivib) ! find the previous bead of bead iutry

            if (icbdir(imolty).eq.1.and.fprev(1,1).gt.iutry) goto 13

            if (fprev(1,1).gt.iring(imolty) .or.lplace(imolty,fprev(1,1))) then
               kickout = kickout + 1
               goto 100
            end if
         end if

         ffrom(1,1) = iutry
         count = 0
         do iv = 1, invtry
            if (iv.ne.ivib.and. .not.lplace(imolty,ijvib(imolty,iutry,iv))) then
               count = count + 1
               ! take down the bead number that is connected to iutry and is about to be regrown next
               flist(1,1,count) = ijvib(imolty,iutry,iv)
               lexshed(flist(1,1,count)) = .false.
            end if
         end do

         if (count.eq.0) then
            kickout = kickout + 1
            goto 100
         end if
         fnum(1) = 1
         fnuma(1,1) = count
      end if

! find all branches going to maxgrow or end of molecule
      do iw = 2, nunit(imolty)
         count = 0
         do j = 1, fnum(iw-1)
            do ja = 1, fnuma(iw-1,j)
               counta = 0
               lcount = .false.
               iu = flist(iw-1,j,ja)
               if (.not. (lfix(iu).and.lfixed).and. .not.lrigi(imolty,iu)) then
                  do iv = 1, invib(imolty,iu)
                     ju = ijvib(imolty,iu,iv)
                     if (ju.ne.ffrom(iw-1,j).and. .not.(lplace(imolty,ju).and. iu.le.iring(imolty))) then
                        if (lfixed) then
                           if (ju.gt.iring(imolty)) then
                              counta = counta + 1
                              flist(iw,count+1,counta) = ju
                              lexshed(ju) = .false.
                              lcount = .true.
                           end if
                        else
                           counta = counta + 1
                           flist(iw,count+1,counta) = ju
                           lexshed(ju) = .false.
                           lcount = .true.
                        end if
                     end if
                  end do
                  if (lcount) then
                     count = count + 1
                     fprev(iw,count) = ffrom(iw-1,j)
                     ffrom(iw,count) = iu
                     fnuma(iw,count) = counta
                  end if
               end if
            end do
         end do

         if (count.eq.0) then
! we hit the end, lets get out of here
            index = iw - 1
            goto 110
         end if
         fnum(iw) = count
         if (iw.eq.findex) then
            lfixed = .true.
            do j = 1, fnum(findex)
               do ja = 1, fnuma(findex,j)
                  if (flist(findex,j,ja).le.iring(imolty)) then
                     lfix(flist(findex,j,ja)) = .true.
                  end if
               end do
            end do
         end if
      end do
! index = fmaxgrow
 110  continue

! don't allow 1-bead regrowths, it won't do anything
      if (index.lt.2) then
         kickout = kickout + 1
         goto 100
      end if

      if (index.lt.findex) then
         kickout = kickout + 1
         goto 100
      end if

! Paul -- don't allow the growth to stop at a branch point
      do j = 1, fnum(index)
         if (fnuma(index, j) .gt. 1) then
            if (movetype .eq. 5) then
               call err_exit(__FILE__,__LINE__,'the safe-swatch closing bead cannot be at a branch point',myid+1)
            else
                kickout = kickout + 1
                goto 100
            end if
         end if
      end do

      lfixed = .false.

! lets set logic so rosenbluth can read it
      count = 0
      do iw = 1, index

         if (iw.eq.findex) then
            lfixed = .true.
         end if

         do j = 1, fnum(iw)
            kickout = 0

            if (lfixed) then
               if (flist(iw,j,1).le.iring(imolty)) goto 122
            end if
            count = count + 1
            grownum(count) = fnuma(iw,j)
            growfrom(count) = ffrom(iw,j)
            growprev(count) = fprev(iw,j)

            do ja = 1, fnuma(iw,j)
               lpick(ja) = .false.
            end do
            do ja = 1, fnuma(iw,j)
 115           continue
               if (kickout.gt.150) then
                  call err_exit(__FILE__,__LINE__,'Randomizer in FECMBC kicked you out',myid+1)
               end if
! this is the only random part
               counta = int(random(-1) * dble(fnuma(iw,j))) + 1
               if (lpick(counta)) then
                  kickout = kickout + 1
                  goto 115
               end if
               lpick(counta) = .true.
               growlist(count,counta) = flist(iw,j,ja)
            end do
! reset our logic keep have counta coicide
            do ja = 1, fnuma(iw,j)
               flist(iw,j,ja) = growlist(count,ja)
            end do
 122        continue
         end do
      end do
      islen = count

! set all fixed points to true
      do iw = findex, index
         do count = 1, fnum(iw)
            do counta = 1, fnuma(iw,count)
               iu = flist(iw,count,counta)
               if (iu.le.iring(imolty)) then
                  lexshed(iu) = .true.
               end if
            end do
         end do
      end do

! find ends, ipasts, and inexts
      count = 0
      do j = 1, fnum(findex)
         do ja = 1, fnuma(findex,j)
            if (flist(findex,j,ja).le.iring(imolty) .and..not.lplace(imolty,flist(findex,j,ja))) then
               count = count + 1
               iend(count) = flist(findex,j,ja)
               counta = 0
               do iv = 1, invib(imolty,iend(count))
                  iu = ijvib(imolty,iend(count),iv)
                  if (iu.ne.ffrom(findex,j)) then
                     counta = counta + 1
                     ipast(iend(count),counta) = iu
                     countb = 0
                     do iv1 = 1, invib(imolty,iu)
                        ju = ijvib(imolty,iu,iv1)
                        if (ju.ne.iend(count)) then
                           countb = countb + 1
                           inext(iu,countb) = ju
                        end if
                     end do
                     nextnum(iu) = countb
                  end if
               end do
               pastnum(iend(count)) = counta
            end if
         end do
         endnum = count
      end do

! now that we found iends and ipasts,
! determine which beads to close with each iend

      do j = 1, igrow
         fcount(j) = 0
      end do

      do count = 1, endnum
         fintnum = 1
         fint(1) = iend(count)
         do iw = findex, 2, -1
            counta = 0
            do k = 1, fintnum
               do iv = 1, invib(imolty,fint(k))
                  ju = ijvib(imolty,fint(k),iv)
                  do j = 1, fnum(iw)
                     if (ju.eq.ffrom(iw,j)) then
                        fcount(ju) = fcount(ju) + 1
                        fclose(ju,fcount(ju)) = iend(count)
                        counta = counta + 1
                        fint(counta) = ju
                     end if
                  end do
               end do
            end do
         end do
         fintnum = counta
      end do

! define iwbef and ibef
      count = 0
      if (lcrank) then
         do j = 1, fnum(findex-1)
            counta = 0
            do ja = 1, fnuma(findex-1,j)
               iu = flist(findex-1,j,ja)
               if (fcount(iu).ne.0) then
                  counta = counta + 1
                  ibef(j,counta) = iu
               end if
            end do
            befnum(j) = counta
         end do
         wbefnum = fnum(findex-1)
      else
         do j = 1, fnum(findex - 1)
            iu = ffrom(findex-1,j)
            if (fcount(iu).ne.0) then
               count = count + 1
               iwbef(count) = iu
               counta = 0
               do ja = 1, fnuma(findex-1,j)
                  ju = flist(findex-1,j,ja)
                  if (fcount(ju).ne.0) then
                     counta = counta + 1
                     ibef(count,counta) = ju
                  end if
               end do
               befnum(count) = counta
            end if
         end do
         wbefnum = count
      end if

! write(io_output,*) '**********************************'
! write(io_output,*) iutry,ivib,findex,fprev(1,1)

! set up place move logic -----------

      do j = 1, nunit(imolty)
         lpnow(j) = .false.
         pnum(j) = 0
      end do

      nplace = 0

      counta = 0

      iw = 1
      do j = 1, fnum(iw)
         iufrom = ffrom(iw,j)
         if (iufrom.le.iring(imolty)) then
            counta = counta + 1
            do ivib = 1, invib(imolty,iufrom)
               iu = ijvib(imolty,iufrom,ivib)

               if (lplace(imolty,iu)) then
                  pnum(counta) = pnum(counta) + 1
                  if (pnum(counta).eq.1) then
                     nplace = nplace + 1
                  end if
                  lexshed(iu) = .false.
                  pprev(nplace) = fprev(iw,j)
                  pfrom(nplace) = iufrom
                  iplace(nplace,pnum(counta)) = iu
                  lpnow(iu) = .true.
               end if
            end do
            if (pnum(counta).eq.0) counta = counta - 1
         end if
      end do

      do iw = 1, findex - 1
         do j = 1, fnum(iw)
            iuprev = ffrom(iw,j)
            do ja = 1, fnuma(iw,j)
               iufrom = flist(iw,j,ja)
               if (iufrom.le.iring(imolty)) then
                  counta = counta + 1
                  do ivib = 1, invib(imolty,iufrom)
                     iu = ijvib(imolty,iufrom,ivib)

                     if (lplace(imolty,iu)) then
                        pnum(counta) = pnum(counta) + 1
                        if (pnum(counta).eq.1) then
                           nplace = nplace + 1
                        end if
                        lexshed(iu) = .false.
                        pprev(nplace) = iuprev
                        pfrom(nplace) = iufrom
                        iplace(nplace,pnum(counta)) = iu
                        lpnow(iu) = .true.
                     end if
                  end do
                  if (pnum(counta).eq.0) counta = counta - 1
               end if
            end do
         end do
      end do

! end place move logic setup ---------

! begin part for rig logic

      do iu = 1, nunit(imolty)
         lfind(iu) = .false.
      end do

      counta = 0
      do iw = 1, islen
         iufrom = growfrom(iw)
         do count = 1, grownum(iw)
            iu = growlist(iw,count)

            if (lrigi(imolty,iu)) then
               counta = counta + 1
               rfrom(counta) = iu
               rprev(counta) = iufrom

               lfind(iu) = .true.
               lfind(iufrom) = .true.

               ja = 0
               do ivib = 1, invib(imolty,iu)
                  ju = ijvib(imolty,iu,ivib)
                  if (ju.ne.iufrom.and..not.lpnow(ju)) then
                     ja = ja + 1
                     rlist(counta,ja) = ju
                     lfind(ju) = .true.
                     lexshed(ju) = .false.
                  end if
               end do

               if (ja.eq.0) then
                  call err_exit(__FILE__,__LINE__,'PROBLEM WITH RIG LOGIC IN SAFESCHEDULE',myid+1)
               else
                  rnum(counta) = ja
               end if
            end if

         end do
      end do

      if (counta.eq.0) then
         llrig = .false.
      else
         do iu = 1, islen
            lsave(iw) = .false.
         end do
         llrig = .true.
         nrigi = counta
      end if

      if (llrig) then
! cycle through the rest of the rigid beads to set lexshed to false

         do iw = 1, nrigi
            do count = 1, rnum(iw)
               iu = rlist(iw,count)
               inum(count) = iu
            end do

            num = rnum(iw)
 5          continue
            ja = 0

            do count = 1, num
               iu = inum(count)

               do iv = 1, invib(imolty,iu)

                  ju = ijvib(imolty,iu,iv)

                  if (.not.lfind(ju)) then

                     ja = ja + 1

                     lfind(ju) = .true.
                     inuma(ja) = ju

                     lexshed(ju) = .false.
                  end if
               end do
            end do
            if (ja.gt.10) then
               write(io_output,*) 'ja',ja
               call err_exit(__FILE__,__LINE__,'need to set max greater in safeschedule',myid+1)
            end if

            num = ja

            if (ja.ne.0) then
               do j = 1, ja
                  inum(j) = inuma(j)
               end do
               goto 5
            end if
         end do

      end if

      return

! take out return for diagnostics ************

      do iw = 1, islen
         write(io_output,*) growfrom(iw),(growlist(iw,count),count=1 ,grownum(iw))

      end do

      write(io_output,*) '---------------------------------------------'

      do iw = 1, nplace
         write(io_output,*) pfrom(iw),(iplace(iw,count),count=1,pnum(iw))
      end do

      write(io_output,*) '--------------------------------------------'

      do iw = 1, nrigi
         write(io_output,*) rfrom(iw),(rlist(iw,count),count=1,rnum(iw))
      end do

      call err_exit(__FILE__,__LINE__,'',myid+1)

!     ***************************************************

      return
  end subroutine safeschedule

!***********************************************************
!> \brief Takes three or four points and determines a point
!> that is a certain length from all of them
!>
!> This wonderful subroutine can also find a unit vector
!> connected to two other unit vectors with all the angles
!> between them given.
!>
!> \param iinit if iinit=1 it finds two possibilities to close 3 beads
!> if iinit=2 it finds one possibility to close 4 beads
!> if iinit=3 it finds a vector connected two others given
!> \attention ALL LENGTHS MUST BE THE SAME
!> \author This mess was unfortunately created by Collin Wick
!> on December 1999, BUT IT DOES WORK
!***********************************************************
  subroutine close(iinit,rx,ry,rz,bondl,angle,lterm)
      logical::lterm
      integer::iinit
      real::x,y,z,rx,ry,rz,length,lengtha,lengthb,xa,ya,za,theta,thetac,ux,uy,uz,bondl,avar,bvar,cvar,rxa,rya,rza,lengthc,angle&
       ,rxf,ryf,rzf,var,dvar,a,bb,c

      dimension rx(6),ry(6),rz(6),x(4),y(4),z(4),ux(3),uy(3),uz(3)
      dimension angle(3)
! ---------------------------------------------------------------------
! Determines a point from three others with equal bond lengths

#ifdef __DEBUG__
      write(io_output,*) 'START CLOSE in ',myid,'. iinit:',iinit
#endif

      if (iinit.ne.3) then

         if (iinit.eq.2) then
            rxf = rx(4)
            ryf = ry(4)
            rzf = rz(4)
         end if

! determine lengths from all coordinates
         lengtha = sqrt( (rx(2)-rx(1))**2 + (ry(2)-ry(1))**2 + (rz(2)-rz(1))**2)
         lengthb = sqrt( (rx(2)-rx(3))**2 + (ry(2)-ry(3))**2 + (rz(2)-rz(3))**2)
         lengthc = sqrt( (rx(3)-rx(1))**2 + (ry(3)-ry(1))**2 + (rz(3)-rz(1))**2)

! detemine the angle at 1
         thetac = - 0.5E0_dp * ( lengthb**2 - lengtha**2 - lengthc**2 ) / ( lengtha * lengthc )
         theta = acos(thetac)

! now lets determine the in-plane coordinates
         x(1) = 0
         y(1) = 0
         x(2) = lengtha
         y(2) = 0
         x(3) = lengthc * thetac
         y(3) = lengthc * sin( theta )

! determine the in-plane coordinates with the same length
         ya = 0.5E0_dp * ( (x(2)**2 - x(3)**2 - y(3)**2) / (-y(3)) + (x(2)**2 * (x(3)-x(2))) / ((- y(3))*x(2)) )

         xa = 0.5E0_dp * x(2)

! determine length to new position
         length = sqrt( (x(1)-xa)**2 + (y(1)-ya)**2 )

         if (length.gt.bondl) then
! the distances are too far to be able to close
            lterm = .true.
            return
         end if

! determine perpendicular distance from here to final point
         lengthc = sqrt( bondl**2 - length**2 )

! find vectors from 1 to 2, 2 to 3, and 3 to 1 in real::space
         x(1) = (rx(2) - rx(1))
         y(1) = (ry(2) - ry(1))
         z(1) = (rz(2) - rz(1))

         x(2) = (rx(3) - rx(2))
         y(2) = (ry(3) - ry(2))
         z(2) = (rz(3) - rz(2))

         x(3) = (rx(1) - rx(3))
         y(3) = (ry(1) - ry(3))
         z(3) = (rz(1) - rz(3))

! find middle point from 1 to 2

         rx(4) = 0.5E0_dp * x(1) + rx(1)
         ry(4) = 0.5E0_dp * y(1) + ry(1)
         rz(4) = 0.5E0_dp * z(1) + rz(1)

! find distance from this point to a

         length = sqrt( length**2 - (0.5E0_dp*lengtha)**2 )

! cross 1 with 2 to find perpendicular vector

         ux(1) = y(1)*z(2) - z(1)*y(2)
         uy(1) = z(1)*x(2) - x(1)*z(2)
         uz(1) = x(1)*y(2) - y(1)*x(2)

! cross this vector with 1 to find vector to final spot

         ux(2) = y(1)*uz(1) - z(1)*uy(1)
         uy(2) = z(1)*ux(1) - x(1)*uz(1)
         uz(2) = x(1)*uy(1) - y(1)*ux(1)

! normalize this vector and give it the appropriate distance

         lengtha = sqrt( ux(2)**2 + uy(2)**2 + uz(2)**2)

         ux(2) = length * ux(2) / lengtha
         uy(2) = length * uy(2) / lengtha
         uz(2) = length * uz(2) / lengtha

! now determine the two possible points here

         rx(5) = rx(4) + ux(2)
         ry(5) = ry(4) + uy(2)
         rz(5) = rz(4) + uz(2)

         rx(6) = rx(4) - ux(2)
         ry(6) = ry(4) - uy(2)
         rz(6) = rz(4) - uz(2)

! with the two possibilities, determine which is closer to 3

         lengtha = sqrt((rx(3)-rx(5))**2 + (ry(3)-ry(5))**2 + (rz(3) - rz(5))**2)
         lengthb = sqrt((rx(3)-rx(6))**2 + (ry(3)-ry(6))**2 + (rz(3) - rz(6))**2)

         if (lengtha.lt.lengthb) then
! number 5 is right
            rxa = rx(5)
            rya = ry(5)
            rza = rz(5)
         else
! number 6 is right
            rxa = rx(6)
            rya = ry(6)
            rza = rz(6)
         end if

! find vector from 1 to a

         xa = rxa - rx(1)
         ya = rya - ry(1)
         za = rza - rz(1)

! cross vector a with 1

         ux(1) = ya*z(1) - za*y(1)
         uy(1) = za*x(1) - xa*z(1)
         uz(1) = xa*y(1) - ya*x(1)

! normalize this and give it the length to our final point

         length = sqrt(ux(1)**2 + uy(1)**2 + uz(1)**2)

         ux(1) = ux(1) * lengthc / length
         uy(1) = uy(1) * lengthc / length
         uz(1) = uz(1) * lengthc / length

! finally lets find our final points and send them back

         rx(1) = rxa + ux(1)
         ry(1) = rya + uy(1)
         rz(1) = rza + uz(1)

         rx(2) = rxa - ux(1)
         ry(2) = rya - uy(1)
         rz(2) = rza - uz(1)

         if (iinit.eq.2) then

! for four bead closes, only one of these will work

            x(1) = rxf - rx(1)
            y(1) = ryf - ry(1)
            z(1) = rzf - rz(1)
            lengtha = sqrt(x(1)**2 + y(1)**2 + z(1)**2)

            x(2) = rxf - rx(2)
            y(2) = ryf - ry(2)
            z(2) = rzf - rz(2)
            lengthb = sqrt(x(2)**2 + y(2)**2 + z(2)**2)

            write(io_output,*) lengtha,lengthb,bondl

            if ((lengtha-bondl).lt.(lengthb-bondl)) then
               rx(1) = rx(1)
               ry(1) = ry(1)
               rz(1) = rz(1)
               if (abs(lengtha-bondl).gt.0.0001) then
                  lterm = .true.
                  return
               end if
            else
               if (abs(lengthb-bondl).gt.0.0001) then
                  lterm = .true.
                  return
               end if

               rx(1) = rx(2)
               ry(1) = ry(2)
               rz(1) = rz(2)
            end if
         end if

      else
! lets find our open vector with angles from prev and closed

! angle 1 is from vector 1 to the new
! angle 2 is from vector 2 to the new
! angle 3 is between vectors 1 and 2

         avar = cos(angle(2)) - rx(2)*cos(angle(1))/rx(1)
         bvar = rx(2) * rz(1) / rx(1) - rz(2)
         cvar = ry(2) - rx(2) * ry(1) / rx(1)

         avar = avar / cvar
         bvar = bvar / cvar

         cvar = (ry(1) * bvar + rz(1)) / rx(1)
         dvar = (cos(angle(1)) - ry(1) * avar) / rx(1)

         a = cvar**2 + bvar**2 + 1.0E0_dp
         bb = 2.0E0_dp * (avar*bvar - cvar*dvar)
         c = dvar**2 + avar**2 - 1.0E0_dp

         var = bb**2 - 4.0E0_dp * a * c

         if (var.lt.0) then
            lterm = .true.
            return
         end if

         rz(3) = (-bb + sqrt(var)) / (2.0E0_dp * a )
         rz(4) = (-bb - sqrt(var)) / (2.0E0_dp * a )

         rx(3) = dvar - rz(3) * cvar
         rx(4) = dvar - rz(4) * cvar

         ry(3) = rz(3) * bvar + avar
         ry(4) = rz(4) * bvar + avar

         return

         avar = cos(angle(1))*(rx(2)*(rx(2)*rz(1) - rx(1)*rz(2)) + ry(2)*(ry(2)*rz(1) - ry(1)*rz(2)))&
          + cos(angle(2))*(rx(1)*(rx(1)*rz(2) - rx(2)*rz(1)) + ry(1)*(ry(1)*rz(2) - ry(2)*rz(1)))

         var = rx(2)**2*(ry(1)**2 + rz(1)**2) + ry(2)**2*(rx(1)**2 + rz(1)**2) + rz(2)**2*(rx(1)**2 + ry(1)**2)&
          - 2.0E0_dp*(rx(1)*rx(2)*ry(1)*ry(2) + rx(1)*rx(2)*rz(1)*rz(2) + ry(1)*ry(2)*rz(1)*rz(2))&
          - (rx(2)**2 + ry(2)**2 + rz(2)**2)*cos(angle(1))**2 - (rx(1)**2 + ry(1)**2 + rz(1)**2)&
          *cos(angle(2))**2 + 2.0E0_dp*(rx(1)*rx(2) + ry(1)*ry(2) + rz(1)*rz(2)) *cos(angle(1))*cos(angle(2))

         if (var.lt.0) then
            var = abs(var)
            lterm = .true.

! return
         end if

         bvar = (rx(1)*ry(2) - rx(2)*ry(1))*sqrt(var)

         cvar = rx(2)**2*(ry(1)**2 + rz(1)**2) + (ry(2)*rz(1) - ry(1)*rz(2))**2 - 2.0E0_dp*rx(1)*rx(2)*(ry(1)*ry(2)&
          + rz(1)*rz(2)) + rx(1)**2*(ry(2)**2 + rz(2)**2)

         rz(3) = (avar + bvar) / cvar
         rz(4) = (avar - bvar) / cvar

         avar = rx(2)*rz(1) - rx(1)*rz(2)
         bvar = rx(1)*cos(angle(2)) - rx(2)*cos(angle(1))
         cvar = rx(1)*ry(2) - rx(2)*ry(1)

         ry(3) = (rz(3)*avar + bvar) / cvar
         ry(4) = (rz(4)*avar + bvar) / cvar

         rx(3) = (cos(angle(1)) - ry(1)*ry(3) - rz(1)*rz(3)) / rx(1)
         rx(4) = (cos(angle(1)) - ry(1)*ry(4) - rz(1)*rz(4)) / rx(1)

      end if

#ifdef __DEBUG__
      write(io_output,*) 'END CLOSE in ',myid,'. iinit:',iinit
#endif

      return
  end subroutine close

  subroutine rigfix(lnew,i,ibox,imolty,lterm,wrig)
    logical::lnew,ovrlap,lterm,lovra,lfind,lshit
    integer::iw,i,ibox,imolty,iufrom,iuprev,ntogrow,count,iu,counta,ilist,ja,max,num,inum,j,ju,iv,nlist,ichoi,ichtor,ip,itor&
     ,it,jut2,jut3,jut4,iwalk,glist,ifrom,inuma

    parameter(max=10)

    real::xub,yub,zub,lengtha,lengthb,dum,xfix,yfix,zfix,phia,bendang,thetac,phidisp,phi,rlength,vdha,vtorsion,phitors,bf_tor&
     ,ran_tor,bs,rxpa,rypa,rzpa,bsuma,vtrya,vtrintraa,vtrexta,vtrelecta,vtrewalda,vtrorienta,vtrintera,bsum,rbf,wrig&
     ,vtrelecta_intra,vtrelecta_inter

    dimension ilist(numax),inum(max),xfix(numax),yfix(numax),zfix(numax),lfind(numax),phia(numax),bendang(numax),rlength(numax)&
     ,vtorsion(nchtor_max),phitors(nchtor_max),bf_tor(nchtor_max),rxpa(numax,nchmax),rypa(numax,nchmax),rzpa(numax,nchmax)&
     ,vtrelecta(nchmax),vtrewalda(nchmax),bsuma(nchmax),vtrya(nchmax),vtrintraa(nchmax),vtrorienta(nchmax),vtrexta(nchmax)&
     ,glist(numax),lovra(nchmax),vtrintera(nchmax),ifrom(numax),inuma(max),vtrelecta_intra(nchmax),vtrelecta_inter(nchmax)
!     ----------------------------------------------------------
#ifdef __DEBUG__
      write(io_output,*) 'START RIGFIX in ',myid
#endif

      wrig = 1.0E0_dp
      do j = 1, nunit(imolty)
         lfind(j) = .false.
      end do

      ichoi = nchoi(imolty)
      ichtor = nchtor(imolty)

      do iw = 1, nrigi
         iufrom = rfrom(iw)
         iuprev = rprev(iw)
         ntogrow = rnum(iw)

         lfind(iufrom) = .true.

! we must first set up cone for old configuration
         xfix(iufrom) = rxu(i,iufrom) - rxu(i,iuprev)
         yfix(iufrom) = ryu(i,iufrom) - ryu(i,iuprev)
         zfix(iufrom) = rzu(i,iufrom) - rzu(i,iuprev)

         lengthb = sqrt(xfix(iufrom)**2 + yfix(iufrom)**2  + zfix(iufrom)**2)

         xub = xfix(iufrom) / lengthb
         yub = yfix(iufrom) / lengthb
         zub = zfix(iufrom) / lengthb

         call cone(1,xub,yub,zub,dum,dum)

! now we must cycle through all sites that we want to be rigid
         do count = 1, ntogrow
            iu = rlist(iw,count)

            xfix(iu) = rxu(i,iu) - rxu(i,iufrom)
            yfix(iu) = ryu(i,iu) - ryu(i,iufrom)
            zfix(iu) = rzu(i,iu) - rzu(i,iufrom)

            lengtha = sqrt(xfix(iu)**2 + yfix(iu)**2 + zfix(iu)**2)

            rlength(iu) = lengtha

            thetac = -(xfix(iu)*xfix(iufrom) + yfix(iu)*yfix(iufrom) + zfix(iu)*zfix(iufrom)) / (lengtha*lengthb)

            if (abs(thetac).gt.1.0E0_dp) call err_exit(__FILE__,__LINE__,'screwup in rigfix',myid+1)

            bendang(iu) = acos(thetac)

            xub = xfix(iu) / lengtha
            yub = yfix(iu) / lengtha
            zub = zfix(iu) / lengtha

! determine the phi value associated with this
            call cone(3,xub,yub,zub,bendang(iu),phia(iu))

            lfind(iu) = .true.
            inum(count) = iu
         end do

         counta = 0
         num = ntogrow
 5       continue
         ja = 0
         do count = 1, num
            iu = inum(count)

            do iv = 1, invib(imolty,iu)

               ju = ijvib(imolty,iu,iv)

               if (.not.lfind(ju)) then
                  ja = ja + 1
                  counta = counta + 1
                  ilist(counta) = ju

                  ifrom(counta) = iu
                  lfind(ju) = .true.
                  inuma(ja) = ju

                  xfix(ju) = rxu(i,ju) - rxu(i,iufrom)
                  yfix(ju) = ryu(i,ju) - ryu(i,iufrom)
                  zfix(ju) = rzu(i,ju) - rzu(i,iufrom)

                  lengtha = sqrt(xfix(ju)**2 + yfix(ju)**2 + zfix(ju)**2)

                  rlength(ju) = lengtha

                  thetac = -(xfix(ju)*xfix(iufrom) + yfix(ju) *yfix(iufrom) + zfix(ju)*zfix(iufrom)) / (lengtha*lengthb)

                  bendang(ju) = acos(thetac)

                  xub = xfix(ju) / lengtha
                  yub = yfix(ju) / lengtha
                  zub = zfix(ju) / lengtha
! determine the phi value associated with this
                  call cone(3,xub,yub,zub,bendang(ju),phia(ju))

               end if
            end do
         end do
         num = ja

         if (ja.ne.0) then
            do j = 1, ja
               inum(j) = inuma(j)
            end do
            goto 5
         end if
         nlist = counta

! now that we determined the phi values for the old configuration,
! let's set up cone for the new configuration

         if (lnew) then
            xub = xvec(iuprev,iufrom)
            yub = yvec(iuprev,iufrom)
            zub = zvec(iuprev,iufrom)

            lengthb = distij(iuprev,iufrom)

            xub = xub / lengthb
            yub = yub / lengthb
            zub = zub / lengthb

            call cone(1,xub,yub,zub,dum,dum)
         end if

         do ip = 1, ichoi
            bsum_tor(ip) = 0.0E0_dp

            do itor = 1, ichtor
               vdha = 0.0E0_dp

               lshit = .false.
               if (lnew.and.ip.eq.17.and.itor.eq.6) then
                  lshit = .true.
               end if

               if (.not.lnew.and.ip.eq.1.and.itor.eq.1) then
                  lshit = .true.
               end if

               do count = 1, ntogrow
                  iu = rlist(iw,count)
                  if (.not.lnew.and.itor.eq.1.and.ip.eq.1) then
                     xub = rxu(i,iu) - rxu(i,iufrom)
                     yub = ryu(i,iu) - ryu(i,iufrom)
                     zub = rzu(i,iu) - rzu(i,iufrom)
                  else
                     phidisp = random(-1) * twopi

                     phi = phia(iu) + phidisp

                     call cone(2,xub,yub,zub,bendang(iu),phi)

                     lengtha = rlength(iu)

                     xub = xub * lengtha
                     yub = yub * lengtha
                     zub = zub * lengtha

                  end if

                  xvec(iufrom,iu) = xub
                  yvec(iufrom,iu) = yub
                  zvec(iufrom,iu) = zub

                  xvec(iu,iufrom) = -xub
                  yvec(iu,iufrom) = -yub
                  zvec(iu,iufrom) = -zub

                  do it = 1, intor(imolty,iu)
                     jut2 = ijtor2(imolty,iu,it)
                     jut3 = ijtor3(imolty,iu,it)
                     jut4 = ijtor4(imolty,iu,it)

                     if (jut2.eq.iufrom.and.jut3.eq.iuprev .and..not.lplace(imolty,jut4)) then
! check to see if jut4 exists
                        if (.not. lexist(jut4)) then
                           write(io_output,*) 'iu,jut2,jut3,jut4',iu ,jut2,jut3,jut4
                           call err_exit(__FILE__,__LINE__,'trouble, jut4 does not exist in rigfix',myid+1)
                        end if
                        vdha = vdha + vtorso(xvec(iu,jut2),yvec(iu,jut2),zvec(iu,jut2),xvec(jut2,jut3),yvec(jut2,jut3)&
                         ,zvec(jut2,jut3),xvec(jut3,jut4),yvec(jut3,jut4),zvec(jut3,jut4),ittor(imolty,iu,it))
                     end if
                  end do
               end do

               vtorsion(itor) = vdha
               phitors(itor) = phidisp
               bf_tor(itor) = exp(-beta*vdha)
               bsum_tor(ip) = bsum_tor(ip) + bf_tor(itor)
            end do

! pick a torsion at random
            if (lnew .or. ip .ne. 1) then
               ran_tor = random(-1)*bsum_tor(ip)
               bs = 0.0E0_dp
               do itor = 1, ichtor
                  bs = bs + bf_tor(itor)
                  if (ran_tor .lt. bs) then
                     vtgtr(ip) = vtorsion(itor)
                     phidisp = phitors(itor)
                     goto 100
                  end if
               end do
            else
               vtgtr(ip) = vtorsion(1)
               phidisp = phitors(1)
            end if

 100        continue

            bsum_tor(ip) = bsum_tor(ip) / dble(ichtor)

            do count = 1, ntogrow
               iu = rlist(iw,count)
               if (lnew .or. ip.ne.1) then
                  phi = phia(iu) + phidisp

                  call cone(2,xub,yub,zub,bendang(iu),phi)

                  lengtha = rlength(iu)

                  xub = xub * lengtha
                  yub = yub * lengtha
                  zub = zub * lengtha

                  if (lnew) then
                     rxpa(iu,ip) = xub + rxnew(iufrom)
                     rypa(iu,ip) = yub + rynew(iufrom)
                     rzpa(iu,ip) = zub + rznew(iufrom)
                  else
                     rxpa(iu,ip) = xub + rxu(i,iufrom)
                     rypa(iu,ip) = yub + ryu(i,iufrom)
                     rzpa(iu,ip) = zub + rzu(i,iufrom)
                  end if
               else
                  rxpa(iu,ip) = rxu(i,iu)
                  rypa(iu,ip) = ryu(i,iu)
                  rzpa(iu,ip) = rzu(i,iu)
               end if
            end do

! we must determine the positions of the remaining sites
            do counta = 1, nlist
               iu = ilist(counta)

               if (lnew.or.ip.ne.1) then
                  phi = phia(iu) + phidisp

                  call cone(2,xub,yub,zub,bendang(iu),phi)

                  lengtha = rlength(iu)

                  xub = xub * lengtha
                  yub = yub * lengtha
                  zub = zub * lengtha

                  if (lnew) then
                     rxpa(iu,ip) = xub + rxnew(iufrom)
                     rypa(iu,ip) = yub + rynew(iufrom)
                     rzpa(iu,ip) = zub + rznew(iufrom)
                  else
                     rxpa(iu,ip) = xub + rxu(i,iufrom)
                     rypa(iu,ip) = yub + ryu(i,iufrom)
                     rzpa(iu,ip) = zub + rzu(i,iufrom)
                  end if
               else
                  rxpa(iu,ip) = rxu(i,iu)
                  rypa(iu,ip) = ryu(i,iu)
                  rzpa(iu,ip) = rzu(i,iu)
               end if

            end do

         end do

! now calculate the intramolecular energies of the beads

! initialize rosenbluth weight
         do ip = 1, ichoi
            bsuma(ip) = 1.0E0_dp
            vtrya(ip) = 0.0E0_dp
            lovra(ip) = .false.
            vtrintraa(ip) = 0.0E0_dp
            vtrexta(ip)   = 0.0E0_dp
            vtrintera(ip) = 0.0E0_dp
            vtrelecta(ip) =  0.0E0_dp
            vtrelecta_intra(ip) =  0.0E0_dp
            vtrelecta_inter(ip) =  0.0E0_dp
            vtrewalda(ip) = 0.0E0_dp
            vtrorienta(ip) = 0.0E0_dp
         end do

         do count = 1, ntogrow
            iu = rlist(iw,count)

            do ip = 1, ichoi
               rxp(count,ip) = rxpa(iu,ip)
               ryp(count,ip) = rypa(iu,ip)
               rzp(count,ip) = rzpa(iu,ip)
            end do
            glist(count) = iu
         end do

         call boltz(lnew,.false.,ovrlap,i,i,imolty,ibox,ichoi ,iufrom,ntogrow,glist,dum)

         if (ovrlap) then
            lterm = .true.
            return
         end if

! propagate rosenbluth weigth and energies
         do ip = 1, ichoi
            if (lovr(ip)) then
               lovra(ip) = .true.
            end if
            bsuma(ip) = bsuma(ip) * bfac(ip)
            vtrya(ip) = vtrya(ip) + vtr(ivTot,ip)
            vtrintraa(ip) = vtrintraa(ip) + vtr(ivIntraLJ,ip)
            vtrexta(ip)   = vtrexta(ip) + vtr(ivExt,ip)
            vtrintera(ip) = vtrintera(ip) + vtr(ivInterLJ,ip)
            vtrelecta(ip) =  vtrelecta(ip) + vtr(ivElect,ip)
            vtrelecta_intra(ip) =  vtrelecta_intra(ip) + vtrelect_intra(ip)
            vtrelecta_inter(ip) =  vtrelecta_inter(ip) + vtrelect_inter(ip)
            vtrewalda(ip) = vtrewalda(ip) + vtr(ivEwald,ip)
            vtrorienta(ip) = vtrorienta(ip) + vtrorient(ip)
         end do

! now run through all other sites
         do counta = 1, nlist
            count = 1
            iu = ilist(counta)

            do ip = 1, ichoi
               rxp(count,ip) = rxpa(iu,ip)
               ryp(count,ip) = rypa(iu,ip)
               rzp(count,ip) = rzpa(iu,ip)
            end do

            glist(count) = iu

            call boltz(lnew,.false.,ovrlap,i,i,imolty,ibox,ichoi,ifrom(counta),1,glist,dum)

            if (ovrlap) then
               lterm = .true.
               return
            end if

! propagate rosenbluth weigth and energies
            do ip = 1, ichoi
               if (lovr(ip)) then
                  lovra(ip) = .true.
               end if
               bsuma(ip) = bsuma(ip) * bfac(ip)
               vtrya(ip) = vtrya(ip) + vtr(ivTot,ip)
               vtrintraa(ip) = vtrintraa(ip) + vtr(ivIntraLJ,ip)
               vtrexta(ip)   = vtrexta(ip) + vtr(ivExt,ip)
               vtrintera(ip) = vtrintera(ip) + vtr(ivInterLJ,ip)
               vtrelecta(ip) =  vtrelecta(ip) + vtr(ivElect,ip)
               vtrelecta_intra(ip) =  vtrelecta_intra(ip) + vtrelect_intra(ip)
               vtrelecta_inter(ip) =  vtrelecta_inter(ip) +  vtrelect_inter(ip)
               vtrewalda(ip) = vtrewalda(ip) + vtr(ivEwald,ip)
               vtrorienta(ip) = vtrorienta(ip) + vtrorient(ip)
            end do
         end do

         bsum = 0.0E0_dp
! add up rosenbluth weight
         do ip = 1, ichoi
            if (.not. lovra(ip)) then
               bsum = bsum + bsuma(ip) * bsum_tor(ip)
            end if
         end do

         if (lnew) then
            wrig = wrig * bsum

            if (wrig .lt. softlog) then
               lterm = .true.
               return
            end if

            rbf = bsum * random(-1)
            bs = 0.0E0_dp
            do ip = 1, ichoi
               if ( .not. lovra(ip) ) then
                  bs = bs + bsuma(ip) * bsum_tor(ip)
                  if (rbf .lt. bs ) then
                     iwalk = ip
                     goto 120
                  end if
               end if
            end do
 120        continue
         else
            wrig = wrig * bsum
            if (wrig .lt. softlog) then
               lterm = .true.
               write(io_output,*) 'RIGFIX OLD REJECTED'
               return
            end if
         end if

! now we must add up energies and record new positions
         if (lnew) then
            vnew(ivTot) = vnew(ivTot) + vtrya(iwalk) + vtgtr(iwalk)
            vnew(ivTorsion) = vnew(ivTorsion) + vtgtr(iwalk)
            vnew(ivExt)   = vnew(ivExt)   + vtr(ivExt,iwalk)
            vnew(ivIntraLJ) = vnew(ivIntraLJ) + vtrintraa(iwalk)
            vnew(ivInterLJ) = vnew(ivInterLJ) + vtrintera(iwalk)
            vnew(ivElect) = vnew(ivElect) + vtrelecta(iwalk)
            vnew(ivEwald) = vnew(ivEwald) + vtrewalda(iwalk)
            vneworient = vneworient + vtrorienta(iwalk)
         else
            vold(ivTot) = vold(ivTot) + vtrya(1) + vtgtr(1)
            vold(ivTorsion) = vold(ivTorsion) + vtgtr(1)
            vold(ivExt)   = vold(ivExt)   + vtr(ivExt,1)
            vold(ivIntraLJ) = vold(ivIntraLJ) + vtrintraa(1)
            vold(ivInterLJ) = vold(ivInterLJ) + vtrintera(1)
            vold(ivElect) = vold(ivElect) + vtrelecta(1)
            vold(ivEwald) = vold(ivEwald) + vtrewalda(1)
            voldorient = voldorient + vtrorienta(1)
         end if

         do count = 1, ntogrow
            iu = rlist(iw,count)

            lexist(iu) = .true.

            if (lnew) then
               rxnew(iu) = rxpa(iu,iwalk)
               rynew(iu) = rypa(iu,iwalk)
               rznew(iu) = rzpa(iu,iwalk)
            end if
            ju = iufrom

            if (lnew) then
               xvec(iu,ju) = rxnew(ju) - rxnew(iu)
               yvec(iu,ju) = rynew(ju) - rynew(iu)
               zvec(iu,ju) = rznew(ju) - rznew(iu)
            else
               xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
               yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
               zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
            end if

            distij(iu,ju) = sqrt(xvec(iu,ju)**2 + yvec(iu,ju)**2 + zvec(iu,ju)**2)

            distij(ju,iu) = distij(iu,ju)

            xvec(ju,iu) = -xvec(iu,ju)
            yvec(ju,iu) = -yvec(iu,ju)
            zvec(ju,iu) = -zvec(iu,ju)

         end do

         do counta = 1, nlist
            iu = ilist(counta)

            lexist(iu) = .true.

            if (lnew) then
               rxnew(iu) = rxpa(iu,iwalk)
               rynew(iu) = rypa(iu,iwalk)
               rznew(iu) = rzpa(iu,iwalk)
            end if

            ju = ifrom(counta)

            if (lnew) then
               xvec(iu,ju) = rxnew(ju) - rxnew(iu)
               yvec(iu,ju) = rynew(ju) - rynew(iu)
               zvec(iu,ju) = rznew(ju) - rznew(iu)
            else
               xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
               yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
               zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
            end if

            distij(iu,ju) = sqrt(xvec(iu,ju)**2 + yvec(iu,ju)**2 + zvec(iu,ju)**2)

            distij(ju,iu) = distij(iu,ju)

            xvec(ju,iu) = -xvec(iu,ju)
            yvec(ju,iu) = -yvec(iu,ju)
            zvec(ju,iu) = -zvec(iu,ju)

         end do
      end do

#ifdef __DEBUG__
      write(io_output,*) 'END RIGFIX in ',myid
#endif

      return
  end subroutine rigfix

  !> \brief Performs group CBMC
  !>
  !> \param lnew true for new configurations
  !> \param lterm true if early terminated
  !> \param i perform rosenbluth growth for chain i
  !> \param imolty molecule type of chain i
  !> \param ibox box number of chain i
  subroutine group_cbmc_grow(lnew,lterm,i,imolty,ibox)
      use util_mp,only:mp_set_displs,mp_allgather
      use energy_intramolecular,only:U_bonded

      logical::lnew,lterm
      integer::i,imolty,ibox

      logical::l_reach_end,ovrlap,lexist_temp(numax),lexshed_original(numax)
      integer::gcbmc_imolty,i_grow_unit,unit_num_local,repeat_unit_list(numax),glist(numax)
      integer::iufrom,iuprev,idir,unit_index,repeat_unit_imolty,iw,igrow,i_repeat_unit,repeat_unit_nunit
      integer::ichoi,ibead,jbead,itor,ichtor,ip,iu,iv,ju,it,jut2,jut3,jut4
      integer::first_bead,second_bead,bead_count
      integer::growfrom_mol(numax),growprev_mol(numax),grownum_mol(numax),growlist_mol(numax,numax) !< growlist of the original molecule
      real::accum_prob,rand_tor,rand_nonb,bsum
      real :: vnew_torsion,phi,maxlen,length,dchain,vvib,vbend,vtg
      real :: xtarget(3), ytarget(3), ztarget(3)
      real :: repeat_unit_energy(nEnergy)
      real :: repeat_unit_rx(numax), repeat_unit_ry(numax), v_repeat_unit(nEnergy)
      real :: repeat_unit_rz(numax)
      real :: repeat_unit_rxp(numax), repeat_unit_ryp(numax), repeat_unit_rzp(numax), bead1_coord(3),bead2_coord(3)
      real, allocatable :: rxnew_temp(:), rynew_temp(:), rznew_temp(:)

      ! MPI
      integer :: rcounts(numprocs),displs(numprocs),my_start,my_end,blocksize,my_itrial,rid
      real, allocatable :: my_bf_tor(:), my_vtorsion(:), my_phitors(:), my_vtorsion_procs(:)
      real, allocatable :: bf_tor(:), vtorsion(:),phitors(:)

      iw = 1
      l_reach_end = .false.
      lterm = .false.
      igrow = nugrow(imolty)
      gcbmc_imolty = gcbmc_mol_list(imolty)
      ichtor = nchtor(imolty)
      ichoi = nchoi(imolty)
      unit_num_local = gcbmc_unit_num(gcbmc_imolty)
      allocate(rxnew_temp(igrow), rynew_temp(igrow), rznew_temp(igrow))
      allocate(my_bf_tor(ichtor), my_vtorsion(ichtor), my_phitors(ichtor), my_vtorsion_procs(ichtor),&
         bf_tor(ichtor), vtorsion(ichtor), phitors(ichtor))

      growprev_mol = growprev
      growfrom_mol = growfrom
      grownum_mol = grownum
      growlist_mol = growlist

      rxnew_temp(1:igrow) = rxnew(1:igrow)
      rynew_temp(1:igrow) = rynew(1:igrow)
      rznew_temp(1:igrow) = rznew(1:igrow)
      lexshed_original(1:igrow) = lexshed(1:igrow)
      lexist_temp(1:igrow) = lexshed(1:igrow)
      lexist(1:igrow) = lexist_temp(1:igrow)

      if (lnew) then
          weight = 1.0E0_dp
      else
          weiold = 1.0E0_dp
      end if

      ! calculate the bond vectors that already exist and store in *vec_temp
      do iu = 1, igrow
          do iv = 1, invib(imolty,iu)
              ju = ijvib(imolty,iu,iv)
              if (lexist(iu) .and. lexist(ju)) then
                  xvec(iu, ju) = rxnew(ju) - rxnew(iu)
                  yvec(iu, ju) = rynew(ju) - rynew(iu)
                  zvec(iu, ju) = rznew(ju) - rznew(iu)
                  distij(iu,ju) = sqrt(xvec(iu,ju)**2 + yvec(iu,ju)**2 + zvec(iu,ju)**2)
              end if
          end do
      end do

      ! MPI
      if (numprocs.gt.1) then
          rid=myid
      else
          rid=-1
      end if
      blocksize = ichtor/numprocs
      rcounts = blocksize
      blocksize = ichtor - blocksize * numprocs
      if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
      call mp_set_displs(rcounts,displs,blocksize,numprocs)
      my_start = displs(myid+1) + 1
      my_end = my_start + rcounts(myid+1) - 1

      ! loop over all segments to grow
      do while (.true.)
         iufrom = growfrom_mol(iw)
         iuprev = growprev_mol(iw)

         ! determine which segment to grow (i_grow_unit)
         ! TO BE ADDED: the growth of the first and last segment alone
         do unit_index = 1, unit_num_local
             if ((iufrom .eq. gcbmc_unit_list(gcbmc_imolty, unit_index, 2)) .and. &
                 (iuprev .eq. gcbmc_unit_list(gcbmc_imolty, unit_index, 1))) then
                 idir = 1
                 i_grow_unit = unit_index
                 goto 925
             else if ((unit_index .gt. 1) .and. (iufrom .eq. gcbmc_unit_list(gcbmc_imolty, unit_index, 1)) .and. &
                 (iuprev .eq. gcbmc_unit_list(gcbmc_imolty, unit_index, 2))) then
                 idir = -1
                 i_grow_unit = unit_index - 1
                 goto 925
             else if (iufrom .eq. 1 .and. iuprev .eq. 0) then
                 idir = 1
                 i_grow_unit = 1
                 goto 925
             else if (iufrom .eq. igrow .and. iuprev .eq. 0) then
                 idir = -1
                 i_grow_unit = unit_num_local
                 goto 925
             end if
         end do

         ! this is not the bead to start the segment, continue
         iw = iw + 1
         cycle

925      continue

         if ((idir .eq. -1 .and. i_grow_unit .eq. 1) .or. &
            (idir .eq. 1 .and. i_grow_unit .eq. unit_num_local)) then
             l_reach_end = .true.
         end if

         repeat_unit_imolty = gcbmc_unit_moltype(gcbmc_imolty, i_grow_unit)
         repeat_unit_nunit = nugrow(repeat_unit_imolty)

         ! localize the correspondence between ibead in repeat_unit list and the original molecule
         do unit_index = 1, repeat_unit_nunit
             ibead = gcbmc_unit_list(gcbmc_imolty, i_grow_unit, unit_index)
             repeat_unit_list(unit_index) = ibead
         end do

         !---choose a repeat unit from the reservoir (or the old one)
         if (lnew) then
             ! pick a repeat unit from the reservoir RANDOMLY
             dchain = real(temtyp(repeat_unit_imolty),dp)
             i_repeat_unit = int( dchain*random(-1) + 1)
             i_repeat_unit = parall(repeat_unit_imolty,i_repeat_unit)
             repeat_unit_rx(1:repeat_unit_nunit) = rxu(i_repeat_unit, 1:repeat_unit_nunit)
             repeat_unit_ry(1:repeat_unit_nunit) = ryu(i_repeat_unit, 1:repeat_unit_nunit)
             repeat_unit_rz(1:repeat_unit_nunit) = rzu(i_repeat_unit, 1:repeat_unit_nunit)
         else
             ! use the existing conformation
             do unit_index = 1, repeat_unit_nunit
                 ibead = repeat_unit_list(unit_index)
                 repeat_unit_rx(unit_index) = rxu(i, ibead)
                 repeat_unit_ry(unit_index) = ryu(i, ibead)
                 repeat_unit_rz(unit_index) = rzu(i, ibead)
             end do

             ! energy calculation step1: temporarily store the bead positions of repeat unit 1
             i_repeat_unit = parall(repeat_unit_imolty,1)
             rxnew(1:repeat_unit_nunit) = rxu(i_repeat_unit, 1:repeat_unit_nunit)
             rynew(1:repeat_unit_nunit) = ryu(i_repeat_unit, 1:repeat_unit_nunit)
             rznew(1:repeat_unit_nunit) = rzu(i_repeat_unit, 1:repeat_unit_nunit)

             ! step2: use the repeat unit coordinates from original molecule to replace repeat unit 1
             do unit_index = 1, repeat_unit_nunit
                 ibead = repeat_unit_list(unit_index)
                 rxu(i_repeat_unit, unit_index) = rxu(i, ibead)
                 ryu(i_repeat_unit, unit_index) = ryu(i, ibead)
                 rzu(i_repeat_unit, unit_index) = rzu(i, ibead)
             end do
         end if

         ! compute the bonded energy for this molecule
         do ibead = 1, repeat_unit_nunit
             rxuion(ibead, 1) = rxu(i_repeat_unit, ibead)
             ryuion(ibead, 1) = ryu(i_repeat_unit, ibead)
             rzuion(ibead, 1) = rzu(i_repeat_unit, ibead)
             qquion(ibead, 1) = qqu(i_repeat_unit, ibead)
         end do

         call U_bonded(i_repeat_unit,repeat_unit_imolty,vvib,vbend,vtg)

         call energy(i_repeat_unit,repeat_unit_imolty,v_repeat_unit,1,gcbmc_box_num,1,&
             repeat_unit_nunit,.true.,ovrlap,.false.,.true.,.false.,.false.)

         repeat_unit_energy = v_repeat_unit
         repeat_unit_energy(ivTot) = repeat_unit_energy(ivTot) + vvib + vbend + vtg
         repeat_unit_energy(ivStretching) = vvib
         repeat_unit_energy(ivBending) = vbend
         repeat_unit_energy(ivTorsion) = vtg

         if (.not. lnew) then
             ! recover the coordinates
             rxu(i_repeat_unit, 1:repeat_unit_nunit) = rxnew(1:repeat_unit_nunit)
             ryu(i_repeat_unit, 1:repeat_unit_nunit) = rynew(1:repeat_unit_nunit)
             rzu(i_repeat_unit, 1:repeat_unit_nunit) = rznew(1:repeat_unit_nunit)

             ! recover r*new
             rxnew(1:repeat_unit_nunit) = rxnew_temp(1:repeat_unit_nunit)
             rynew(1:repeat_unit_nunit) = rynew_temp(1:repeat_unit_nunit)
             rznew(1:repeat_unit_nunit) = rznew_temp(1:repeat_unit_nunit)
         end if

         !---start to grow the original molecule
         ! for later translation and rotation of the repeat unit so that first two beads match ifrom and iprev
         if (lnew) then
             xtarget(1) = rxnew(iuprev)
             ytarget(1) = rynew(iuprev)
             ztarget(1) = rznew(iuprev)
             xtarget(2) = rxnew(iufrom)
             ytarget(2) = rynew(iufrom)
             ztarget(2) = rznew(iufrom)
         else
             xtarget(1) = rxu(i, iuprev)
             ytarget(1) = ryu(i, iuprev)
             ztarget(1) = rzu(i, iuprev)
             xtarget(2) = rxu(i, iufrom)
             ytarget(2) = ryu(i, iufrom)
             ztarget(2) = rzu(i, iufrom)
         end if

         ! for later dihedral_rigrot the repeat unit
         bead1_coord(1) = xtarget(1)
         bead1_coord(2) = ytarget(1)
         bead1_coord(3) = ztarget(1)
         bead2_coord(1) = xtarget(2)
         bead2_coord(2) = ytarget(2)
         bead2_coord(3) = ztarget(2)

         ! find the two overlapping beads for the align_lines and dihedral_rigrot
         do unit_index = 1, repeat_unit_nunit
             ibead = repeat_unit_list(unit_index)
             if (ibead .eq. iufrom) second_bead = unit_index
             if (ibead .eq. iuprev) first_bead = unit_index
         end do

         call align_lines(first_bead, second_bead, repeat_unit_nunit, repeat_unit_rx, &
                    repeat_unit_ry, repeat_unit_rz, xtarget, ytarget, ztarget)

         do ip = 1, ichoi
             bsum_tor(ip) = 0.0E0_dp
             my_itrial = 0

             do itor = my_start, my_end
                 my_itrial = my_itrial + 1

                 if (lnew .or. ip .gt. 1 .or. itor .gt. 1) then
                     ! new conformation
                     phi = random(rid) * twopi
                 else
                     ! old conformation
                     phi = 0.0E0_dp
                 end if

                 call dihedral_rigrot(bead1_coord, bead2_coord, repeat_unit_nunit, repeat_unit_rx, &
                    repeat_unit_ry, repeat_unit_rz, phi, repeat_unit_rxp, repeat_unit_ryp, repeat_unit_rzp)

                 ! update bond vectors for the calculation of the torsional potential
                 do iu = 1, repeat_unit_nunit
                     ibead = repeat_unit_list(iu)
                     lexist(ibead) = .true.
                     do iv = 1, invib(repeat_unit_imolty, iu)
                         ju = ijvib(repeat_unit_imolty, iu, iv)
                         jbead = repeat_unit_list(ju)
                         xvec(ibead, jbead) = repeat_unit_rxp(ju) - repeat_unit_rxp(iu)
                         yvec(ibead, jbead) = repeat_unit_ryp(ju) - repeat_unit_ryp(iu)
                         zvec(ibead, jbead) = repeat_unit_rzp(ju) - repeat_unit_rzp(iu)
                     end do
                 end do

                 ! find the new torsional angle and compute the torsional potential
                 vnew_torsion = 0.0E0_dp

                 ! the middle two beads are iufrom and iuprev
                 do iv = 1, invib(imolty, iuprev)
                     iu = ijvib(imolty, iuprev, iv)
                     if (iu .ne. iufrom) then
                         do it = 1, intor(imolty, iu)
                             jut2 = ijtor2(imolty, iu, it)
                             jut3 = ijtor3(imolty, iu, it)
                             if (jut2 .eq. iuprev .and. jut3 .eq. iufrom) then
                                 jut4 = ijtor4(imolty, iu, it)

                                 ! if jut4 does NOT exist, do NOT count it
                                 if (.not. lexist(jut4)) then
                                     cycle
                                 else
                                    vnew_torsion = vnew_torsion + vtorso(xvec(jut4,jut3),yvec(jut4,jut3),&
                                       zvec(jut4,jut3),xvec(jut3,jut2),yvec(jut3,jut2),zvec(jut3,jut2),&
                                       xvec(jut2,iu),yvec(jut2,iu),zvec(jut2,iu),ittor(imolty,iu,it))
                                 end if
                             end if
                         end do
                     end if
                 end do

                 ! compute boltzmann factor
                 my_bf_tor(my_itrial) = exp (-vnew_torsion * beta)

                 ! store vtorsion and phi for this trial
                 my_phitors(my_itrial) = phi
                 my_vtorsion(my_itrial) = vnew_torsion
                 bsum_tor(ip) = bsum_tor(ip) + my_bf_tor(my_itrial)

                 ! recover lexist
                 lexist = lexist_temp
             end do

             ! select one of the trial phi
             if (lnew .or. ip .gt. 1) then
                 rand_tor = random(rid) * bsum_tor(ip)
                 accum_prob = 0.0E0_dp

                 do itor = 1, rcounts(myid+1)
                     accum_prob = accum_prob + my_bf_tor(itor)
                     if (rand_tor .lt. accum_prob) then
                         phi = my_phitors(itor)
                         my_vtorsion_procs(ip) = my_vtorsion(itor)
                         exit
                     end if
                 end do
             else
                 phi = my_phitors(1)
                 my_vtorsion_procs(ip) = my_vtorsion(1)
             end if

             if (numprocs .gt. 1) then
                 call mp_allgather(bsum_tor(ip),bf_tor,groupid)
                 bsum_tor(ip)=sum(bf_tor(1:numprocs))
                 call mp_allgather(my_vtorsion_procs(ip),vtorsion,groupid)
                 call mp_allgather(phi,phitors,groupid)

                 if (lnew .or. ip .gt. 1) then
                     do itor = 1, numprocs
                         rand_tor = random(-1) * bsum_tor(ip)
                         accum_prob = 0.0E0_dp
                         if (rand_tor .lt. accum_prob) then
                             phi = phitors(itor)
                             my_vtorsion_procs(ip) = vtorsion(itor)
                             exit
                         end if
                     end do
                 else
                     phi = phitors(1)
                     my_vtorsion_procs(ip) = vtorsion(1)
                 end if
             end if

             bsum_tor(ip) = bsum_tor(ip) / real(ichtor,dp)

             ! update the bead coordinates and existence
             if (lnew .or. ip .gt. 1) then
                 call dihedral_rigrot(bead1_coord, bead2_coord, repeat_unit_nunit, repeat_unit_rx, &
                    repeat_unit_ry, repeat_unit_rz, phi, repeat_unit_rxp, repeat_unit_ryp, repeat_unit_rzp)
             else
                 do unit_index = 1, repeat_unit_nunit
                     ibead = repeat_unit_list(unit_index)
                     repeat_unit_rxp(unit_index) = rxu(i, ibead)
                     repeat_unit_ryp(unit_index) = ryu(i, ibead)
                     repeat_unit_rzp(unit_index) = rzu(i, ibead)
                 end do
             end if

             bead_count = 0
             do unit_index = 1, repeat_unit_nunit
                 ibead = repeat_unit_list(unit_index)
                 if (ibead .ne. iufrom .and. ibead .ne. iuprev) then !as long as they are not overlapping beads
                     bead_count = bead_count + 1
                     rxp(bead_count, ip) = repeat_unit_rxp(unit_index)
                     ryp(bead_count, ip) = repeat_unit_ryp(unit_index)
                     rzp(bead_count, ip) = repeat_unit_rzp(unit_index)
                 end if
             end do
         end do

         ! now that we have phi, update the max len for boltz
         maxlen = 0.0E0_dp
         bead_count = 0
         do unit_index = 1, repeat_unit_nunit
             ibead = repeat_unit_list(unit_index)
             if (ibead .ne. iufrom .and. ibead .ne. iuprev) then
                 bead_count = bead_count + 1
                 glist(bead_count) = ibead

                 do ip = 1, ichoi
                     length = sqrt((rxp(bead_count, ip) - rxu(i, iufrom))**2 + &
                        (ryp(bead_count, ip) - ryu(i, iufrom))**2 + &
                        (rzp(bead_count, ip) - rzu(i, iufrom))**2)
                     if (length .gt. maxlen) maxlen = length
                 end do
             end if
         end do

         ! now that we have the trial site, compute non-bonded energy
         call boltz(lnew,.false.,ovrlap,i,i,imolty,ibox,ichoi,iufrom,repeat_unit_nunit-2,glist,maxlen)

         if (ovrlap) then
             lterm = .true.
             return
         end if

         ! now that we have performed all the growths, it's time to sum up and calculate the weight
         bsum = 0.0E0_dp
         do ip = 1, ichoi   !< the weight includes torsional and LJ/qq
             bsum = bsum + bfac(ip) * bsum_tor(ip)
         end do

         if (lnew) then
             weight = weight * bsum
             if ( weight .lt. softlog ) then
                 lterm=.true.
                 return
             end if

             ! select one trial site as the final one
             accum_prob = 0.0E0_dp
             rand_nonb = random(-1) * bsum

             do ip = 1, ichoi
                 if ( .not. lovr(ip) ) then
                     accum_prob = accum_prob + bfac(ip) * bsum_tor(ip)
                     if (rand_nonb .lt. accum_prob) then
                         exit
                     end if
                 end if
             end do
         else
             weiold = weiold * bsum
             if (weiold .lt. softlog) write(io_output,*)  '###old weight too low in group-CBMC growth'
         end if

         ! update the trial energy and weight
         if (lnew) then
             vnew(ivTot) = vnew(ivTot) + repeat_unit_energy(ivTot) + my_vtorsion_procs(ip) + vtr(ivTot,ip)
             vnew(ivStretching) = vnew(ivStretching) + repeat_unit_energy(ivStretching)
             vnew(ivBending) = vnew(ivBending) + repeat_unit_energy(ivBending)
             vnew(ivTorsion) = vnew(ivTorsion) + repeat_unit_energy(ivTorsion) + my_vtorsion_procs(ip)
             vnew(ivExt)   = vnew(ivExt)   + repeat_unit_energy(ivExt) + vtr(ivExt,ip)
             vnew(ivIntraLJ) = vnew(ivIntraLJ) + repeat_unit_energy(ivIntraLJ) + vtr(ivIntraLJ,ip)
             vnew(ivInterLJ) = vnew(ivInterLJ) + vtr(ivInterLJ,ip)
             vnew(ivElect) = vnew(ivElect) + repeat_unit_energy(ivElect) + vtr(ivElect,ip)
             vnew(ivEwald) = vnew(ivEwald) + repeat_unit_energy(ivEwald) + vtr(ivEwald,ip)
         else
             vold(ivTot) = vold(ivTot) + repeat_unit_energy(ivTot) + my_vtorsion_procs(1) + vtr(ivTot,1)
             vold(ivStretching) = vold(ivStretching) + repeat_unit_energy(ivStretching)
             vold(ivBending) = vold(ivBending) + repeat_unit_energy(ivBending)
             vold(ivTorsion) = vold(ivTorsion) + repeat_unit_energy(ivTorsion) + my_vtorsion_procs(1)
             vold(ivExt)   = vold(ivExt)  + repeat_unit_energy(ivExt) + vtr(ivExt,1)
             vold(ivIntraLJ) = vold(ivIntraLJ) + repeat_unit_energy(ivIntraLJ) + vtr(ivIntraLJ,1)
             vold(ivInterLJ) = vold(ivInterLJ) + vtr(ivInterLJ,1)
             vold(ivElect) = vold(ivElect) + repeat_unit_energy(ivElect) + vtr(ivElect,1)
             vold(ivEwald) = vold(ivEwald) + repeat_unit_energy(ivEwald) + vtr(ivEwald,1)
         end if

         ! update the coordinates and vectors
         bead_count = 0
         do unit_index = 1, repeat_unit_nunit
             ibead = repeat_unit_list(unit_index)
             if (ibead .ne. iufrom .and. ibead .ne. iuprev) then
                 bead_count = bead_count + 1
                 lexist(ibead) = .true.
                 if (lnew) then
                     rxnew(ibead) = rxp(bead_count, ip)
                     rynew(ibead) = ryp(bead_count, ip)
                     rznew(ibead) = rzp(bead_count, ip)
                 end if

                 do iv = 1, invib(imolty, ibead)
                     ju = ijvib(imolty, ibead, iv)
                     if (lexist(ju)) then
                         if (lnew) then
                             xvec(ibead, ju) = rxnew(ju) - rxnew(ibead)
                             yvec(ibead, ju) = rynew(ju) - rynew(ibead)
                             zvec(ibead, ju) = rznew(ju) - rznew(ibead)
                         else
                             xvec(ibead, ju) = rxu(i, ju) -rxu(i, ibead)
                             yvec(ibead, ju) = ryu(i, ju) -ryu(i, ibead)
                             zvec(ibead, ju) = rzu(i, ju) -rzu(i, ibead)
                         end if
                         distij(ibead, ju) = sqrt(xvec(ibead, ju)**2 + yvec(ibead, ju)**2 + zvec(ibead, ju)**2)
                         xvec(ju, ibead) = - xvec(ibead, ju)
                         yvec(ju, ibead) = - yvec(ibead, ju)
                         zvec(ju, ibead) = - zvec(ibead, ju)
                         distij(ju, ibead) = distij(ibead, ju)
                     end if
                 end do
             end if
         end do

         ! recover some variables
         rxnew_temp(1:igrow) = rxnew(1:igrow)
         rynew_temp(1:igrow) = rynew(1:igrow)
         rznew_temp(1:igrow) = rznew(1:igrow)
         lexist_temp(1:igrow) = lexist(1:igrow)

         if (l_reach_end) then
            lexshed(1:igrow) = lexshed_original(1:igrow)
            exit
         else
            ! now that we have grown this segment, move along
            iw = iw + 1
         end if

     end do
  end subroutine group_cbmc_grow

  subroutine allocate_cbmc()
    integer::jerr
    if (allocated(vtr)) deallocate(vtr,vtrorient,vtrelect_intra,vtrelect_inter,bfac,rxp,ryp,rzp,vwellipswot,vwellipswnt,vipswot,vipswnt,lovr,vtvib,vtgtr,vtbend,bsum_tor,stat=jerr)
    allocate(vtr(nEnergy,nchmax),vtrorient(nchmax)&
     ,vtrelect_intra(nchmax),vtrelect_inter(nchmax),bfac(nchmax)&
     ,rxp(numax,nchmax),ryp(numax,nchmax),rzp(numax,nchmax)&
     ,vwellipswot(nchmax),vwellipswnt(nchmax),vipswot(nchmax)&
     ,vipswnt(nchmax),lovr(nchmax),vtvib(nchmax),vtgtr(nchmax)&
     ,vtbend(nchmax),bsum_tor(nchmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_cbmc: allocation failed',jerr)
    end if
  end subroutine allocate_cbmc

  subroutine init_cbmc(io_input,lprint)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::jerr,i
    integer,allocatable::avbmc_version(:)
    namelist /mc_cbmc/ rcutin,pmcb,pmcbmt,pmall,nchoi1,nchoi,nchoir,nchoih,nchoig,nchtor,nchbna,nchbnb,icbdir,icbsta&
     ,rbsmax,rbsmin,avbmc_version,first_bead_to_swap,pmbias,pmbsmt,pmbias2&
     ,pmfix,pmgroup,lrig,lpresim,iupdatefix

    if (allocated(lexshed)) deallocate(lexshed,llplace,lpnow,lsave,first_bead_to_swap,bncb,bscb,fbncb,fbscb,iend,ipast,pastnum,fclose,fcount,iwbef,ibef,befnum,xx,yy,zz,distij,nextnum,inext,kforceb,equilb,flength,vequil,vkforce,rlist,rfrom,rprev,rnum,iplace,pfrom,pnum,pprev,avbmc_version,stat=jerr)
    allocate(lexshed(numax),llplace(ntmax),lpnow(numax),lsave(numax),bncb(ntmax,nbxmax,numax),bscb(ntmax,2,nbxmax,numax),fbncb(ntmax,nbxmax,numax)&
     ,fbscb(ntmax,2,nbxmax,numax),iend(numax),ipast(numax,numax),pastnum(numax),fclose(numax,numax),fcount(numax),iwbef(numax)&
     ,ibef(numax,numax),befnum(numax),xx(numax),yy(numax),zz(numax),distij(numax,numax),nextnum(numax),inext(numax,numax)&
     ,kforceb(numax,numax),equilb(numax,numax),flength(numax,numax),vequil(numax,numax),vkforce(numax,numax),rlist(numax,numax)&
     ,rfrom(numax),rprev(numax),rnum(numax),iplace(numax,numax),pfrom(numax),pnum(numax),pprev(numax),avbmc_version(nmolty)&
     ,first_bead_to_swap(nmolty),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_cbmc: allocation failed',jerr)

    llplace=.FALSE.
    lsave=.FALSE.
    bncb = 0.0E0_dp
    bscb = 0.0E0_dp
    fbncb = 0.0E0_dp
    fbscb = 0.0E0_dp
    counthist=0

    ! defaults for namelist mc_cbmc
    do i=1,nmolty
       pmcbmt(i)=real(i,dp)/nmolty
    end do
    pmall=0.0_dp
    nchoi1=32
    nchoi=16
    nchoir=16
    nchoih=1
    nchoig=16
    nchtor=100
    nchbna=1000
    nchbnb=1000
    icbdir=0
    icbsta=0
    first_bead_to_swap=1
    avbmc_version=0
    pmbias=0.0_dp
    pmbsmt=0.0_dp
    pmbias2=0.0_dp
    pmfix=0.0_dp
    pmgroup=0.0_dp
    lrig=.false.

    !> read namelist mc_cbmc
    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=mc_cbmc,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_cbmc',jerr)
    end if

    call mp_bcast(rcutin,1,rootid,groupid)
    call mp_bcast(pmcb,1,rootid,groupid)
    call mp_bcast(pmcbmt,nmolty,rootid,groupid)
    call mp_bcast(pmall,nmolty,rootid,groupid)
    call mp_bcast(nchoi1,nmolty,rootid,groupid)
    call mp_bcast(nchoi,nmolty,rootid,groupid)
    call mp_bcast(nchoir,nmolty,rootid,groupid)
    call mp_bcast(nchoih,nmolty,rootid,groupid)
    call mp_bcast(nchoig,nmolty,rootid,groupid)
    call mp_bcast(nchtor,nmolty,rootid,groupid)
    call mp_bcast(nchbna,nmolty,rootid,groupid)
    call mp_bcast(nchbnb,nmolty,rootid,groupid)
    call mp_bcast(icbdir,nmolty,rootid,groupid)
    call mp_bcast(icbsta,nmolty,rootid,groupid)
    call mp_bcast(rbsmax,1,rootid,groupid)
    call mp_bcast(rbsmin,1,rootid,groupid)
    call mp_bcast(first_bead_to_swap,nmolty,rootid,groupid)
    call mp_bcast(avbmc_version,nmolty,rootid,groupid)
    call mp_bcast(pmbias,nmolty,rootid,groupid)
    call mp_bcast(pmbsmt,nmolty,rootid,groupid)
    call mp_bcast(pmbias2,nmolty,rootid,groupid)
    call mp_bcast(pmfix,nmolty,rootid,groupid)
    call mp_bcast(pmgroup,nmolty,rootid,groupid)
    call mp_bcast(lrig,nmolty,rootid,groupid)
    call mp_bcast(lpresim,1,rootid,groupid)
    call mp_bcast(iupdatefix,1,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_CBMC','------------------------------------------'
       write(io_output,'(A,F6.3,A)') 'CBMC inner cutoff (rcutin): ',rcutin,' [Ang]'
       write(io_output,'(2(A,F6.3),A)') 'AVBMC outer cutoff (rbsmax): ',rbsmax,' [Ang], inner cutoff (rbsmin): ',rbsmin,' [Ang]'
       write(io_output,'(A,L2)') 'lpresim: ',lpresim
       write(io_output,'(A,I0)') 'iupdatefix: ',iupdatefix
       write(io_output,'(A,G16.9)') 'pmcb: ',pmcb
       write(io_output,'(/,A)') 'molecule type: nchoi1  nchoi nchoir nchoih nchtor nchbna nchbnb icbdir icbsta first_bead_to_swap'
       do i=1,nmolty
          write(io_output,'(I13,A,10(1X,I6))') i,':',nchoi1(i),nchoi(i),nchoir(i),nchoih(i),nchtor(i),nchbna(i),nchbnb(i)&
           ,icbdir(i),icbsta(i),first_bead_to_swap(i)
       end do
       write(io_output,'(/,A)')&
        'molecule type:    pmcbmt         pmall  avbmc_version    pmbias        pmbsmt       pmbias2         pmfix   lrig'
       do i=1,nmolty
          write(io_output,'(I13,A,2(1X,G13.6),1X,I10,4(1X,G13.6),1X,L2)') i,':',pmcbmt(i),pmall(i),avbmc_version(i),pmbias(i)&
           ,pmbsmt(i),pmbias(2),pmfix(i),lrig(i)
       end do
    end if

    vol_eff = (4.0E0_dp/3.0E0_dp)*onepi*(rbsmax*rbsmax*rbsmax-rbsmin*rbsmin*rbsmin)
    nchmax=0
    nchtor_max=0
    nchbn_max=0
    lavbmc1=.false.
    lavbmc2=.false.
    lavbmc3=.false.
    lbias=.false.
    lneighbor = .false.
    do i=1,nmolty
       if (nchoi1(i).gt.nchmax) then
          nchmax=nchoi1(i)
       end if
       if (nchoi(i).gt.nchmax) then
          nchmax=nchoi(i)
       end if
       if (nchoir(i).gt.nchmax) then
          nchmax=nchoir(i)
       end if
       if (nchoih(i).ne.1.and.nunit(i).eq.nugrow(i)) call err_exit(__FILE__,__LINE__,'nchoih must be 1 (one) if nunit = nugrow'&
        ,myid+1)
       if (nchtor(i).gt.nchtor_max) then
          nchtor_max=nchtor(i)
       end if
       if (nchbna(i).gt.nchbn_max) then
          nchbn_max=nchbna(i)
       end if
       if (nchbnb(i).gt.nchbn_max) then
          nchbn_max=nchbnb(i)
       end if
       if (abs(icbsta(i)).gt.nunit(i)) call err_exit(__FILE__,__LINE__,'icbsta > nunit for molecule '//integer_to_string(i),myid+1)

       if (avbmc_version(i).eq.1) then
          lavbmc1(i)=.true.
       else if (avbmc_version(i).eq.2) then
          lavbmc2(i)=.true.
       else if (avbmc_version(i).eq.3) then
          lavbmc3(i)=.true.
       end if
       if (lavbmc1(i).or.lavbmc2(i).or.lavbmc3(i)) then
          lbias(i) = .true.
       end if
       if ((lavbmc2(i).or.lavbmc3(i)).and.(.not.lgaro)) lneighbor = .true.

       if (lring(i).and.pmfix(i).lt.1.and..not.lrig(i)) call err_exit(__FILE__,__LINE__,'a ring can only be used with safe-cbmc'&
        ,myid+1)
    end do

    if (any(lbias(1:nmolty)).and.rbsmax.lt.rbsmin) call err_exit(__FILE__,__LINE__,'rbsmax should be greater than rbsmin',myid+1)

    call allocate_cbmc()
    call read_safecbmc(io_input,lprint)
    if (ANY(pmgroup(1:nmolty).gt.0)) call read_groupcbmc(io_input,lprint)
  end subroutine init_cbmc

  subroutine read_groupcbmc(io_input,lprint)
    use var_type,only:default_path_length,default_string_length
    use util_string,only:uppercase
    use util_files,only:get_iounit,readLine
    use util_mp,only:mp_bcast
    integer,intent(in)::io_input
    logical,intent(in)::lprint
    character(LEN=default_string_length)::line_in
    integer::imol,unit_index,jerr,gcbmc_molty,repeat_unit_nunit

    ! Look for section GROUP_CBMC
    if (myid.eq.rootid .AND.ANY(pmgroup(1:nmolty).gt.0) ) then
        REWIND(io_input)
        CYCLE_READ_GROUPCBMC:DO
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section GROUP_CBMC not found',jerr)

            if (UPPERCASE(line_in(1:10)).eq.'GROUP_CBMC') then
                exit cycle_read_groupcbmc
            end if
        END DO CYCLE_READ_GROUPCBMC

        call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
        if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section GROUP_CBMC',jerr)
        if (UPPERCASE(line_in(1:14)).eq.'END GROUP_CBMC') &
            call err_exit(__FILE__,__LINE__,'Section GROUP_CBMC not complete!',myid+1)
        read(line_in,*) gcbmc_box_num

        call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
        if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section GROUP_CBMC',jerr)
        read(line_in,*) gcbmc_mol_num

        allocate(gcbmc_mol_list(nmolty))
        allocate(gcbmc_unit_num(gcbmc_mol_num))
        allocate(gcbmc_unit_moltype(gcbmc_mol_num, maxval(nunit)))
        allocate(gcbmc_unit_list(gcbmc_mol_num, maxval(nunit), maxval(nunit)))

        do imol=1, gcbmc_mol_num
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            read(line_in,*) gcbmc_molty, gcbmc_unit_num(imol)
            gcbmc_mol_list(gcbmc_molty) = imol

            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            read(line_in,*) gcbmc_unit_moltype(imol,1:gcbmc_unit_num(imol))

            ! corresponding list for each segment
            do unit_index = 1, gcbmc_unit_num(imol)
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                repeat_unit_nunit = nunit(gcbmc_unit_moltype(imol, unit_index))
                read(line_in,*) gcbmc_unit_list(imol, unit_index, 1:repeat_unit_nunit)
            end do
        end do
    end if

    call mp_bcast(gcbmc_box_num,1,rootid,groupid)
    call mp_bcast(gcbmc_mol_num,1,rootid,groupid)
    call mp_bcast(gcbmc_mol_list,nmolty,rootid,groupid)
    call mp_bcast(gcbmc_unit_num,gcbmc_mol_num,rootid,groupid)
    call mp_bcast(gcbmc_unit_moltype,gcbmc_mol_num*maxval(nunit),rootid,groupid)
    call mp_bcast(gcbmc_unit_list,gcbmc_mol_num*maxval(nunit)*maxval(nunit),rootid,groupid)

  end subroutine read_groupcbmc

  subroutine read_safecbmc(io_input,lprint)
    use var_type,only:default_path_length,default_string_length
    use util_string,only:uppercase
    use util_files,only:get_iounit,readLine
    use util_mp,only:mp_bcast
    integer,intent(in)::io_input
    LOGICAL,INTENT(IN)::lprint
    character(LEN=default_path_length),parameter::file_safecbmc='fort.23'
    character(LEN=default_string_length)::line_in
    integer::io_safecbmc,jerr,imol,i,j,bin,bdum

    ! Looking for section SAFE_CBMC
    if (myid.eq.rootid.and.ANY(lrig(1:nmolty))) then
       REWIND(io_input)
       CYCLE_READ_SAFECBMC:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section SAFE_CBMC not found',jerr)

          if (UPPERCASE(line_in(1:9)).eq.'SAFE_CBMC') then
             exit cycle_read_safecbmc
          end if
       END DO CYCLE_READ_SAFECBMC
    end if

    do imol=1,nmolty+1
       if (imol.ne.nmolty+1) then
          if (.not.lrig(imol)) cycle
       else if (ALL(.not.lrig(1:nmolty))) then
          exit
       end if

       if (myid.eq.rootid) then
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SAFE_CBMC',jerr)

          if (UPPERCASE(line_in(1:13)).eq.'END SAFE_CBMC') then
             if (imol.ne.nmolty+1) call err_exit(__FILE__,__LINE__,'Section SAFE_CBMC not complete!',myid+1)
             exit
          else if (imol.eq.nmolty+1) then
             call err_exit(__FILE__,__LINE__,'Section SAFE_CBMC has more than nmolty records!',myid+1)
          end if

          read(line_in,*) nrig(imol)
          if (lprint) then
             write(io_output,'(2(A,I0))') '   Molecule type ',imol,': nrig = ',nrig(imol)
          end if
       else if (imol.eq.nmolty+1) then
          exit
       end if

       call mp_bcast(nrig(imol),1,rootid,groupid)

       if (nrig(imol).gt.0) then
          ! read in specific points to keep rigid in growth
          if (myid.eq.rootid) then
             do i = 1, nrig(imol)
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SAFE_CBMC',jerr)
                read(line_in,*) irig(imol,i),frig(imol,i)

                if (lprint) then
                   write(io_output,'(3(A,I0))') '      rigid part ',i,': irig = ',irig(imol,i),', frig = ',frig(imol,i)
                end if
             end do
          end if

          call mp_bcast(irig(imol,1:nrig(imol)),nrig(imol),rootid,groupid)
          call mp_bcast(frig(imol,1:nrig(imol)),nrig(imol),rootid,groupid)

          lrigi(imol,irig(imol,:)) = .true.
       else
          ! we will pick irig at random in each case if nrig = 0
          if (myid.eq.rootid) then
             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SAFE_CBMC',jerr)
             read(line_in,*) nrigmin(imol),nrigmax(imol)
          end if

          call mp_bcast(nrigmin(imol),1,rootid,groupid)
          call mp_bcast(nrigmax(imol),1,rootid,groupid)

          if (lprint) then
             write(io_output,'(2(A,I0))') '      nrigmin: ',nrigmin(imol),', nrigmax: ',nrigmax(imol)
          end if
       end if
    end do
! -------------------------------------------------------------------

    if (myid.eq.rootid.AND.ANY(pmfix(1:nmolty).gt.0)) then
       io_safecbmc=get_iounit()
       open(unit=io_safecbmc,access='sequential',action='read',file=file_safecbmc,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open safecbmc file '//trim(file_safecbmc),jerr)

       do imol = 1, nmolty
          if (pmfix(imol).gt.0) then !> \bug Allow only one molecule type to have pmfix > 0
             read(io_safecbmc,*) counttot
             ! read in from fort.23 the bead-bead distribution
             do i = 1, iring(imol)
                do j = 1, iring(imol)
                   if (i.eq.j) cycle
                   do bin = 1, maxbin
                      read(io_safecbmc,*) bdum,probf(i,j,bin)
                   end do
                end do
             end do
             exit
          end if
       end do
    end if

    call mp_bcast(counttot,1,rootid,groupid)
    call mp_bcast(probf,60*60*maxbin,rootid,groupid) ! Paul
    hist = 0._dp

    if (myid.eq.rootid.AND.ANY(pmfix(1:nmolty).gt.0)) close(io_safecbmc)
  end subroutine read_safecbmc

  subroutine opt_safecbmc()
    integer::imolty,j,k,bin
    real::histrat,histtot

    do imolty = 1, nmolty
       if (pmfix(imolty).gt.0.0001_dp) then
          ! readjust fixed end point data to assure optimum efficiency ***
          counttot = counttot + counthist
          histrat = real(counthist,dp) / real(counttot,dp)

          ! reset counthist
          counthist = 0
          do j = 1, iring(imolty)
             do k = 1, iring(imolty)
                if (j.eq.k) cycle
                histtot = 0
                do bin = 1, maxbin
                   histtot = histtot + hist(j,k,bin)
                end do

                if (histtot.eq.0) cycle

                ! normalize and multiply hist by its weighting using above condition
                do bin = 1, maxbin
                   hist(j,k,bin) = hist(j,k,bin)*histrat/histtot
                   ! add weighed hist to hist
                   probf(j,k,bin) = probf(j,k,bin) + hist(j,k,bin)
                   ! reset hist to zero for next iteration
                   hist(j,k,bin) = 0
                end do
                ! renormalize new distribution
                histtot = 0
                do bin = 1, maxbin
                   histtot = histtot + probf(j,k,bin)
                end do
                do bin = 1, maxbin
                   probf(j,k,bin) = probf(j,k,bin) / histtot
                end do
             end do
          end do
       end if
    end do
  end subroutine opt_safecbmc

!> \brief write some information about config performance
  subroutine output_cbmc_stats(io_output)
    integer,intent(in)::io_output
    integer::i,inb,ii,ibox
    real::pscb1,pscb2

    write(io_output,*)
    write(io_output,*) '### Configurational-bias ###'
    write(io_output,*)
    do i = 1, nmolty
       do ibox = 1, nbox
          write(io_output,"(A,I0,A)",advance='no') 'molecule typ = ',i,'    '
          write(io_output,'(A10,A,I0)',advance='no')molecname(i),' in box ',ibox
          write(io_output,*)
          write(io_output,*) '    length  attempts  succ.growth  accepted' ,'   %su.gr.    %accep.'
          do inb = 1, nunit(i)
             if ( bncb(i,ibox,inb) .gt. 0.0E0_dp ) then
                pscb1 = bscb(i,1,ibox,inb) * 100.0E0_dp / bncb(i,ibox,inb)
                pscb2 = bscb(i,2,ibox,inb) * 100.0E0_dp / bncb(i,ibox,inb)
                write(io_output,'(i9,3f10.1,2f10.2)') inb, bncb(i,ibox,inb), bscb(i,1,ibox,inb), bscb(i,2,ibox,inb), pscb1, pscb2
             end if
          end do
          if (pmfix(i).gt.0.0E0_dp) then
             write(io_output,*) ' SAFE-CBMC move '
             write(io_output,*) '    length  attempts  succ.growth  ', 'accepted   %su.gr.    %accep.'
             do inb = 1, nunit(i)
                if (fbncb(i,ibox,inb) .gt. 0.0E0_dp ) then
                   pscb1 = fbscb(i,1,ibox,inb) * 100.0E0_dp  / fbncb(i,ibox,inb)
                   pscb2 = fbscb(i,2,ibox,inb) * 100.0E0_dp  / fbncb(i,ibox,inb)
                   write(io_output,'(i9,3f10.1,2f10.2)') inb,fbncb(i,ibox,inb), fbscb(i,1,ibox,inb), fbscb(i,2,ibox,inb) , pscb1, pscb2
                end if
             end do
          end if
       end do
    end do
    write(io_output,*)
  end subroutine output_cbmc_stats

  subroutine output_safecbmc()
    use var_type,only:default_path_length
    use util_files,only:get_iounit
    character(LEN=default_path_length),parameter::file_safecbmc='fort.21'
    integer::io_safecbmc,jerr,j,k,bin,imolty
    real::histtot

    if (lpresim.or.ANY(pmfix(1:nmolty).gt.0)) then
       io_safecbmc=get_iounit()
       open(unit=io_safecbmc,access='sequential',action='write',file=file_safecbmc,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open safecbmc file '//trim(file_safecbmc),jerr)
    end if

    ! normalize and write out presim results in fort.22 **
    if (lpresim) then
       if (counttot.eq.0) then
          write(io_safecbmc,*) counthist
       else
          write(io_safecbmc,*) counttot
       end if

       do j = 1, iring(1)
          do k = 1, iring(1)
             if (j.eq.k) cycle
             histtot = 0
             do bin = 1, maxbin
                hist(j,k,bin) = hist(j,k,bin) + 1.0E0_dp
                histtot = histtot + hist(j,k,bin)
             end do

             do bin = 1, maxbin
                hist(j,k,bin) = hist(j,k,bin) / histtot
                write(io_safecbmc,*) bin,hist(j,k,bin)
             end do
          end do
       end do
    end if

    ! put new distribution back into a file
    do imolty = 1, nmolty
       if (pmfix(imolty).gt.0) then
          if (counttot.eq.0) then
             write(io_safecbmc,*) counthist
          else
             write(io_safecbmc,*) counttot
          end if
          do j = 1, iring(imolty)
             do k = 1, iring(imolty)
                if (j.eq.k) cycle
                do bin = 1, maxbin
                   write(io_safecbmc,*) bin,probf(j,k,bin)
                end do
             end do
          end do
       end if
    end do

    if (lpresim.or.ANY(pmfix(1:nmolty).gt.0)) close(io_safecbmc)
  end subroutine output_safecbmc

  subroutine read_checkpoint_cbmc(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) bncb,bscb,fbncb,fbscb
    call mp_bcast(bncb,ntmax*nbxmax*numax,rootid,groupid)
    call mp_bcast(bscb,ntmax*2*nbxmax*numax,rootid,groupid)
    call mp_bcast(fbncb,ntmax*nbxmax*numax,rootid,groupid)
    call mp_bcast(fbscb,ntmax*2*nbxmax*numax,rootid,groupid)
  end subroutine read_checkpoint_cbmc

  subroutine write_checkpoint_cbmc(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) bncb,bscb,fbncb,fbscb
  end subroutine write_checkpoint_cbmc
end MODULE moves_cbmc
