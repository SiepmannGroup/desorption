module transfer_swatch
  use util_runtime,only:err_exit
  use util_random,only:random
  use sim_system
  use sim_cell
  use energy_kspace,only:recip
  use energy_pairwise,only:energy,coru
  use moves_cbmc,only:rosenbluth,schedule,safeschedule,explct,lexshed,align_lines
  implicit none
  private
  save
  public::swatch,init_swatch,output_swatch_stats,read_checkpoint_swatch,write_checkpoint_swatch

  integer,allocatable::bnswat(:,:,:),bnswat_empty(:,:,:),bsswat(:,:,:) !< accumulators for swatch performance
contains
!> Added intrabox move for two particles within one box
!> in combined move that shares the same parameters.
!> Will also accept rigid (lrigid) molecules.  Contains
!> several critical bug fixes as well.
!> \since 9-25-02 JMS
  subroutine swatch()
    use sim_particle,only:ctrmas,update_coord_in_tree
    use energy_intramolecular,only:U_bonded
    use transfer_shared,only:lopt_bias,update_bias,gcmc_setup,gcmc_cleanup,gcmc_exchange

    logical::lterm
    real::rpair,tweight,tweiold,vnbox(nEnergy,nbxmax),rx_1(numax),ry_1(numax),rz_1(numax),rxut(4,numax),ryut(4,numax),rzut(4,numax)&
     ,waddold,waddnew,dvol,vola,volb,rho,dinsta,wnlog,wolog,wdlog,wswat,v(nEnergy),vdum,delen,deleo,dicount,vrecipn,vrecipo,vdum2&
     ,cwtorfn,cwtorfo&
     ,total_NBE,total_tor,total_bend,total_vib,vtgn,vbendn,vvibn,Rosenbluth_normalization(ntmax),xtarget(2),ytarget(2),ztarget(2)&
     ,vnewtemp(nEnergy),voldtemp(nEnergy) !< temporarily store vnew and vold between regular CBMC and SAFE-CBMC
    integer::ipair,iparty,type_a,type_b,imolta,imoltb,ipairb,boxa,boxb,imola,imolb,ibox,iboxal,iboxbl,izz,from(2*numax)&
     ,prev(2*numax),orgaia,orgaib,orgbia,orgbib,bdmol_a,bdmol_b,s_type,o_type,new,old,ifirst,iprev,imolty,igrow,islen,iunit,iboxnew&
     ,iboxold,self,other,iunita,iunitb,fromsafe(2*numax),prevsafe(2*numax),indexsafe(2*numax),iindex,first_bead,second_bead
    integer::iii,j
    integer::oldchain,newchain,oldunit,newunit,iins
    integer::ic,icbu,jj,mm,imt,jmt,imolin,imolrm
    integer::icallrose
    real::tmpvvib,tmpvbend,tmpvtg
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef __DEBUG__
    write(io_output,*) 'start SWATCH in ',myid
#endif

    ! randomly select chains to switch between boxa boxb
    ! select a pair type to switch
    if (nswaty.gt.1) then
       rpair = random(-1)
       do ipair=1,nswaty
          if (rpair.lt.pmsatc(ipair)) then
             iparty = ipair
             exit
          end if
       end do
    else
       iparty = 1
    end if

    ! select the molecules from the pair
    imolta = nswatb(iparty,1)
    imoltb = nswatb(iparty,2)
    type_a = 1
    type_b = 2

    ! choose box A and box B at random
    if (nswtcb(iparty).gt.1) then
       rpair = random(-1)
       do ipair = 1, nswtcb(iparty)
          if (rpair.lt.pmswtcb(iparty,ipair)) then
             ipairb = ipair
             exit
          end if
       end do
    else
       ipairb = 1
    end if

    if (random(-1).lt.0.5E0_dp) then
       boxa=box3(iparty,ipairb)
       boxb=box4(iparty,ipairb)
    else
       boxa=box4(iparty,ipairb)
       boxb=box3(iparty,ipairb)
    end if

    ! add one attempt to the count for iparty
    bnswat(iparty,ipairb,boxa) = bnswat(iparty,ipairb,boxa) + 1

    ! check if the particle types are in their boxes
    if (ncmt(boxa,imolta).ne.0) then
       ! get particle from box a
       iboxal = int(real(ncmt(boxa,imolta),dp)*random(-1)) + 1
       imola = parbox(iboxal,boxa,imolta)
       if (moltyp(imola).ne.imolta) call err_exit(__FILE__,__LINE__,'screwup',myid+1)
       if (nboxi(imola).ne.boxa) call err_exit(__FILE__,__LINE__,'problem in swatch',myid+1)
    else if (lgrand.and.boxa.ne.1) then
       call gcmc_setup(imolta,boxa,imola,iboxal)
    end if

    if (ncmt(boxb,imoltb).ne.0) then
       ! get particle from box b
       iboxbl = int(real(ncmt(boxb,imoltb),dp)*random(-1)) + 1
       imolb = parbox(iboxbl,boxb,imoltb)
       if (moltyp(imolb).ne.imoltb) call err_exit(__FILE__,__LINE__,'screwup',myid+1)
       if (nboxi(imolb).ne.boxb) call err_exit(__FILE__,__LINE__,'problem in swatch',myid+1)
    else if (lgrand.and.boxb.ne.1) then
       call gcmc_setup(imoltb,boxb,imolb,iboxbl)
    end if

    if ((ncmt(boxa,imolta).eq.0).or.(ncmt(boxb,imoltb).eq.0)) then
       bnswat_empty(iparty,ipairb,boxa) = bnswat_empty(iparty,ipairb,boxa) + 1
       if (lgrand) then
          if (boxa.ne.1) then
             call gcmc_cleanup(imolta,boxa)
          else if (boxb.ne.1) then
             call gcmc_cleanup(imoltb,boxb)
          end if
       end if
       return
    end if

!************************
! Begin Growth Setups ***
!************************

    ! assign from and prev for each moltyp
    do izz = 1,ncut(iparty,type_a)
       from(type_a+2*(izz-1)) = gswatc(iparty,type_a,1+2*(izz-1))
       prev(type_a+2*(izz-1)) = gswatc(iparty,type_a,2+2*(izz-1))
    end do

    do izz = 1,ncut(iparty,type_b)
       from(type_b+2*(izz-1)) = gswatc(iparty,type_b,1+2*(izz-1))
       prev(type_b+2*(izz-1)) = gswatc(iparty,type_b,2+2*(izz-1))
    end do

    do izz = 1,ncutsafe(iparty,type_a)
       fromsafe(type_a+2*(izz-1)) = gswatcsafe(iparty,type_a,1+3*(izz-1))
       prevsafe(type_a+2*(izz-1)) = gswatcsafe(iparty,type_a,2+3*(izz-1))
       indexsafe(type_a+2*(izz-1)) = gswatcsafe(iparty,type_a,3+3*(izz-1))
    end do

    do izz = 1,ncutsafe(iparty,type_b)
       fromsafe(type_b+2*(izz-1)) = gswatcsafe(iparty,type_b,1+3*(izz-1))
       prevsafe(type_b+2*(izz-1)) = gswatcsafe(iparty,type_b,2+3*(izz-1))
       indexsafe(type_b+2*(izz-1)) = gswatcsafe(iparty,type_b,3+3*(izz-1))
    end do

    ! store number of units in iunita and iunitb
    iunita = nunit(imolta)
    iunitb = nunit(imoltb)

    ! store number of each type in the boxes
    orgaia = ncmt(boxa,imolta)
    orgbia = ncmt(boxa,imoltb)
    orgaib = ncmt(boxb,imolta)
    orgbib = ncmt(boxb,imoltb)

    ! initialize trial weights
    tweight = 1.0E0_dp
    tweiold = 1.0E0_dp

    ! set the trial energies to zero
    vnbox = 0.0E0_dp


    if (boxa.eq.boxb) then
       ! store position 1, a's original site
       do izz = 1, nunit(imolta)
          rx_1(izz) = rxu(imola,izz)
          ry_1(izz) = ryu(imola,izz)
          rz_1(izz) = rzu(imola,izz)
       end do

       ! store same bead coordinates for molecules a and b
       do icbu = 1, nsampos(iparty)
          bdmol_a = splist(iparty,icbu,type_a)
          rxut(3,icbu) = rxu(imola,bdmol_a)
          ryut(3,icbu) = ryu(imola,bdmol_a)
          rzut(3,icbu) = rzu(imola,bdmol_a)

          bdmol_b = splist(iparty,icbu,type_b)
          rxut(4,icbu) = rxu(imolb,bdmol_b)
          ryut(4,icbu) = ryu(imolb,bdmol_b)
          rzut(4,icbu) = rzu(imolb,bdmol_b)
       end do

!**************************
! Liswinc Determination ***
!**************************
       !> \bug this part will not work with rigid swatches
       ! only need this for molecule b, since b is grown after a
       ! but a is grown when some of b doesn't exist
       ! imolty = imoltb
       ! igrow = nugrow(imoltb)

       ! ifirst = from(type_b)
       ! iprev = prev(type_b)

       ! ! determine which beads aren't in the same positions
       ! call schedule(igrow,imolty,islen,ifirst,iprev,3)

       ! if (ncut(iparty,type_b) .gt. 1) then
       !    do izz = 2,ncut(iparty,type_b)
       !       ifirst = from(type_b+2*(izz-1))
       !       iprev = prev(type_b+2*(izz-1))
       !       call schedule(igrow,imolty,islen,ifirst,iprev,5)
       !    end do
       ! end if

       ! ! assign growth schedule for molecule b
       ! do izz = 1,nunit(imolty)
       !    liswinc(izz,imolty) = lexshed(izz)
       !    ! write(io_output,*) izz, liswinc(izz,imolty)
       ! end do
    end if

!***************
! ROSENBLUTH ***
!***************
    ! compute the rosenbluth weights for each molecule type in each box
    oldnew: do ic = 1,2
       if (ic.eq.1) then
          ! if (boxa.eq.boxb) liswatch = .true.
          self = imola
          other = imolb
          s_type = type_a
          o_type = type_b
          iboxnew = boxb
          iboxold = boxa
          imolty = imolta
          igrow = nugrow(imolta)
          iunit = nunit(imolta)
       else
          ! liswatch = .false.
          self = imolb
          other = imola
          s_type = type_b
          o_type = type_a
          iboxnew = boxa
          iboxold = boxb
          imolty = imoltb
          igrow = nugrow(imoltb)
          iunit = nunit(imoltb)
       end if

       ! initialize vnewtemp and voldtemp
       vnewtemp = 0.0E0_dp
       voldtemp = 0.0E0_dp

       ! store the beads that are identical
       do icbu = 1, nsampos(iparty)
          new = splist(iparty,icbu,s_type)
          old = splist(iparty,icbu,o_type)
          rxnew(new) = rxu(other,old)
          rynew(new) = ryu(other,old)
          rznew(new) = rzu(other,old)
       end do

       ! If there exists regular CBMC part, do it; otherwise do SAFE-CBMC part
       if (ncut(iparty,s_type) .gt. 0) then

          ! set up growth schedule
          ifirst = from(s_type)
          iprev = prev(s_type)

          !*** rigid molecule add on ***
          if (lrigid(imolty)) then
             !cc--- JLR 11-24-09
             !cc--- adding some if statements for rigid swaps:
             !cc--- if we have swapped the whole rigid part we will
             !cc--- not compute the vectors from ifirst in old position
             !if (nsampos(iparty).lt.iunit) then
             ! if ((rindex(imolty).eq.0).or.(ifirst.lt.riutry(imolty,1))) then
             if (nsampos(iparty).ge.3) then
                call align_planes(iparty,self,other,s_type,o_type,rxnew,rynew,rznew)
             else
                ! calculate new vector from initial bead
                ! BX: If lrigid=.true., nsampos.eq.1 and nsampos.eq.2 should be treated in another way.
                first_bead = ifirst
                second_bead = 0
                do j = 1,iunit
                   rxnew(j) = rxnew(ifirst) - (rxu(self,ifirst) - rxu(self,j))
                   rynew(j) = rynew(ifirst) - (ryu(self,ifirst) - ryu(self,j))
                   rznew(j) = rznew(ifirst) - (rzu(self,ifirst) - rzu(self,j))
                end do
                if(nsampos(iparty).eq.2) then
                   second_bead = iprev
                   old = from(o_type)
                   xtarget(1) = rxu(other,old)
                   ytarget(1) = ryu(other,old)
                   ztarget(1) = rzu(other,old)
                   old = prev(o_type)
                   xtarget(2) = rxu(other,old)
                   ytarget(2) = ryu(other,old)
                   ztarget(2) = rzu(other,old)
                   call align_lines(first_bead,second_bead,iunit,rxnew,rynew,rznew,xtarget,ytarget,ztarget)
                end if
             end if
             ! end if
          !end if
          !> \bug problem if nsampos.eq.iunit, but rindex.gt.0
             call schedule(igrow,imolty,islen,ifirst,iprev,4,0)
             !cc---END JLR 11-24-09
          else
             call schedule(igrow,imolty,islen,ifirst,iprev,3,0)
          end if

          if (lgrand.and.boxa.ne.boxb) then
             Rosenbluth_normalization(imolty)=(real(nchoi(imolty),dp)**islen)*real(nchoih(imoltb),dp)
             if (lrigid(imolty).and.nsampos(iparty).lt.3) Rosenbluth_normalization(imolty)=Rosenbluth_normalization(imolty)*real(nchoir(imolta),dp)
          end if

          !> \bug Need to check: I wonder if this works with lrigid?
          ! Adding in multiple end regrowths
          if (ncut(iparty,s_type).gt.1) then
             do izz = 2,ncut(iparty,s_type)
                ifirst = from(s_type+2*(izz-1))
                iprev = prev(s_type+2*(izz-1))
                call schedule(igrow,imolty,islen,ifirst,iprev,5,0)

                if (lgrand.and.boxa.ne.boxb) Rosenbluth_normalization(imolty)=Rosenbluth_normalization(imolty)*(real(nchoi(imolty),dp)**islen)
             end do
          end if

          ! Paul -- SAFE-SWATCH: if there are beads needed to be regrown using SAFE-CBMC
          ! set lexshed to be false to avoid the calculation of intramolecular LJ interactions
          if (ncutsafe(iparty,s_type) .gt. 0) then
             lexshed = .false.
             do icbu = 1, nsampos(iparty)
                jj = splist(iparty,icbu,s_type)
                lexshed(jj) = .true.
             end do
          end if

!******************
! CBMC - new growth
!******************
          if (boxa.eq.boxb) then
             ! moving molecules for rosenbluth
             if (ic.eq.1) then
                nboxi(other) = 0
                ! putting molecule b in position 1
                ! do izz = 1, nsampos(iparty)
                !    bdmol_b = splist(iparty,izz,type_b)
                !    rxu(other,bdmol_b) = rxut(3,izz)
                !    ryu(other,bdmol_b) = ryut(3,izz)
                !    rzu(other,bdmol_b) = rzut(3,izz)
                ! end do
             else
                ! putting molecule a into its (fully grown) trial position 2
                do izz = 1, nunit(moltyp(other))
                   rxu(other,izz) = rxut(1,izz)
                   ryu(other,izz) = ryut(1,izz)
                   rzu(other,izz) = rzut(1,izz)
                end do
             end if
          end if

          if (lgrand.and.iboxnew.ne.1) goto 1000

          ! grow molecules
          ! changing for lrigid to include waddnew
          waddnew = 1.0E0_dp

          ! JLR 11-24-09 New stuff for rigid swatch
          ! Different logic/calls to rosenbluth for rigid molecules
          if (lrigid(imolty)) then
             if (nsampos(iparty).ge.iunit) then
                !molecule is all there
                !don't regrow anything in rosenbluth
                icallrose = 4
                !> \bug why don't do rigrot if ifirst.ge.riutry (ifirst is part of rigid beads)?
                ! else if ((nsampos(iparty).ge.3).or.(rindex(imolty).gt.0.and.ifirst.ge.riutry(imolty,1))) then
             else if (nsampos(iparty).ge.3) then
                ! rigid part is grown, don't do rigrot in rosebluth
                icallrose = 3
             else if (nsampos(iparty).eq.2) then !BX
                ! rigid part is not grown, do rigrot
                icallrose = 5
             else
                icallrose = 0
             end if
          else
             ! flexible molecule call rosenbluth in normal fashion
             icallrose = 2
          end if
          ! END JLR 11-24-09

          ! grow new chain conformation
          !> \bug why call with other and self instead of self and self as for interbox swatch?
          if (boxa.eq.boxb) then
             call rosenbluth(.true.,lterm,self,self,imolty,islen,boxa,igrow,waddnew,.false.,vdum2,icallrose,first_bead,second_bead)
          else
             call rosenbluth(.true.,lterm,other,self,imolty,islen,iboxnew,igrow,waddnew,.false.,vdum2,icallrose,first_bead,second_bead)
          end if

          if (boxa.eq.boxb) then
             ! moving molecules back
             if (ic.eq.1) then
                nboxi(other) = iboxnew
                ! do izz = 1, nsampos(iparty)
                !    bdmol_b = splist(iparty,izz,type_b)
                !    rxu(other,bdmol_b) = rxut(4,izz)
                !    ryu(other,bdmol_b) = ryut(4,izz)
                !    rzu(other,bdmol_b) = rzut(4,izz)
                ! end do
             else
                do izz = 1, nunit(moltyp(other))
                   rxu(other,izz) = rx_1(izz)
                   ryu(other,izz) = ry_1(izz)
                   rzu(other,izz) = rz_1(izz)
                end do
             end if
          end if

          ! termination of cbmc attempt due to walk termination
          if (lterm) then
             ! if (boxa.eq.boxb) liswatch = .false.
             if (lgrand) then
                if (boxa.ne.1) then
                   call gcmc_cleanup(imolta,boxa)
                else if (boxb.ne.1) then
                   call gcmc_cleanup(imoltb,boxb)
                end if
             end if
             return
          end if

          ! propagate new rosenbluth weight
          tweight = tweight*weight*waddnew

      !*** end rigid add on ***

      ! make a copy of original coords

      do izz = 1, nunit(moltyp(self))
         rxut(ic,izz) = rxu(self,izz)
         ryut(ic,izz) = ryu(self,izz)
         rzut(ic,izz) = rzu(self,izz)
      end do

      ! Copy the rosenbluth-generated coords to "real (r*u)" coords

      do izz = 1, nunit(moltyp(self))
         rxu(self,izz) = rxnew(izz)
         ryu(self,izz) = rynew(izz)
         rzu(self,izz) = rznew(izz)
      end do

      ! Calculated U_bonded

      call U_bonded(self,imolty,tmpvvib,tmpvbend,tmpvtg)

      ! Copy back the "real (r*u)" coords

      do izz = 1, nunit(moltyp(self))
         rxu(self,izz) = rxut(ic,izz)
         ryu(self,izz) = ryut(ic,izz)
         rzu(self,izz) = rzut(ic,izz)
      end do

          ! save the new coordinates
          do jj = 1,igrow
             rxut(ic,jj) = rxnew(jj)
             ryut(ic,jj) = rynew(jj)
             rzut(ic,jj) = rznew(jj)
          end do

          ! Corrections for switched beads, and DC-CBMC
          ! Assign all of the grown new and old beads to rxuion
          ! with rxuion: new = 2
          iii = 2
          do j=1,igrow
             rxuion(j,iii) = rxnew(j)
             ryuion(j,iii) = rynew(j)
             rzuion(j,iii) = rznew(j)
             qquion(j,iii) = qqu(self,j)
          end do

          ! added from new-combined code
          if (iunit.ne.igrow) then
             ! for explicit-hydrogen model, put on the hydrogens
             ! use phony number iins and call explct to add constrained hydrogens
             iins = nchain + 1
             moltyp(iins) = imolty
             do j=1,igrow
                rxu(iins,j) = rxnew(j)
                ryu(iins,j) = rynew(j)
                rzu(iins,j) = rznew(j)
             end do
             call explct(iins,vdum,.false.,.false.)
             do j = igrow + 1, iunit
                rxuion(j,iii) = rxu(iins,j)
                ryuion(j,iii) = ryu(iins,j)
                rzuion(j,iii) = rzu(iins,j)
                qquion(j,iii) = qqu(self,j)
                rxut(ic,j) = rxu(iins,j)
                ryut(ic,j) = ryu(iins,j)
                rzut(ic,j) = rzu(iins,j)
             end do
          end if

          ! DC-CBMC correction only if there is no SAFE-SWATCH!
          ! Begin DC-CBMC, explicit-hydrogen and
          ! switched bead corrections for NEW configuration
          ! Calculate the true site-site energy

          if (ncutsafe(iparty,s_type) .eq. 0) then
             if (boxa.eq.boxb) then
                if (ic.eq.1) then
                   ! exclude molecule b from energy calculation, put in other box *
                   nboxi(imolb) = 0
                else
                   ! put molecule a into position 2 (fully grown trial position)
                   ! for the energy of b's new position (second time around)
                   do jj = 1,nunit(imolta)
                      rxu(imola,jj) = rxut(1,jj)
                      ryu(imola,jj) = ryut(1,jj)
                      rzu(imola,jj) = rzut(1,jj)
                   end do
                end if
             end if

             !> \bug energy problem for cases which involve the change of the bending type
             !> and torsional type for those units swatched!!!!!!
             if (boxa.eq.boxb) then
                call energy(self,imolty,v,iii,boxa,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)
             else
                call energy(other,imolty,v,iii,iboxnew,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)
             end if

             if (boxa.eq.boxb) then
                ! return to normal
                if (ic.eq.1) then
                   ! return b to original box
                   nboxi(imolb) = boxb
                else
                   ! return a to position 1
                   do jj = 1, nunit(imolta)
                      rxu(imola,jj) = rx_1(jj)
                      ryu(imola,jj) = ry_1(jj)
                      rzu(imola,jj) = rz_1(jj)
                   end do
                end if
             end if

             if (lterm) then
             ! call err_exit(__FILE__,__LINE__,'interesting screwup in CBMC swatch',myid+1)
                if (lgrand) then
                   if (boxa.ne.1) then
                      call gcmc_cleanup(imolta,boxa)
                   else if (boxb.ne.1) then
                      call gcmc_cleanup(imoltb,boxb)
                   end if
                end if
                return
             end if

             ! add on the changes in energy
             delen = v(ivTot) - ( vnew(ivTot) - (tmpvvib+tmpvbend+tmpvtg) )
             tweight = tweight*exp(-beta*delen)

             vnew(ivTot) = vnew(ivTot) + delen
             vnew(ivInterLJ) = v(ivInterLJ)
             vnew(ivIntraLJ) = v(ivIntraLJ)
             vnew(ivExt) = v(ivExt)
             vnew(ivElect) = v(ivElect)
             vnew(ivEwald)= v(ivEwald)
             ! End DC-CBMC and switched bead Corrections for NEW configuration

             ! save the trial energies
             vnbox(ivTot,iboxnew)   = vnbox(ivTot,iboxnew)  + vnew(ivTot)
             vnbox(ivInterLJ,iboxnew)  = vnbox(ivInterLJ,iboxnew) + vnew(ivInterLJ)
             vnbox(ivIntraLJ,iboxnew)  = vnbox(ivIntraLJ,iboxnew) + vnew(ivIntraLJ)
             vnbox(ivStretching,iboxnew)  = vnbox(ivStretching,iboxnew) + tmpvvib
             vnbox(ivBending,iboxnew) = vnbox(ivBending,iboxnew) + tmpvbend
             vnbox(ivTorsion,iboxnew) = vnbox(ivTorsion,iboxnew) + tmpvtg
             vnbox(ivExt,iboxnew)  = vnbox(ivExt,iboxnew) + vnew(ivExt)
             vnbox(ivElect,iboxnew) = vnbox(ivElect,iboxnew)+ vnew(ivElect)
             vnbox(ivEwald,iboxnew) = vnbox(ivEwald,iboxnew)+ vnew(ivEwald)
          else
             ! if there is still SAFE-CBMC move, temporarily store energies in vnewtemp
             vnewtemp = vnew
          end if

          if (lgrand.and.iboxold.ne.1) cycle oldnew

!******************
! CBMC - old growth
!******************
          ! rigid add on
1000      waddold = 1.0E0_dp

          ! grow old chain conformation
          ! JLR 11-24-09 New stuff for rigid swatch
          if (lrigid(imolty)) then
             if (nsampos(iparty).ge.iunit) then
                !molecule is all there
                !don't regrow anything in rosenbluth
                icallrose = 4
                !> \bug why don't do rigrot if ifirst.ge.riutry (ifirst is part of rigid beads)?
                ! else if ((nsampos(iparty).ge.3).or.(rindex(imolty).gt.0.and.ifirst.ge.riutry(imolty,1))) then
             else if (nsampos(iparty).ge.3) then
                ! rigid part is grown, don't do rigrot in rosebluth
                icallrose = 3
             else if (nsampos(iparty).eq.2) then !BX
                ! rigid part is not grown, do rigrot
                icallrose = 5
             else
                icallrose = 0
             end if
          else
             ! flexible molecule call rosenbluth in normal fashion
             icallrose = 2
          end if

          call rosenbluth(.false.,lterm,self,self,imolty,islen,iboxold,igrow,waddold,.false.,vdum2,icallrose,first_bead,second_bead)
          ! END JLR 11-24-09

          ! termination of old walk due to problems generating orientations
          if (lterm) then
             write(io_output,*) 'SWATCH: old growth rejected'
             ! if (boxa.eq.boxb) liswatch = .false.
             if (lgrand) then
                if (boxa.ne.1) then
                   call gcmc_cleanup(imolta,boxa)
                else if (boxb.ne.1) then
                   call gcmc_cleanup(imoltb,boxb)
                end if
             end if
             return
          end if

          ! propagate old rosenbluth weight
          tweiold = tweiold*weiold*waddold

          ! end rigid add on

      ! Calculated U_bonded

      call U_bonded(self,imolty,tmpvvib,tmpvbend,tmpvtg)

          ! store the old grown beads and explict placed beads positions
          ! 1 = old conformation
          iii = 1
          do j = 1,iunit
             rxuion(j,1) = rxu(self,j)
             ryuion(j,1) = ryu(self,j)
             rzuion(j,1) = rzu(self,j)
             qquion(j,1) = qqu(self,j)
          end do

          ! Begin corrections for DC-CBMC and switched beads for OLD configuration
          ! correct the acceptance rules
          ! calculate the full rcut site-site energy
          ! only when there is no SAFE-SWATCH

          if (ncutsafe(iparty,s_type) .eq. 0) then
             ! excluding molecule b for first loop
             if (ic.eq.1.and.boxa.eq.boxb) then
                nboxi(imolb) = 0
             end if

             ! get total energy
             call energy(self,imolty,v,iii,iboxold,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)

             if (ic.eq.1.and.boxa.eq.boxb) then
                ! return b to current box
                nboxi(imolb) = boxb
             end if

             if (lterm) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf SWATCH',myid+1)

              deleo = v(ivTot) - ( vold(ivTot) - (tmpvvib+tmpvbend+tmpvtg) )
             tweiold = tweiold*exp(-beta*deleo)

             vold(ivTot) = vold(ivTot) + deleo
             vold(ivIntraLJ) = v(ivIntraLJ)
             vold(ivInterLJ) = v(ivInterLJ)
             vold(ivExt) = v(ivExt)
             vold(ivElect) = v(ivElect)
             vold(ivEwald)= v(ivEwald)
             ! End Correction for DC-CBMC and switched beads for OLD configuration

             ! save the trial energies
             vnbox(ivTot,iboxold)   = vnbox(ivTot,iboxold)  - vold(ivTot)
             vnbox(ivInterLJ,iboxold)  = vnbox(ivInterLJ,iboxold) - vold(ivInterLJ)
             vnbox(ivIntraLJ,iboxold)  = vnbox(ivIntraLJ,iboxold) - vold(ivIntraLJ)
             vnbox(ivStretching,iboxold) = vnbox(ivStretching,iboxold) - tmpvvib
             vnbox(ivBending,iboxold) = vnbox(ivBending,iboxold)    - tmpvbend
             vnbox(ivTorsion,iboxold) = vnbox(ivTorsion,iboxold)    - tmpvtg
             vnbox(ivExt,iboxold)  = vnbox(ivExt,iboxold) - vold(ivExt)
             vnbox(ivElect,iboxold) = vnbox(ivElect,iboxold)- vold(ivElect)
             vnbox(ivEwald,iboxold) = vnbox(ivEwald,iboxold)- vold(ivEwald)
          else
             ! if there is still SAFE-CBMC move, temporarily store energies in voldtemp
             voldtemp = vold
          end if

      end if ! if there is regular SAFE-CBMC part


!*******************************
! SAFE-CBMC - new and old growth
!*******************************
       ! Paul -- now grow SAFE-CBMC part
       waddnew = 1.0E0_dp

       do izz = 1, ncutsafe(iparty,s_type)

          if (boxa.eq.boxb) then
             ! moving molecules for rosenbluth
             if (ic.eq.1) then
                nboxi(other) = 0
             else
                ! putting molecule a into its (fully grown) trial position 2
                do jj = 1, nunit(moltyp(other))
                   rxu(other,jj) = rxut(1,jj)
                   ryu(other,jj) = ryut(1,jj)
                   rzu(other,jj) = rzut(1,jj)
                end do
             end if
          end if

          ifirst = fromsafe(s_type+2*(izz-1))
          iprev = prevsafe(s_type+2*(izz-1))
          iindex = indexsafe(s_type+2*(izz-1))
          call safeschedule(igrow,imolty,islen,ifirst,iindex+1,5,iprev)

          if (boxa.eq.boxb) then
              call rosenbluth(.true.,lterm,self,self,imolty,islen,boxa,igrow,waddnew,.true.,cwtorfn,2)
          else
              call rosenbluth(.true.,lterm,other,self,imolty,islen,iboxnew,igrow,waddnew,.true.,cwtorfn,2)
          end if

          if (boxa.eq.boxb) then
              ! moving molecules back
              if (ic.eq.1) then
                 nboxi(other) = iboxnew
              else
                 do jj = 1, nunit(moltyp(other))
                    rxu(other,jj) = rx_1(jj)
                    ryu(other,jj) = ry_1(jj)
                    rzu(other,jj) = rz_1(jj)
                 end do
              end if
          end if

          ! termination of cbmc attempt due to walk termination
          if (lterm) then
              ! if (boxa.eq.boxb) liswatch = .false.
              if (lgrand) then
                 if (boxa.ne.1) then
                    call gcmc_cleanup(imolta,boxa)
                 else if (boxb.ne.1) then
                    call gcmc_cleanup(imoltb,boxb)
                 end if
              end if
              return
           end if

           ! propagate new rosenbluth weight
           tweight = tweight*weight*waddnew/cwtorfn

           ! save the new coordinates
           do jj = 1,igrow
              rxut(ic,jj) = rxnew(jj)
              ryut(ic,jj) = rynew(jj)
              rzut(ic,jj) = rznew(jj)
           end do

           ! Corrections for switched beads, and DC-CBMC
           ! Assign all of the grown new and old beads to rxuion
           ! with rxuion: new = 2
           iii = 2
           do j=1,igrow
              rxuion(j,iii) = rxnew(j)
              ryuion(j,iii) = rynew(j)
              rzuion(j,iii) = rznew(j)
              qquion(j,iii) = qqu(self,j)
           end do

           ! added from new-combined code
           if (iunit.ne.igrow) then
              ! for explicit-hydrogen model, put on the hydrogens
              ! use phony number iins and call explct to add constrained hydrogens
              iins = nchain + 1
              moltyp(iins) = imolty
              do j=1,igrow
                 rxu(iins,j) = rxnew(j)
                 ryu(iins,j) = rynew(j)
                 rzu(iins,j) = rznew(j)
              end do
              call explct(iins,vdum,.false.,.false.)
              do j = igrow + 1, iunit
                 rxuion(j,iii) = rxu(iins,j)
                 ryuion(j,iii) = ryu(iins,j)
                 rzuion(j,iii) = rzu(iins,j)
                 qquion(j,iii) = qqu(self,j)
                 rxut(ic,j) = rxu(iins,j)
                 ryut(ic,j) = ryu(iins,j)
                 rzut(ic,j) = rzu(iins,j)
              end do
           end if

           ! accumulates energies in vnewtemp
           vnewtemp = vnewtemp + vnew

           ! Begin DC-CBMC, explicit-hydrogen and
           ! switched bead corrections for NEW configuration
           ! Calculate the true site-site energy
           ! Only when it is the last step of SAFE-SWATCH!

           if (izz .eq. ncutsafe(iparty,s_type)) then
              if (boxa.eq.boxb) then
                 if (ic.eq.1) then
                    ! exclude molecule b from energy calculation, put in other box *
                    nboxi(imolb) = 0
                 else
                    ! put molecule a into position 2 (fully grown trial position)
                    ! for the energy of b's new position (second time around)
                    do jj = 1,nunit(imolta)
                       rxu(imola,jj) = rxut(1,jj)
                       ryu(imola,jj) = ryut(1,jj)
                       rzu(imola,jj) = rzut(1,jj)
                    end do
                 end if
              end if

              !> \bug energy problem for cases which involve the change of the bending type
              !> and torsional type for those units swatched!!!!!!
              if (boxa.eq.boxb) then
                 call energy(self,imolty,v,iii,boxa,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)
              else
                 call energy(other,imolty,v,iii,iboxnew,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)
              end if

              if (boxa.eq.boxb) then
                 ! return to normal
                 if (ic.eq.1) then
                    ! return b to original box
                    nboxi(imolb) = boxb
                 else
                    ! return a to position 1
                    do jj = 1, nunit(imolta)
                       rxu(imola,jj) = rx_1(jj)
                       ryu(imola,jj) = ry_1(jj)
                       rzu(imola,jj) = rz_1(jj)
                    end do
                 end if
              end if

              if (lterm) then
                 ! call err_exit(__FILE__,__LINE__,'interesting screwup in CBMC swatch',myid+1)
                 if (lgrand) then
                    if (boxa.ne.1) then
                       call gcmc_cleanup(imolta,boxa)
                    else if (boxb.ne.1) then
                       call gcmc_cleanup(imoltb,boxb)
                    end if
                 end if
                 return
              end if

              ! use vnewtemp which is accumulated from CBMC and SAFE-CBMC
              vnew = vnewtemp

              ! add on the changes in energy
              delen = v(ivTot) - ( vnew(ivTot) - (vnew(ivStretching) + vnew(ivBending) + vnew(ivTorsion)) )
              tweight = tweight*exp(-beta*delen)

              vnew(ivTot) = vnew(ivTot) + delen
              vnew(ivInterLJ) = v(ivInterLJ)
              vnew(ivIntraLJ) = v(ivIntraLJ)
              vnew(ivExt) = v(ivExt)
              vnew(ivElect) = v(ivElect)
              vnew(ivEwald)= v(ivEwald)
              ! End DC-CBMC and switched bead Corrections for NEW configuration

              ! save the trial energies
              vnbox(ivTot,iboxnew)   = vnbox(ivTot,iboxnew)  + vnew(ivTot)
              vnbox(ivInterLJ,iboxnew)  = vnbox(ivInterLJ,iboxnew) + vnew(ivInterLJ)
              vnbox(ivIntraLJ,iboxnew)  = vnbox(ivIntraLJ,iboxnew) + vnew(ivIntraLJ)
              vnbox(ivStretching,iboxnew)  = vnbox(ivStretching,iboxnew) + vnew(ivStretching)
              vnbox(ivTorsion,iboxnew)   = vnbox(ivTorsion,iboxnew)  + vnew(ivTorsion)
              vnbox(ivExt,iboxnew)  = vnbox(ivExt,iboxnew) + vnew(ivExt)
              vnbox(ivBending,iboxnew)  = vnbox(ivBending,iboxnew) + vnew(ivBending)
              vnbox(ivElect,iboxnew) = vnbox(ivElect,iboxnew)+ vnew(ivElect)
              vnbox(ivEwald,iboxnew) = vnbox(ivEwald,iboxnew)+ vnew(ivEwald)
           end if

           if (lgrand.and.iboxold.ne.1) cycle oldnew

           ! SAFE-CBMC old growth
           ! rigid add on
           waddold = 1.0E0_dp

           call rosenbluth(.false.,lterm,self,self,imolty,islen,iboxold,igrow,waddold,.true.,cwtorfo,2)

           ! termination of old walk due to problems generating orientations
           if (lterm) then
               write(io_output,*) 'SAFE-SWATCH: old growth rejected'
               ! if (boxa.eq.boxb) liswatch = .false.
               if (lgrand) then
                   if (boxa.ne.1) then
                       call gcmc_cleanup(imolta,boxa)
                   else if (boxb.ne.1) then
                       call gcmc_cleanup(imoltb,boxb)
                   end if
               end if
               return
           end if

           ! propagate old rosenbluth weight
           tweiold = tweiold*weiold*waddold/cwtorfo
           ! end rigid add on

           ! store the old grown beads and explict placed beads positions
           ! 1 = old conformation
           iii = 1
           do j = 1,iunit
              rxuion(j,1) = rxu(self,j)
              ryuion(j,1) = ryu(self,j)
              rzuion(j,1) = rzu(self,j)
              qquion(j,1) = qqu(self,j)
           end do

           ! accumulates energies in voldtemp
           voldtemp = voldtemp + vold

           ! Begin corrections for DC-CBMC and switched beads for OLD configuration
           ! correct the acceptance rules
           ! calculate the full rcut site-site energy
           ! only when this is the last SAFE-SWATCH

           if (izz .eq. ncutsafe(iparty,s_type)) then
              ! excluding molecule b for first loop
              if (ic.eq.1.and.boxa.eq.boxb) then
                 nboxi(imolb) = 0
              end if

              ! get total energy
              call energy(self,imolty,v,iii,iboxold,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)

              if (ic.eq.1.and.boxa.eq.boxb) then
                 ! return b to current box
                 nboxi(imolb) = boxb
              end if

              if (lterm) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf SWATCH',myid+1)

              ! use voldtemp which is accumulated from CBMC and SAFE-CBMC
              vold = voldtemp

              deleo = v(ivTot) - ( vold(ivTot) - (vold(ivStretching) + vold(ivBending) + vold(ivTorsion)) )
              tweiold = tweiold*exp(-beta*deleo)

              vold(ivTot) = vold(ivTot) + deleo
              vold(ivIntraLJ) = v(ivIntraLJ)
              vold(ivInterLJ) = v(ivInterLJ)
              vold(ivExt) = v(ivExt)
              vold(ivElect) = v(ivElect)
              vold(ivEwald)= v(ivEwald)
              ! End Correction for DC-CBMC and switched beads for OLD configuration

              ! save the trial energies
              vnbox(ivTot,iboxold)   = vnbox(ivTot,iboxold)  - vold(ivTot)
              vnbox(ivInterLJ,iboxold)  = vnbox(ivInterLJ,iboxold) - vold(ivInterLJ)
              vnbox(ivIntraLJ,iboxold)  = vnbox(ivIntraLJ,iboxold) - vold(ivIntraLJ)
              vnbox(ivStretching,iboxold)  = vnbox(ivStretching,iboxold) - vold(ivStretching)
              vnbox(ivTorsion,iboxold)   = vnbox(ivTorsion,iboxold)  - vold(ivTorsion)
              vnbox(ivExt,iboxold)  = vnbox(ivExt,iboxold) - vold(ivExt)
              vnbox(ivBending,iboxold)  = vnbox(ivBending,iboxold) - vold(ivBending)
              vnbox(ivElect,iboxold) = vnbox(ivElect,iboxold)- vold(ivElect)
              vnbox(ivEwald,iboxold) = vnbox(ivEwald,iboxold)- vold(ivEwald)
           end if
       end do
    end do oldnew

    ! Perform the Ewald sum reciprical space corrections
    if (lewald) then
       ! added into tweight even though it really contains new-old
       ! Box A
       if (.not.lideal(boxa)) then
          ! store the reciprocal space vector
          if (boxa.eq.boxb) call recip(boxa,vdum,vdum,3)

          ! Position 1 ***
          oldchain = imola
          newchain = imolb
          oldunit = nunit(imolta)
          newunit = nunit(imoltb)

          do j = 1,oldunit
             rxuion(j,1) = rxu(oldchain,j)
             ryuion(j,1) = ryu(oldchain,j)
             rzuion(j,1) = rzu(oldchain,j)
             qquion(j,1) = qqu(oldchain,j)
          end do
          moltion(1) = imolta
          do j = 1,newunit
             rxuion(j,2) = rxut(2,j)
             ryuion(j,2) = ryut(2,j)
             rzuion(j,2) = rzut(2,j)
             qquion(j,2) = qqu(newchain,j)
          end do
          moltion(2) = imoltb

          call recip(boxa,vrecipn,vrecipo,1)

          delen = vrecipn - vrecipo
          tweight = tweight * exp(-beta*delen)

          vnbox(ivEwald,boxa) = vnbox(ivEwald,boxa) + delen
          vnbox(ivTot,boxa) = vnbox(ivTot,boxa) + delen

          ! update the reciprocal space terms *
          if (boxa.eq.boxb) call recip(boxa,vdum,vdum,2)
       end if
       if (.not.lideal(boxb)) then
          ! Box B
          oldchain = imolb
          newchain = imola
          oldunit = nunit(imoltb)
          newunit = nunit(imolta)

          do j = 1,oldunit
             rxuion(j,1) = rxu(oldchain,j)
             ryuion(j,1) = ryu(oldchain,j)
             rzuion(j,1) = rzu(oldchain,j)
             qquion(j,1) = qqu(oldchain,j)
          end do
          moltion(1) = imoltb
          do j = 1,newunit
             rxuion(j,2) = rxut(1,j)
             ryuion(j,2) = ryut(1,j)
             rzuion(j,2) = rzut(1,j)
             qquion(j,2) = qqu(newchain,j)
          end do
          moltion(2) = imolta

          call recip(boxb,vrecipn,vrecipo,1)

          delen = vrecipn - vrecipo
          tweight = tweight * exp(-beta*delen)

          vnbox(ivEwald,boxb) = vnbox(ivEwald,boxb) + delen
          vnbox(ivTot,boxb) = vnbox(ivTot,boxb) + delen
       end if
    end if
    ! End Ewald-sum Corrections

    if (ltailc.and.boxa.ne.boxb) then
       ! add tail corrections
       if (lpbcz) then
          if (lsolid(boxa) .and. .not. lrect(boxa)) then
             vola = (hmat(boxa,1) * (hmat(boxa,5) * hmat(boxa,9) - hmat(boxa,8)*hmat(boxa,6))+ hmat(boxa,4)*(hmat(boxa,8) * hmat(boxa,3)-hmat(boxa,2)* hmat(boxa,9))+hmat(boxa,7) * (hmat(boxa,2)*hmat(boxa,6)- hmat(boxa,5)*hmat(boxa,3)))
          else
             vola=boxlx(boxa)*boxly(boxa)*boxlz(boxa)
          end if

          if (lsolid(boxb) .and. .not. lrect(boxb)) then
             volb = (hmat(boxb,1) * (hmat(boxb,5) * hmat(boxb,9) - hmat(boxb,8)*hmat(boxb,6))+ hmat(boxb,4)*(hmat(boxb,8) * hmat(boxb,3)-hmat(boxb,2)* hmat(boxb,9))+hmat(boxb,7) * (hmat(boxb,2)*hmat(boxb,6)- hmat(boxb,5)*hmat(boxb,3)))
          else
             volb=boxlx(boxb)*boxly(boxb)*boxlz(boxb)
          end if
       else
          vola=boxlx(boxa)*boxly(boxa)
          volb=boxlx(boxb)*boxly(boxb)
       end if

       ! for new BOXINS with inserted particle
       do mm = 1,2
          dinsta = 0.0E0_dp
          if ( mm .eq. 1 ) then
             ibox   = boxa
             imolin = imoltb
             dvol   = vola
             imolrm = imolta
          else
             ibox   = boxb
             imolin = imolta
             dvol   = volb
             imolrm = imoltb
          end if
          ! JLR 11-24-09 don't do tail corrections for ideal box
          if (.not.lideal(ibox)) then
             ! new logic for tail correction (same answer) MGM 3-25-98
             do jmt = 1, nmolty
                rho = real( ncmt(ibox,jmt) ,dp)
                if ( jmt .eq. imolin ) rho = rho + 1.0E0_dp
                if ( jmt .eq. imolrm ) rho = rho - 1.0E0_dp
                rho = rho / dvol
                do imt = 1, nmolty
                   dicount = ncmt(ibox,imt)
                   if ( imt .eq. imolin ) dicount = dicount + 1
                   if ( imt .eq. imolrm ) dicount = dicount - 1
                   dinsta = dinsta +  dicount * coru(imt,jmt,rho,ibox)
                end do
             end do
             dinsta = dinsta - vbox(ivTail,ibox)

             tweight=tweight*exp(-beta*dinsta)

             vnbox(ivTail,ibox) = dinsta
          end if
          ! END JLR 11-24-09
       end do
    else
       vnbox(ivTail,boxa) = 0.0E0_dp
       vnbox(ivTail,boxb) = 0.0E0_dp
    end if

    wnlog = log10( tweight )
    wolog = log10( tweiold )
    wdlog = wnlog - wolog
    if ( wdlog .lt. -softcut ) then
       ! write(io_output,*) '### underflow in wratio calculation ###'
       if (.not.lideal(boxa).and.boxa.eq.boxb) call recip(boxa,vdum,vdum,4)
       if (lgrand) then
          if (boxa.ne.1) then
             call gcmc_cleanup(imolta,boxa)
          else if (boxb.ne.1) then
             call gcmc_cleanup(imoltb,boxb)
          end if
       end if
       return
    end if

    wswat = tweight / tweiold

    if (boxa.ne.boxb) then
       if (lgrand) then
          if (boxa.eq.1) then
             wswat = wswat*Rosenbluth_normalization(imolta)/Rosenbluth_normalization(imoltb)*B(imoltb)/B(imolta)*orgaia/real(orgbia+1,dp)
          else
             wswat = wswat*Rosenbluth_normalization(imoltb)/Rosenbluth_normalization(imolta)*B(imolta)/B(imoltb)*orgbib/real(orgaib+1,dp)
          end if
       else
          wswat = wswat*real(orgaia*orgbib,dp)/real((orgbia+1)*(orgaib+1),dp)*exp(beta*(eta2(boxa,imolta)+eta2(boxb,imoltb)-eta2(boxa,imoltb)-eta2(boxb,imolta)))
       end if

       if (lopt_bias(imolta)) call update_bias(log(wswat*2.0)/beta/2.0,boxa,boxb,imolta)
       if (lopt_bias(imoltb)) call update_bias(log(wswat*2.0)/beta/2.0,boxb,boxa,imoltb)
    end if

    if (random(-1).le.wswat) then
       ! we can now accept !!!!!
       bsswat(iparty,ipairb,boxa) = bsswat(iparty,ipairb,boxa) + 1

       total_tor =0.0_dp
       total_bend=0.0_dp
       total_vib =0.0_dp
!       write(*,*) rindex(imolta), rindex(imoltb), nsampos(iparty),nunit(imolta), nunit(imoltb)
       if (lrigid(imolta).and.((rindex(imolta).gt.0.and.nsampos(iparty).eq.nunit(imolta)).or.pm_atom_tra.gt.0)) then
          call U_bonded(imola,imolta,vvibn,vbendn,vtgn)
          total_tor =vtgn
          total_bend=vbendn
          total_vib =vvibn
       end if
       if (lrigid(imoltb).and.((rindex(imoltb).gt.0.and.nsampos(iparty).eq.nunit(imoltb)).or.pm_atom_tra.gt.0)) then
          call U_bonded(imolb,imoltb,vvibn,vbendn,vtgn)
          total_tor  = total_tor - vtgn
          total_bend = total_bend - vbendn
          total_vib  = total_vib - vvibn
       end if
       total_NBE = total_tor + total_bend + total_vib

       do jj=1,2
          if (jj.eq.1) then
             if (lgrand.and.boxa.ne.1) cycle
             ic = boxa
          else if (jj.eq.2) then
             if ((boxa.eq.boxb).or.(lgrand.and.boxb.ne.1)) exit
             ic = boxb
             total_tor =-total_tor
             total_bend=-total_bend
             total_vib =-total_vib
             total_NBE =-total_NBE
          end if
!          write(*,*) 'total tor diff', total_tor
          vbox(ivTot,ic) = vbox(ivTot,ic) + vnbox(ivTot,ic)   + vnbox(ivTail,ic) - total_NBE
          vbox(ivInterLJ,ic) = vbox(ivInterLJ,ic) + vnbox(ivInterLJ,ic)  + vnbox(ivTail,ic)
          vbox(ivIntraLJ,ic) = vbox(ivIntraLJ,ic) + vnbox(ivIntraLJ,ic)
          vbox(ivStretching,ic) = vbox(ivStretching,ic) + vnbox(ivStretching,ic)  - total_vib
          vbox(ivTorsion,ic) = vbox(ivTorsion,ic) + vnbox(ivTorsion,ic)   - total_tor
          vbox(ivExt,ic) = vbox(ivExt,ic) + vnbox(ivExt,ic)
          vbox(ivBending,ic) = vbox(ivBending,ic) + vnbox(ivBending,ic)  - total_bend
          vbox(ivTail,ic) = vbox(ivTail,ic) + vnbox(ivTail,ic)
          vbox(ivElect,ic) = vbox(ivElect,ic) + vnbox(ivElect,ic) + vnbox(ivEwald,ic)
       end do

       ! Update coordinates in kdtree for A
       if ((.not. lcutcm) .and. lkdtree .and. (lkdtree_box(boxa) .or. lkdtree_box(boxb))) then
           do ic = 1, iunita
               rxu_update(ic) = rxut(1, ic)
               ryu_update(ic) = ryut(1, ic)
               rzu_update(ic) = rzut(1, ic)
           end do

           call update_coord_in_tree(imola, iunita, boxa, boxb, .true., .false.)

           do ic = 1, iunitb
               rxu_update(ic) = rxut(2, ic)
               ryu_update(ic) = ryut(2, ic)
               rzu_update(ic) = rzut(2, ic)
           end do

           call update_coord_in_tree(imolb, iunitb, boxb, boxa, .true., .false.)
       end if

       ! Update coordinates in r*u arrays
       do ic = 1,iunita
          rxu(imola,ic) = rxut(1,ic)
          ryu(imola,ic) = ryut(1,ic)
          rzu(imola,ic) = rzut(1,ic)
       end do
       do ic = 1,iunitb
          rxu(imolb,ic) = rxut(2,ic)
          ryu(imolb,ic) = ryut(2,ic)
          rzu(imolb,ic) = rzut(2,ic)
       end do

       ! update book keeping
       nboxi(imola) = boxb
       nboxi(imolb) = boxa

       if (boxa.ne.boxb) then
          parbox(orgbia+1,boxa,imoltb)= imolb
          parbox(orgaib+1,boxb,imolta)= imola
          parbox(iboxal,boxa,imolta) = parbox(orgaia,boxa,imolta)
          parbox(iboxbl,boxb,imoltb) = parbox(orgbib,boxb,imoltb)
          parbox(orgaia,boxa,imolta) = 0
          parbox(orgbib,boxb,imoltb) = 0

          ncmt(boxa,imolta) = orgaia - 1
          ncmt(boxa,imoltb) = orgbia + 1
          ncmt(boxb,imolta) = orgaib + 1
          ncmt(boxb,imoltb) = orgbib - 1
          if (lgrand) then
             if (boxa.ne.1) then
                parall(imolta,temtyp(imolta))=nchain
                call gcmc_exchange(imolb,orgaib+1)
                imola=imolb
             else
                parall(imoltb,temtyp(imoltb))=nchain
                call gcmc_exchange(imola,orgbia+1)
                imolb=imola
             end if
          end if
       end if

       if ( lewald ) then
          ! update reciprocal-space sum
          if (.not.lideal(boxa)) call recip(boxa,vdum,vdum,2)
          if (boxa.ne.boxb.and..not.lideal(boxb)) call recip(boxb,vdum,vdum,2)
       end if

       ! update center of mass
       if (.not.(lgrand.and.boxa.ne.1)) call ctrmas(.false.,boxa,imolb,8,boxb)
       if (.not.(lgrand.and.boxb.ne.1)) call ctrmas(.false.,boxb,imola,8,boxa)
    else
       if (lgrand) then
          if (boxa.ne.1) then
             call gcmc_cleanup(imolta,boxa)
          else if (boxb.ne.1) then
             call gcmc_cleanup(imoltb,boxb)
          end if
       end if
       if (lewald.and.boxa.eq.boxb.and..not.lideal(boxa)) then
          ! recover the reciprocal space vectors
          ! if the move is not accepted
          call recip(boxa,vdum,vdum,4)
       end if
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end SWATCH in ',myid
#endif

    return
  end subroutine swatch

!> \brief Designed with swatching rigid planar PAH molecules in mind
!> it will work for other rigid stuff, but it's not tailored for that
!>
!> 1) Translates so swatch bead 1 of the swathced molecule is in
!> same position as swatch bead 1 of the other molecule \n
!> 2) Rotate the swathched molecule so its swatch bead1--bead2 vector
!> is aligned with the swatch bead1--bead2 of the other molecule \n
!> 3) Rotate the swatched molecule again so its swatch bead2--bead3
!> vector is parallel with the the bead2--bead3 of the other molecule
!>
!> RESULTS: \n
!> A) The planes defined by the three swatched beads in each molecule
!> will become coplanar.  Thus, if both molecules are planar, then they
!> will definitely become coplanar. \n
!> B) If both molecules have the same 1--2--3 "bond" angles and
!> 1--2 "bond" lengths, the 1--2--3 bonds will be exactly aligned
!>     (where 1,2,3 refer to the three beads you swatched). If not,
!> best of luck; I hope you know what you are trying to do :)
!>
!> \b a subscripts refer to molecule that exists, the "other" \n
!> \b b subscripts refer to molecule that being swatched in, the "self"
!>
!> \author Written by Jake L. Rafferty on the fine day of 2.28.07
  subroutine align_planes(isplist,imol_b,imol_a,itype_b,itype_a,xb,yb,zb)
    use const_math,only:onepi,twopi

    ! INPUT VARIABLES ---
    integer::isplist, imol_b, imol_a, itype_b, itype_a

    ! OUTPUT VARIABLES ---
    real::xb(numax), yb(numax), zb(numax)

    ! LOCAL VARIABLES ---
    integer::i, imoltype_b, nunit_b
    integer::ia_bead1, ia_bead2, ia_bead3
    integer::ib_bead1, ib_bead2, ib_bead3

    real::xa(3),   ya(3),   za(3)
    real::xorigin, yorigin, zorigin
    real::xcross,  ycross,  zcross
    real::xtemp,   ytemp,   ztemp
    real::dxa,     dya
    real::dxb,     dyb

    real::rnorm, d
    real::gamma_a, gamma_b
    real::theta, cos_theta, sin_theta

#ifdef __DEBUG_VERBOSE__
    write(io_output,*) 'BEGIN align_planes in ',myid
#endif

    imoltype_b = moltyp(imol_b)
    nunit_b = nunit(imoltype_b)

    ! Find the three beads on each that are being swatched
    ia_bead1=splist(isplist,1,itype_a)
    ia_bead2=splist(isplist,2,itype_a)
    ia_bead3=splist(isplist,3,itype_a)
    ib_bead1=splist(isplist,1,itype_b)
    ib_bead2=splist(isplist,2,itype_b)
    ib_bead3=splist(isplist,3,itype_b)

    !CC -- FIRST ROTATION

    !CC -- STEP 1, translate bead 1 of both molecules to origin

    ! save position of bead_1 on molecule_a
    xorigin = rxu(imol_a,ia_bead1)
    yorigin = ryu(imol_a,ia_bead1)
    zorigin = rzu(imol_a,ia_bead1)
    ! use this bead as the origin
    xa(1) = 0.0E0_dp
    ya(1) = 0.0E0_dp
    za(1) = 0.0E0_dp
    xa(2) = rxu(imol_a,ia_bead2) - xorigin
    ya(2) = ryu(imol_a,ia_bead2) - yorigin
    za(2) = rzu(imol_a,ia_bead2) - zorigin
    xa(3) = rxu(imol_a,ia_bead3) - xorigin
    ya(3) = ryu(imol_a,ia_bead3) - yorigin
    za(3) = rzu(imol_a,ia_bead3) - zorigin

    ! translate molecule be to the origin
    do i=1,nunit_b
       xb(i) = rxu(imol_b,i) - rxu(imol_b,ib_bead1)
       yb(i) = ryu(imol_b,i) - ryu(imol_b,ib_bead1)
       zb(i) = rzu(imol_b,i) - rzu(imol_b,ib_bead1)
    end do

    ! Get first rotation vector -- the vector orthogonal to
    ! both of the 1--2 vectors, i.e. the cross product

    ! take cross product of 1--2 vectors
    xcross = ya(2)*zb(ib_bead2) - za(2)*yb(ib_bead2)
    ycross = za(2)*xb(ib_bead2) - xa(2)*zb(ib_bead2)
    zcross = xa(2)*yb(ib_bead2) - ya(2)*xb(ib_bead2)
    ! normalize cross product to unit length
    rnorm = sqrt(xcross*xcross + ycross*ycross + zcross*zcross)
    xcross = xcross/rnorm
    ycross = ycross/rnorm
    zcross = zcross/rnorm
    ! find projection of this unit vector onto yz plane
    d = sqrt(ycross*ycross + zcross*zcross)

    !CC -- STEP 2, rotate space about x-axis so that rotation lies in yz plane

    do i=2,3
       ytemp = (ya(i)*zcross - za(i)*ycross)/d
       ztemp = (ya(i)*ycross + za(i)*zcross)/d
       ya(i) = ytemp
       za(i) = ztemp
    end do

    do i=1,nunit_b
       ytemp = (yb(i)*zcross - zb(i)*ycross)/d
       ztemp = (yb(i)*ycross + zb(i)*zcross)/d
       yb(i) = ytemp
       zb(i) = ztemp
    end do

    !CC -- STEP 3, rotate space so that rotation axis is along positive z-axis

    do i=2,3
       xtemp = xa(i)*d - za(i)*xcross
       ztemp = xa(i)*xcross + za(i)*d
       xa(i) = xtemp
       za(i) = ztemp
    end do

    do i=1,nunit_b
       xtemp = xb(i)*d - zb(i)*xcross
       ztemp = xb(i)*xcross + zb(i)*d
       xb(i) = xtemp
       zb(i) = ztemp
    end do

    !CC -- STEP 4, Rotate by theta that about z-axis

    ! first we need the angle between the two 1--2 vectors

    ! angle of vector for molecule_a with x axis
    gamma_a = atan(abs(ya(2)/xa(2)))
    if (ya(2).lt.0.0) then
       if (xa(2).lt.0.0) then
          gamma_a = gamma_a + onepi
       else
          gamma_a = twopi - gamma_a
       end if
    else
       if (xa(2).lt.0.0) gamma_a = onepi - gamma_a
    end if

    ! angle of vector for molecule_b with x axis
    gamma_b = atan(abs(yb(ib_bead2)/xb(ib_bead2)))
    if (yb(ib_bead2).lt.0.0) then
       if (xb(ib_bead2).lt.0.0) then
          gamma_b = gamma_b + onepi
       else
          gamma_b = twopi - gamma_b
       end if
    else
       if (xb(ib_bead2).lt.0.0) gamma_b = onepi - gamma_b
    end if

    theta = gamma_a - gamma_b
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    ! now rotate molecule_b by theta
    do i=1,nunit_b
       xtemp = xb(i)*cos_theta - yb(i)*sin_theta
       ytemp = xb(i)*sin_theta + yb(i)*cos_theta
       xb(i) = xtemp
       yb(i) = ytemp
    end do

    !CC -- INVERT STEP 3

    do i=2,3
       xtemp = xa(i)*d + za(i)*xcross
       ztemp = -xa(i)*xcross + za(i)*d
       xa(i) = xtemp
       za(i) = ztemp
    end do

    do i=1,nunit_b
       xtemp = xb(i)*d + zb(i)*xcross
       ztemp = -xb(i)*xcross + zb(i)*d
       xb(i) = xtemp
       zb(i) = ztemp
    end do

    !CC --- INVERT STEP 2

    do i=2,3
       ytemp = (ya(i)*zcross + za(i)*ycross)/d
       ztemp = (-ya(i)*ycross + za(i)*zcross)/d
       ya(i) = ytemp
       za(i) = ztemp
    end do

    do i=1,nunit_b
       ytemp = (yb(i)*zcross + zb(i)*ycross)/d
       ztemp = (-yb(i)*ycross + zb(i)*zcross)/d
       yb(i) = ytemp
       zb(i) = ztemp
    end do

    !CC --- NOW PERFORM THE SECOND ROTATION
    !CC --- THIS ONE IS ABOUT THE 1--2 VECTOR OF MOLECULE_a

    !CC --- STEP 1, already done in last rotation

    !CC --- STEP 2, rotate space about x-axis so that rotation lies in yz plane

    ! The unit rotation axis
    rnorm = sqrt( xa(2)*xa(2) + ya(2)*ya(2) + za(2)*za(2) )
    xcross = xa(2)/rnorm
    ycross = ya(2)/rnorm
    zcross = za(2)/rnorm
    ! find projection of this unit vector onto yz plane
    d = sqrt(ycross*ycross + zcross*zcross)

    do i=2,3
       ytemp = (ya(i)*zcross - za(i)*ycross)/d
       ztemp = (ya(i)*ycross + za(i)*zcross)/d
       ya(i) = ytemp
       za(i) = ztemp
    end do

    do i=1,nunit_b
       ytemp = (yb(i)*zcross - zb(i)*ycross)/d
       ztemp = (yb(i)*ycross + zb(i)*zcross)/d
       yb(i) = ytemp
       zb(i) = ztemp
    end do

    !CC -- STEP 3, rotate space so that rotation axis is along positive z-axis

    do i=2,3
       xtemp = xa(i)*d - za(i)*xcross
       ztemp = xa(i)*xcross + za(i)*d
       xa(i) = xtemp
       za(i) = ztemp
    end do

    do i=1,nunit_b
       xtemp = xb(i)*d - zb(i)*xcross
       ztemp = xb(i)*xcross + zb(i)*d
       xb(i) = xtemp
       zb(i) = ztemp
    end do

    !CC -- STEP 4, Rotate by theta that about z-axis

    ! This time theta is angle between the two planes
    ! which correponds to the 2--3 vectors

    dxa = xa(3) - xa(2)
    dya = ya(3) - ya(2)

    dxb = xb(ib_bead3) - xb(ib_bead2)
    dyb = yb(ib_bead3) - yb(ib_bead2)
    ! angle of vector for molecule_a with x axis
    gamma_a = atan(abs(dya/dxa))
    if (dya.lt.0.0) then
       if (dxa.lt.0.0) then
          gamma_a = gamma_a + onepi
       else
          gamma_a = twopi - gamma_a
       end if
    else
       if (dxa.lt.0.0) gamma_a = onepi - gamma_a
    end if

    ! angle of vector for molecule_b with x axis
    gamma_b = atan(abs(dyb/dxb))
    if (dyb.lt.0.0) then
       if (dxb.lt.0.0) then
          gamma_b = gamma_b + onepi
       else
          gamma_b = twopi - gamma_b
       end if
    else
       if (dxb.lt.0.0) gamma_b = onepi - gamma_b
    end if

    theta = gamma_a - gamma_b
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    ! now rotate molecule_b by theta
    do i=1,nunit_b
       xtemp = xb(i)*cos_theta - yb(i)*sin_theta
       ytemp = xb(i)*sin_theta + yb(i)*cos_theta
       xb(i) = xtemp
       yb(i) = ytemp
    end do

    !CC -- INVERT STEP 3

    do i=1,nunit_b
       xtemp = xb(i)*d + zb(i)*xcross
       ztemp = -xb(i)*xcross + zb(i)*d
       xb(i) = xtemp
       zb(i) = ztemp
    end do

    !CC --- INVERT STEP 2

    do i=1,nunit_b
       ytemp = (yb(i)*zcross + zb(i)*ycross)/d
       ztemp = (-yb(i)*ycross + zb(i)*zcross)/d
       yb(i) = ytemp
       zb(i) = ztemp
    end do

    !CC --- INVERT STEP 1

    do i=1,nunit_b
       xb(i) = xb(i) + xorigin
       yb(i) = yb(i) + yorigin
       zb(i) = zb(i) + zorigin
    end do

    !CC --- FINALLY WE ARE DONE

#ifdef __DEBUG_VERBOSE__
! Write a check
    ! open(90, file='align.xyz',status='unknown')
    ! nunit_a = nunit(moltyp(imol_a))
    ! write(90,*) nunit_a+nunit_b
    ! write(90,*) 'Plane Alignment Test'
    ! do i=1,nunit_a
    !    write(90,*) 'C ',rxu(imol_a,i),ryu(imol_a,i),rzu(imol_a,i)
    ! end do
    ! do i=1,nunit_b
    !    write(90,*) 'O ',xb(i),yb(i),zb(i)
    ! end do
    ! close(90)
    write(io_output,*) 'END align_planes in ',myid
#endif
  end subroutine align_planes

  subroutine init_swatch(io_input,lprint)
    use var_type,only:default_string_length
    use util_string,only:uppercase
    use util_files,only:readLine
    use util_mp,only:mp_bcast
    integer,intent(in)::io_input
    LOGICAL,INTENT(IN)::lprint
    character(LEN=default_string_length)::line_in
    character(LEN=5)::trash_string ! Paul -- for SAFE-swatch read line
    integer::jerr,i,j,k
    namelist /mc_swatch/ pmswat,nswaty,pmsatc

    if (allocated(bnswat)) deallocate(bnswat,bnswat_empty,bsswat,stat=jerr)
    allocate(bnswat(npamax,npabmax,nbxmax),bnswat_empty(npamax,npabmax,nbxmax),bsswat(npamax,npabmax,nbxmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_swatch: allocation failed',jerr)

    bnswat = 0
    bsswat = 0
    bnswat_empty = 0
    ! liswatch = .false.

    !> read namelist mc_swatch
    nswaty=nmolty*(nmolty-1)/2
    do i=1,nswaty
       pmsatc(i)=real(i,dp)/nswaty
    end do

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=mc_swatch,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_swatch',jerr)
    end if

    call mp_bcast(pmswat,1,rootid,groupid)
    call mp_bcast(nswaty,1,rootid,groupid)
    call mp_bcast(pmsatc,nswaty,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_SWATCH','------------------------------------------'
       write(io_output,'(A,G16.9)') 'pmswat: ',pmswat
       write(io_output,'(A,I0)') '   number of swatch pairs (nswaty): ',nswaty
       do i=1,nswaty
          write(io_output,'(A,G16.9)') '   probability of each swatch pair: ',pmsatc(i)
       end do
    end if

    if (nswaty.gt.npamax) call err_exit(__FILE__,__LINE__,'nswaty gt npamax',myid+1)
    if (pmswat.gt.0.and.lneigh) call err_exit(__FILE__,__LINE__,'Neighbor list currently does not work with CBMC particle identity switch moves!',myid+1)

    ! Looking for section MC_SWATCH
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_SWATCH:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section MC_SWATCH not found',jerr)

          if (UPPERCASE(line_in(1:9)).eq.'MC_SWATCH') then
             do i=1,nswaty+1
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                if (UPPERCASE(line_in(1:13)).eq.'END MC_SWATCH') then
                   if (i.ne.nswaty+1) call err_exit(__FILE__,__LINE__,'Section MC_SWATCH not complete!',myid+1)
                   exit
                else if (i.eq.nswaty+1) then
                   call err_exit(__FILE__,__LINE__,'Section MC_SWATCH has more than nswaty records!',myid+1)
                end if

                ! moltyp1<->moltyp2 nsampos 2xncut
                ! Paul -- SAFE-swatch part
                if (scan(line_in(1:5),'sS') >= 1) then
                    read(line_in,*) trash_string,nswatb(i,1:2),nsampos(i),ncut(i,1:2),ncutsafe(i,1:2)
                    if (ncut(i,1) .eq. 0 .and. ncutsafe(i,1) .eq. 0) then
                        call err_exit(__FILE__,__LINE__,'ncut and ncutsafe cannot be zero at the same time',myid+1)
                    end if
                else
                    read(line_in,*) nswatb(i,1:2),nsampos(i),ncut(i,1:2)
                    ncutsafe(i,1:2) = 0
                end if
                if (nswatb(i,1).eq.nswatb(i,2)) then
                   ! safety checks on swatch
                   write(io_output,*) 'nswaty ',i,' has identical moltyp'
                   call err_exit(__FILE__,__LINE__,'cannot swatch identical moltyp',myid+1)
                end if

                if (lprint) then
                   write(io_output,'(/,A,2(4X,I0))') '   swatch molecule type pairs:',nswatb(i,1:2)
                   write(io_output,'(A,I0,A,2(2X,I0),A,2(2X,I0))') '   nsampos: ',nsampos(i),', ncut:',ncut(i,1:2),', ncutsafe:',ncutsafe(i,1:2)
                end if

                ! Read gswatc information
                if (ncut(i,1) .gt. 0) then
                    call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                    ! gswatc 2x(ifrom, iprev)
                    read(line_in,*) (gswatc(i,j,1:2*ncut(i,j)),j=1,2)
                end if

                if (ncutsafe(i,1) .gt. 0) then
                    call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                    ! gswatcsafe (ifrom, iprev, index)
                    read(line_in,*) (gswatcsafe(i,j,1:3*ncutsafe(i,j)),j=1,2)
                end if

                if (lprint) then
                   do j=1,2
                      write(io_output,FMT='(A,I0)') '   molecule ',j
                      do k = 1,ncut(i,j)
                         write(io_output,'(3(A,I0))') '   ncut ',k,': grow from ',gswatc(i,j,2*k-1),', prev ',gswatc(i,j,2*k)
                      end do
                      do k = 1, ncutsafe(i, j)
                         write(io_output,'(4(A,I0))') '   ncutsafe ',k,': grow from ',gswatcsafe(i,j,3*k-2),', prev ',gswatcsafe(i,j,3*k-1),', index ',gswatcsafe(i,j,3*k)
                      end do
                   end do
                end if

                do j = 1,nsampos(i)
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                   ! splist
                   read(line_in,*) splist(i,j,1:2)
                   if (lprint) then
                      write(io_output,'(A,2(4X,I0))') '   splist:',splist(i,j,1:2)
                   end if
                end do

                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                ! nswtcb pmswtcb
                read(line_in,*) nswtcb(i),pmswtcb(i,1:nswtcb(i))
                if (lprint) then
                   write(io_output,'(A,I0)') '   number of swatch box pairs: ',nswtcb(i)
                   do j=1,nswtcb(i)
                      write(io_output,'(A,G16.9)') '   probability of the swatch box pair: ',pmswtcb(i,j)
                   end do
                end if

                do j = 1,nswtcb(i)
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                   ! box numbers
                   read(line_in,*) box3(i,j),box4(i,j)
                   if (lprint) then
                      write(io_output,'(A,2(4X,I0))') '   box pair:',box3(i,j),box4(i,j)
                   end if
                   if (pmswat.gt.0.and.licell.and.(box3(i,j).eq.boxlink.or.box4(i,j).eq.boxlink)) call err_exit(__FILE__,__LINE__,'Cell structure currently does not work with CBMC particle identity switch moves!',myid+1)
                end do
             end do
             exit cycle_read_swatch
          end if
       END DO CYCLE_READ_SWATCH
    end if

    call mp_bcast(nswatb,npamax*2,rootid,groupid)
    call mp_bcast(nsampos,nswaty,rootid,groupid)
    call mp_bcast(ncut,npamax*2,rootid,groupid)
    call mp_bcast(ncutsafe,npamax*2,rootid,groupid)
    call mp_bcast(gswatc,npamax*npamax*4,rootid,groupid)
    call mp_bcast(gswatcsafe,npamax*npamax*6,rootid,groupid)
    call mp_bcast(splist,npamax*numax*2,rootid,groupid)
    call mp_bcast(nswtcb,nswaty,rootid,groupid)
    call mp_bcast(pmswtcb,npamax*npabmax,rootid,groupid)
    call mp_bcast(box3,npamax*npabmax,rootid,groupid)
    call mp_bcast(box4,npamax*npabmax,rootid,groupid)
  end subroutine init_swatch

  subroutine output_swatch_stats(io_output)
    integer,intent(in)::io_output
    integer::i,j,ibox,jbox,ii

    write(io_output,'(/,A,/)') '### Molecule swatch     ###'
    do i = 1, nswaty
       if (sum(bnswat(i,:,:)) .gt. 0) then ! only output if swatch attempt > 0
          write(io_output,'(A,I0)') 'pair typ = ',i
          write(io_output,'(A,I0,A)',advance='no') 'moltyps = ',nswatb(i,1),'    '
          write(io_output,'(A10)',advance='no')molecname(nswatb(i,1))
          write(io_output,'(A,I0,A)',advance='no')' and ',nswatb(i,2),'    '
          write(io_output,'(A10)',advance='no')molecname(nswatb(i,2))
          write(io_output,*)
          do j = 1, nswtcb(i)
             do jbox = 1,2
                if (jbox.eq.1) ibox=box3(i,j)
                if (jbox.eq.2) then
                   if (box3(i,j).eq.box4(i,j)) exit
                   ibox=box4(i,j)
                end if
                ! JLR 12-1-09 changing to exclude empty box attempts from swatch rate
                write(io_output,"('between box ',I0,' and ',I0,' into box ',I0,'   uattempts = ',I0,   '  attempts = ',I0,'  accepted = ',I0)") box3(i,j),box4(i,j),ibox,bnswat(i,j,ibox),bnswat(i,j,ibox)-bnswat_empty(i,j,ibox),bsswat(i,j,ibox)
                if (bnswat(i,j,ibox) .gt. 0) then
                   write(io_output,"(' accepted % =',F7.3)") 100.0_dp*real(bsswat(i,j,ibox),dp)/real(bnswat(i,j,ibox)-bnswat_empty(i,j,ibox),dp)
                end if
                ! EN JLR 12-1-09
             end do
          end do
       end if
    end do
  end subroutine output_swatch_stats

  subroutine read_checkpoint_swatch(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) bnswat,bnswat_empty,bsswat
    call mp_bcast(bnswat,npamax*npabmax*nbxmax,rootid,groupid)
    call mp_bcast(bnswat_empty,npamax*npabmax*nbxmax,rootid,groupid)
    call mp_bcast(bsswat,npamax*npabmax*nbxmax,rootid,groupid)
  end subroutine read_checkpoint_swatch

  subroutine write_checkpoint_swatch(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) bnswat,bnswat_empty,bsswat
  end subroutine write_checkpoint_swatch
end module transfer_swatch
