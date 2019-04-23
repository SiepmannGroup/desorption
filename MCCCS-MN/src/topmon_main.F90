MODULE topmon_main
  use sim_system
  use sim_cell
  implicit none
  private
  save
  public::output_version,monola

  ! variables added for GCMC histogram reweighting
  integer,parameter::fmax=1E6,nprop1=11
  logical::lstop=.false.,use_checkpoint=.false.
  integer::blockm,checkpoint_interval=1800,checkpoint_copies=1,nstep=1,nnstep,nnn,acmove,acnp,acipsw,nblock&
   ,io_movie,io_solute,io_cell,io_traj,time_limit=0&
   ,N_add=0,box2add=1,moltyp2add=1 !< defaults for mc_shared
  real::enthalpy,enthalpy2& !< enthalpy (NpT) or internal energy (NVT)
     ,acdvdl,binvir,binvir2
  real,allocatable::acdens(:,:)& !< (ibox,itype): accumulators of box density
   ,molfra(:,:)& !< (ibox,itype): accumulators of mole fraction
   ,acnbox(:,:)& !< accumulators of ncmt
   ,acnbox2(:,:,:)& !< accumulators of ncmt2
   ,acv(:,:),acvsq(:,:),acvkjmol(:,:)& !< (j,ibox): accumulators of vbox
   ,acdipole(:,:),acdipolesq(:,:)& !< (j,ibox): accumulators of dipole values; j = 1: dipolex; 2: dipoley; 3: dipolez; 4: norm of dipole
   ,acboxa(:,:)& !< (ibox,j) accumulators of cell angles for ibox
   ,acboxl(:,:)& !< (ibox,j) accumulators of cell lengths for ibox
   ,acvol(:),acvolsq(:),acvolume(:)& !< accumulators of box volume
   ,acpres(:),acsurf(:),accomp(:)& !< accumulators of pressure,surface tension and compressibility factor
   ,acEnthalpy(:),acEnthalpy1(:)& !< accumulators of enthalpies using calculated pressure (acEnthalpy) and specified enthalpy (acEnthalpy1)
   ,acsolpar(:,:,:),avsolinter(:,:),avsolintra(:,:),avsolbend(:,:),avsoltor(:,:),avsolelc(:,:)&
   ,asetel(:,:)& !< (ibox,itype): accumulators of square end-to-end distance; counter is in mnbox
   ,vstart(:),vend(:),pres(:),compress(:),molvol(:),speden(:)&
   ,aver1(:,:,:),stdev1(:,:,:),errme1(:,:,:),bccold1(:,:,:),baver1(:,:,:,:)& !< accumulators for solubility parameter and heat of vaporization (Neeraj)
   ,aver(:,:),stdev(:,:),errme(:,:),bccold(:,:),baver(:,:,:) !< (j,ibox): properties of ibox; j =
  !< \verbatim
  !< ---------------------------------------------------------
  !< 1                                              = specific density
  !< 2                                              = pressure
  !< 3                      to (2+nEnergy)          = energies
  !< 1+(2+nEnergy)          to (2+nEnergy)+  nmolty = chemical potential
  !< 1+(2+nEnergy)+  nmolty to (2+nEnergy)+2*nmolty = square-end-to-end-length
  !< 1+(2+nEnergy)+2*nmolty to (2+nEnergy)+3*nmolty = number density
  !< 1+(2+nEnergy)+3*nmolty to (2+nEnergy)+4*nmolty = mole fraction
  !< 1+(2+nEnergy)+4*nmolty                         = surface tension
  !< 2+(2+nEnergy)+4*nmolty                         = system volume
  !< 3+(2+nEnergy)+4*nmolty                         = enthalpy inst
  !< 4+(2+nEnergy)+4*nmolty                         = enthalpy ext
  !< ---------------------------------------------------------
  !< \endverbatim
  integer,allocatable::io_box_movie(:),io_box_movie_pdb(:),nminp(:),nmaxp(:),ncmt_list(:,:)&
   ,mnbox(:,:),solcount(:,:),nccold1(:,:,:),nccold(:,:),ndist(:,:) !< GCMC reweighting histograms

  ! for flexible filenames
  character(LEN=20)::string

contains
#include "defines.h"
  !> \brief Output version and build information
  !>
  !> defines.h is generated automatically by CMake using
  !> template defines.h.in. If compiling manually, create
  !> defines.h and fill in relevant parameters.
  subroutine output_version(io_version)
    integer,intent(in)::io_version

    write(io_version,'(/,3A,/,3(2A,/),4A,/)')&
     'MCCCS topmon (branch: ',__BRANCH__,')',&
     'Commit hash: ',__COMMIT_HASH__,&
     'Build on host: ',__BUILD_HOSTNAME__,&
     'Preprocessor definitions: ',__PREPROCESSOR_DEFINITIONS__,&
     'Using ',__Fortran_COMPILER_ID__,' compiler: ',__Fortran_COMPILER__
  end subroutine output_version


!> \brief Main control logic of topmon
!>
!> reads the control-data from unit 4
!> starts and controls the simulation
  subroutine monola(file_in)
    use var_type,only:dp,default_string_length
    use const_math,only:twopi,raddeg
    use const_phys,only:debroglie_factor,N_Avogadro,R_gas,MPa2SimUnits
    use util_math,only:update_average,calculate_statistics
    use util_random,only:random
    use util_string,only:format_n,integer_to_string
    use util_runtime,only:err_exit
    use util_timings,only:time_init,time_date_str,time_now
    use util_files,only:get_iounit
    use util_kdtree,only:allocate_kdtree,construct_kdtree
    use sim_initia,only:setup_molecule_config
    use sim_particle,only:ctrmas,neighbor,neigh_cnt
    use energy_pairwise,only:sumup
    use moves_simple,only:translation,rotation,Atom_translation,output_translation_rotation_stats
    use moves_volume,only:volume_1box,volume_2box,output_volume_stats
    use moves_cbmc,only:config,schedule,output_safecbmc,output_cbmc_stats
    use moves_ee,only:eesetup,eemove,ee_index_swap,expand,numcoeff,output_ee_stats
    use transfer_shared,only:gcmc_setup,gcmc_exchange
    use transfer_swap,only:swap,cnt,output_swap_stats,acchem,bnchem
    use transfer_swatch,only:swatch,output_swatch_stats
    use prop_pressure,only:pressure
    use prop_widom,only:write_deltaG_map,write_prop_widom,write_prop_widom_with_stats

    character(LEN=*),intent(in)::file_in

    character(LEN=default_path_length)::file_flt,file_hist,file_ndis,file_cnt,file_config

    ! descriptions of different kinds of energies
    character(LEN=default_string_length)::vname(nEnergy)=(/' Total energy',' Inter LJ    ',' Tail  LJ    ',' Intra LJ    '&
     ,' Stretch     ',' Bond bending',' Torsion     ',' Coulomb     ',' External pot',' 3-body Garo ',' Fluc Q      '&
     ,'             ','             ','             '/)

    integer::io_flt,io_hist,io_cnt,io_ndis,io_config,i,jerr,ibox,itype,itype2,Temp_nmol,nentry,j,nummol,imolty,ii,ntii&
     ,igrow,steps,itemp,jbox,itel,ig,il,nbl,n,zzz,ichkpt,nnn_1st,nstep_per_cycle
    real::v(nEnergy),press1,surf,comp,time_prev,time_cur,rm,temvol,tmp,vhist,eng_list(fmax),temacd,temspd,debroglie,starviro,dummy&
     ,inside,bvirial,gconst,ostwald,stdost,molfrac,time_st,time_av,time_cycle=0.0E0_dp,nImages
    logical::ovrlap
! ----------------------------------------------------------------
    ! Initialize the timer
    call time_init()

    call readdat(file_in)

    if (allocated(nminp)) deallocate(nminp,nmaxp,ncmt_list,ndist,vstart,vend,pres,molvol,speden,acdens,molfra,acnbox,acnbox2,acv,acvsq,acvkjmol,acdipole,acdipolesq,acboxa,acboxl,acvol,acvolsq,acvolume,acpres,acsurf,acEnthalpy,acEnthalpy1,solcount,acsolpar,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc,mnbox,asetel,nccold1,bccold1,baver1,nccold,bccold,baver,aver1,stdev1,errme1,aver,stdev,errme,compress,accomp,stat=jerr)
    allocate(nminp(ntmax),nmaxp(ntmax),ncmt_list(fmax,ntmax),ndist(0:nmax,ntmax),vstart(nbxmax),vend(nbxmax),pres(nbxmax)&
     ,compress(nbxmax),molvol(nbxmax),speden(nbxmax),acdens(nbxmax,ntmax),molfra(nbxmax,ntmax),acnbox(nbxmax,ntmax),acnbox2(nbxmax,ntmax,20)&
     ,acv(nEnergy,nbxmax),acvsq(nEnergy,nbxmax),acvkjmol(nEnergy,nbxmax),acdipole(4,nbxmax),acdipolesq(4,nbxmax),acboxa(nbxmax,3)&
     ,acboxl(nbxmax,3),acvol(nbxmax),acvolsq(nbxmax),acvolume(nbxmax),acpres(nbxmax),acsurf(nbxmax),accomp(nbxmax),acEnthalpy(nbxmax)&
     ,acEnthalpy1(nbxmax),solcount(nbxmax,ntmax),acsolpar(nprop1,nbxmax,nbxmax),avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax)&
     ,avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax),avsolelc(nbxmax,ntmax),mnbox(nbxmax,ntmax),asetel(nbxmax,ntmax)&
     ,nccold1(nprop1,nbxmax,nbxmax),bccold1(nprop1,nbxmax,nbxmax),baver1(nprop1,nbxmax,nbxmax,blockm),nccold(nprop,nbxmax)&
     ,bccold(nprop,nbxmax),baver(nprop,nbxmax,blockm),aver1(nprop1,nbxmax,nbxmax),stdev1(nprop1,nbxmax,nbxmax)&
     ,errme1(nprop1,nbxmax,nbxmax),aver(nprop,nbxmax),stdev(nprop,nbxmax),errme(nprop,nbxmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'monola: allocation failed',jerr)
    end if

    ! SETTING UP ARRAYS FOR ANALYSYS PURPOSE
    ! JLR 11-11-09
    ! do not call analysis if you set ianalyze to be greater than number of cycles
    ! KM 01/10 remove analysis
    ! if (ianalyze.le.nstep) then
    !    call analysis(0)
    ! end if
    ! END JLR 11-11-09

    ! set the centers of mass if LFOLD = .TRUE.
    if (lfold) then
       do ibox=1,nbox
          call ctrmas(.true.,ibox,0,6)
       end do
    end if

    !> Set up the kd-tree
    if (lkdtree) then
        do ibox = 1, nbox
            if (lkdtree_box(ibox)) call construct_kdtree(ibox, ibox, .true.)
        end do
    end if

    if (licell) then
       ! check that rintramax is really valid
       do i=1,nchain
          if (2.0_dp*rcmu(i).gt.rintramax) call err_exit(__FILE__,__LINE__,'rintramax for linkcell list is too small',myid+1)
       end do

       ! set up initial linkcell
       call build_linked_cell()
    end if

    ! set up thermodynamic integration stuff
    if (lmipsw) then
       call ipswsetup()
    else
       lstagea = .false.
       lstageb = .false.
       lstagec = .false.
    end if

    ! set up expanded ensemble stuff
    if (lexpee) then
       call eesetup
       if (lmipsw) call err_exit(__FILE__,__LINE__,'not for BOTH lexpee AND lmipsw',myid+1)
    else
       leemove=.false.
    end if

    if (lneighbor) then
       neighbor = 0
       if (myid.eq.rootid) then
          file_cnt='fort.21'
          io_cnt=get_iounit()
          open(unit=io_cnt,access='sequential',action='write',file=file_cnt,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open cnt file '//trim(file_cnt),jerr)
          write(io_cnt,*) 'ii:',nnstep,(neigh_cnt(i),i=1,nchain)
       end if
    end if


    ! setup files for histogram reweighting
    if(lgrand) then
       if (myid.eq.rootid) then
          file_flt = "nfl"//run_num(1:len_trim(run_num))//suffix//".dat"
          !write(file_flt,'("nfl",I1.1,A,".dat")') run_num,suffix
          io_flt=get_iounit()
          open(unit=io_flt,access='sequential',action='write',file=file_flt,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open flt file '//trim(file_flt),jerr)

          file_hist = "his"//run_num(1:len_trim(run_num))//suffix//".dat"
          !write(file_hist,'("his",I1.1,A,".dat")') run_num,suffix
          io_hist=get_iounit()
          open(unit=io_hist,access='sequential',action='write',file=file_hist,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open hist file '//trim(file_hist),jerr)
          write(io_hist,'(f8.4,2x,i0,2x,'//format_n(nmolty,'g15.5')//',3(2x,f12.3))')  temp,nmolty,(temp*log(B(i)),i=1,nmolty),boxlx(1),boxly(1),boxlz(1)
       end if

       ! extra zero accumulator for grand-canonical ensemble
       nentry = 0
       nminp = 1E6_dp
       nmaxp = -1E6_dp
       ndist = 0
    end if

    ! 2nd viral coefficient
    if (lvirial) then
       binvir = 0.0E0_dp
       binvir2 = 0.0E0_dp
       ! if ( lvirial2 ) then
       !    call virial2(binvir,binvir2,nvirial,starvir,stepvir)
       !    goto 2000
       ! end if
       ! profile=0.0E0_dp
       ! binstep = 0.05E0_dp
       ! lvirial2 = .false.
    end if

    ! initialize variables
    acmove = 0
    acdens = 0.0E0_dp
    molfra = 0.0E0_dp
    acnbox = 0.0E0_dp
    acnbox2 = 0.0E0_dp
    acv = 0.0E0_dp
    acvsq = 0.0E0_dp
    acvkjmol = 0.0E0_dp
    acdipole = 0.0_dp
    acdipolesq = 0.0_dp
    acboxa = 0.0E0_dp
    acboxl = 0.0E0_dp
    acvol = 0.0E0_dp
    acvolsq = 0.0E0_dp
    acvolume = 0.0E0_dp
    enthalpy = 0.0E0_dp
    enthalpy2 = 0.0E0_dp
    acnp = 0
    acpres = 0.0E0_dp
    acsurf = 0.0E0_dp
    accomp = 0.0E0_dp
    acEnthalpy = 0.0E0_dp
    acEnthalpy1 = 0.0E0_dp
    solcount = 0
    acsolpar=0.0E0_dp
    avsolinter = 0.0E0_dp
    avsolintra = 0.0E0_dp
    avsolbend = 0.0E0_dp
    avsoltor = 0.0E0_dp
    avsolelc = 0.0E0_dp
    mnbox  = 0
    asetel = 0.0E0_dp
    acipsw=0
    acdvdl=0.0E0_dp
    ! accumulators for block averages ---
    nblock = 0
    nccold1 = 0
    bccold1 = 0.0E0_dp
    baver1=0.0E0_dp
    nccold = 0
    bccold = 0.0E0_dp
    baver=0.0E0_dp
    ! accumulators for fluctuating charge performance
    bnflcq = 0.0E0_dp
    bnflcq2 = 0.0E0_dp
    bsflcq = 0.0E0_dp
    bsflcq2 = 0.0E0_dp
    ichkpt = 0

! -----------------------------------------------------------------
    if (use_checkpoint) then
       call read_checkpoint_main('save-stats')
       nnn_1st=nnn+1
    else
       nnn_1st=1
       ! calculate initial energy and check for overlaps ***
       do ibox=1,nbox
          call sumup(ovrlap,v,ibox,lvol=.false.)
          vbox(:,ibox) = v
          vbox(ivIpswb,ibox) = vipsw
          vbox(ivWellIpswb,ibox) = vwellipsw
          if (ovrlap) then
             call err_exit(__FILE__,__LINE__,'overlap in initial configuration',myid+1)
          end if

          vstart(ibox) = vbox(ivTot,ibox)
          if (myid.eq.rootid) then
             write(io_output,*)
             write(io_output,*) 'box  ',ibox,' initial v   = ', vbox(ivTot,ibox)
          end if

          ! calculate initial pressure ***
          call pressure( press1, surf, comp, ibox )
          if (myid.eq.rootid) then
             write(io_output,"(' surf. tension :   box',i2,' =',f14.5)") ibox, surf
             write(io_output,"(' pressure check:   box',i2,' =',f14.2)") ibox, press1
             write(io_output,"(' compress factor:  box',i2,' =',f14.5)") ibox, comp
         end if
       end do

       if (myid.eq.rootid) then
          write(io_output,*)
          write(io_output,*) '+++++ start of markov chain +++++'
          write(io_output,*)
          write(io_output,*)  'Cycle   Total   Energy    Boxlength    Pressure     Compress    Molecules'
       end if
    end if

    if (lstop) then
       nstep_per_cycle = 1
    else
       nstep_per_cycle = nchain
    end if

!************************************************************
! loops over all cycles and all molecules                  **
!************************************************************
    time_prev = time_now()
    do nnn = nnn_1st, nstep
       tmcc = nnstep + nnn
       time_st = time_now()
       do ii = 1, nstep_per_cycle
          acmove = acmove + 1
          ! select a move-type at random ***
          rm = random(-1)

          ! special ensemble dependent moves ###
          if (rm.le.pmvol) then
             ! volume move ---
             if (lnpt) then
                call volume_1box()
             else
                call volume_2box()
             end if
          else if (rm.le.pmswat) then
             ! CBMC switch move ---
             call swatch()
          else if (rm.le.pmswap) then
             ! swap move for linear and branched molecules ---
             call swap()
          else if (rm.le.pmcb) then
             ! configurational bias move ---
             call config()
          else if (rm.le.pmflcq) then
             ! displacement of fluctuating charges ---
             call flucq(2,0)
          else if (rm.le.pmexpc) then
             ! expanded-ensemble move ---
             call expand()
          else if (rm.le.pmexpc1) then
             ! new expanded-ensemble move ---
             ! call expand
             if (random(-1).le.eeratio) then
                call ee_index_swap()
             else
                call eemove()
             end if
          else if (rm.le.pm_atom_tra) then
             call Atom_translation()
          else if (rm.le.pmtra) then
             ! translational move ---
             call translation()
          else
             ! rotation move --
             call rotation()
          end if

          ! accumulate probability of being in an expanded ensemble state
          if (lexpee) then
             ee_prob(mstate) = ee_prob(mstate)+1
          end if

          ! calculate instantaneous values ***
          ! accumulate averages ***
          do ibox=1,nbox
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                temvol = cell_vol(ibox)
             else
                if ( lpbcz ) then
                   temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
                else
                   temvol = boxlx(ibox)*boxly(ibox)
                end if
             end if

             do itype = 1, nmolty
                call update_average(acdens(ibox,itype),ncmt(ibox,itype)/temvol,acmove)
                if (nchbox(ibox).gt.0) then
                   call update_average(molfra(ibox,itype),real(ncmt(ibox,itype),dp)/real(nchbox(ibox),dp),acmove)
                end if

                call update_average(acnbox(ibox,itype),real(ncmt(ibox,itype),dp),acmove)
                if (lexpand(itype)) then
                   do itype2 = 1, numcoeff(itype)
                      call update_average(acnbox2(ibox,itype,itype2),real(ncmt2(ibox,itype,itype2),dp),acmove)
                   end do
                end if
             end do

             call update_average(acv(1:11,ibox),vbox(1:11,ibox),acmove)
             call update_average(acvsq(1:11,ibox),vbox(1:11,ibox)**2,acmove)

             if (lnpt) then
                call update_average(acv(14,ibox),vbox(ivTot,ibox)*temvol,acmove)
             end if

             ! KMB/KEA Energy in kJ/mol
             Temp_nmol = sum(ncmt(ibox,1:nmolty))
             call update_average(acvkjmol(1:11,ibox),vbox(1:11,ibox)/Temp_nmol,acmove)

             ! leftover from Bin, not currently used
             if ( ldielect ) then
                call update_average(acdipole(1,ibox),dipolex(ibox),acmove)
                call update_average(acdipolesq(1,ibox),dipolex(ibox)**2,acmove)
                call update_average(acdipole(2,ibox),dipoley(ibox),acmove)
                call update_average(acdipolesq(2,ibox),dipoley(ibox)**2,acmove)
                call update_average(acdipole(3,ibox),dipolez(ibox),acmove)
                call update_average(acdipolesq(3,ibox),dipolez(ibox)**2,acmove)
                call update_average(acdipole(4,ibox),sqrt(dipolex(ibox)*dipolex(ibox)+dipoley(ibox)*dipoley(ibox)&
                 +dipolez(ibox)*dipolez(ibox)),acmove)
                acdipolesq(4,ibox)=acdipolesq(1,ibox)+acdipolesq(2,ibox)+acdipolesq(3,ibox)
             end if

             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                boxlx(ibox) = cell_length(ibox,1)
                boxly(ibox) = cell_length(ibox,2)
                boxlz(ibox) = cell_length(ibox,3)
                call update_average(acboxa(ibox,1:3),cell_ang(ibox,1:3),acmove)
             end if

             call update_average(acboxl(ibox,1),boxlx(ibox),acmove)
             call update_average(acboxl(ibox,2),boxly(ibox),acmove)
             call update_average(acboxl(ibox,3),boxlz(ibox),acmove)

             call update_average(acvol(ibox),temvol,acmove)
             call update_average(acvolsq(ibox),temvol*temvol,acmove)
             call update_average(acvolume(ibox),temvol,acmove)
          end do

          if (.not.lgibbs) then
             ibox=1
             if (lnpt) then
                tmp=vbox(ivTot,ibox)+express(ibox)*boxlx(ibox)*boxly(ibox)*boxlz(ibox)
             else
                tmp=vbox(ivTot,ibox)
             end if
             call update_average(enthalpy,tmp,acmove)
             call update_average(enthalpy2,tmp*tmp,acmove)
          end if

          ! collect histogram data (added 8/30/99)
          if (lgrand) then
             ibox = 1
             vhist = vbox(ivInterLJ,ibox) + vbox(ivElect,ibox) + vbox(ivFlucq,ibox) !inter+elect+flucq
             if (mod(acmove,ninstf).eq.0.and.myid.eq.rootid) then
                write(io_flt,FMT='(i0,5x,'//format_n(nmolty,'i7')//',5x,g15.6)') acmove,(ncmt(ibox,i),i=1,nmolty),vhist
             end if

             if(mod(acmove,ninsth).eq.0.and.acmove.gt.nequil) then
                do imolty = 1, nmolty
                   nminp(imolty) = min(nminp(imolty),ncmt(ibox,imolty))
                   nmaxp(imolty) = max(nmaxp(imolty),ncmt(ibox,imolty))
                end do
                nentry = nentry + 1

                do imolty=1,nmolty
                   ncmt_list(nentry,imolty) = ncmt(ibox,imolty)
                   ndist(ncmt(ibox,imolty),imolty) = ndist(ncmt(ibox,imolty),imolty) + 1
                end do

                eng_list(nentry) = vhist
             end if

             if (mod(acmove,ndumph).eq.0.and.myid.eq.rootid) then
                do i=1,nentry
                   write(io_hist, * ) (ncmt_list(i,imolty), imolty=1,nmolty), eng_list(i)
                end do
                nentry = 0

                io_ndis=get_iounit()
                do imolty=1,nmolty
                   string = integer_to_string(imolty)
                   file_ndis = "n"//string(1:len_trim(string))//"dis"//run_num(1:len_trim(run_num))//suffix//".dat"
                   !write(file_ndis,'("n",I2.2,"dis",I1.1,A,".dat")') imolty,run_num,suffix
                   open(unit=io_ndis,access='sequential',action='write',file=file_ndis,form='formatted',iostat=jerr&
                    ,status='unknown')
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open ndis file '//trim(file_ndis),jerr)
                   do n=nminp(imolty),nmaxp(imolty)
                      write(io_ndis,*) n,ndist(n,imolty)
                   end do
                end do
                close(io_ndis)
             end if
          end if
       end do
!*********************************************
! ends loop over chains                     **
!*********************************************
       ! perform periodic operations
       call monper()

       time_cur = time_now()
       if (time_cur-time_prev.gt.checkpoint_interval) then
          time_prev = time_cur
          ichkpt = ichkpt+1
          if (ichkpt.gt.checkpoint_copies) ichkpt = 1
          ! write out the restart configurations to save-config file
          if (myid.eq.rootid) then
             call dump('save-config.'//integer_to_string(ichkpt))
             call write_checkpoint_main('save-stats.'//integer_to_string(ichkpt))
          end if
       end if
       if (time_limit.gt.0) then
           time_cycle = time_cycle + (time_cur - time_st)
           time_av = time_cycle / (nnn-nnn_1st+1.0E0_dp)
           if (time_cur.gt.(time_limit-time_av-120.0E0_dp)) exit
       end if
    end do
!********************************************
! ends the loop over cycles                **
!********************************************
    call cnt()

    if (lneighbor.and.myid.eq.rootid) then
       write(io_cnt,*) 'ii:',tmcc,(neigh_cnt(i),i=1,nchain)
       close(io_cnt)
    end if

    call output_safecbmc()

    ! do bin = 1,1000
    !    write(26,*) binstep*(real(bin,dp)-0.5E0_dp),profile(bin)/nstep
    ! end do

    if (myid.eq.rootid) then
       if (io_movie.ge.0) close(io_movie)
       do ibox=1,nbox
          if (io_box_movie(ibox).ge.0) close(io_box_movie(ibox))
          if (io_box_movie_pdb(ibox).ge.0) close(io_box_movie_pdb(ibox))
       end do
       if (io_solute.ge.0) close(io_solute)
       if (io_cell.ge.0) close(io_cell)
       if (ltraj) close(io_traj)

       if (lgrand) then
          close(io_flt)
          close(io_hist)
       end if

! -------------------------------------------------------------------
       ! Output MC move statistics
       write(io_output,*)
       write(io_output,*) '+++++ end of markov chain +++++'

       call output_translation_rotation_stats(io_output)

       if (pmcb .gt. 0.0E0_dp) then
          call output_cbmc_stats(io_output)
       end if

       ! write some information about volume performance ***
       if (lgibbs .or. lnpt) then
          call output_volume_stats(io_output)
       end if

       call output_swap_stats(io_output)
       if (nmolty > 1) then
          call output_swatch_stats(io_output)
       end if

       write(io_output,*)
       write(io_output,*)    '### Charge Fluctuation  ###'
       write(io_output,*)

       do i = 1, nmolty
          do j = 1,nbox
             bnflcq2(i,j) = bnflcq2(i,j) + bnflcq(i,j)
             bsflcq2(i,j) = bsflcq2(i,j) + bsflcq(i,j)
             if (bnflcq2(i,j) .gt. 0.5E0_dp) then
                write(io_output,*) 'molecule typ =',i,'  box =',j
                bsflcq2(i,j) = bsflcq2(i,j)/bnflcq2(i,j)
                write(io_output,"(' attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)") bnflcq2(i,j),bsflcq2(i,j)&
                 ,rmflcq(i,j)
             end if
          end do
       end do

       call output_ee_stats(io_output)

       write(io_output,"(/,'New Biasing Potential')")
       do i=1,nmolty
          write(io_output,"(/,'molecule ',I2,': ',$)") i
          do j=1,nbox
             write(io_output,"(G16.9,1X,$)") eta2(j,i)
          end do
       end do
       write(io_output,*)
    end if
! -------------------------------------------------------------------
    ! Check final energies
    do ibox=1,nbox
       if ( ldielect ) then
          ! store old dipole moment
          call dipole(ibox,2)
          ! sumup the final box dipole
          call dipole(ibox,0)
       end if

       ! checks final value of the potential energy is consistent ***
       call sumup(ovrlap,v,ibox,.false.)
       vend(ibox) = v(ivTot)

       ! need to check
       if (myid.eq.rootid) then
          if ( abs(v(ivTot) - vbox(ivTot,ibox)) .gt. 0.0001) then
             write(io_output,*) '### problem with energy ###  box ',ibox
             write(io_output,*) ' Total energy: ',v(ivTot),vbox(ivTot,ibox),v(ivTot)-vbox(ivTot,ibox)
          end if
          if ( abs(v(ivInterLJ) - vbox(ivInterLJ,ibox)) .gt. 0.000001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Inter mol.en.: ',v(ivInterLJ),vbox(ivInterLJ,ibox),v(ivInterLJ) - vbox(ivInterLJ,ibox)
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                write(io_output,*)'You might check the cutoff wrt box widths'
                write(io_output,*) 'Normal PBC might be failing'
             end if
          end if
          if ( abs(v(ivTail) - vbox(ivTail,ibox)) .gt. 0.000001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Tail corr.en.: ',v(ivTail),vbox(ivTail,ibox),v(ivTail) - vbox(ivTail,ibox)
          end if
          if ( abs(v(ivIntraLJ) - vbox(ivIntraLJ,ibox)) .gt. 0.000001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Intra mol.en.: ',v(ivIntraLJ),vbox(ivIntraLJ,ibox),v(ivIntraLJ) - vbox(ivIntraLJ,ibox)
          end if
          if ( abs(v(ivStretching) - vbox(ivStretching,ibox)) .gt. 0.001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' bond vib. en.: ',v(ivStretching),vbox(ivStretching,ibox),v(ivStretching) - vbox(ivStretching,ibox)
          end if
          if ( abs(v(ivBending) - vbox(ivBending,ibox)) .gt. 0.001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Bond ben.en.: ',v(ivBending),vbox(ivBending,ibox),v(ivBending) - vbox(ivBending,ibox)
          end if
          if ( abs(v(ivTorsion) - vbox(ivTorsion,ibox)) .gt. 0.001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Torsion.en.: ',v(ivTorsion),vbox(ivTorsion,ibox),v(ivTorsion) - vbox(ivTorsion,ibox)
          end if
          if ( abs(v(ivExt) - vbox(ivExt,ibox)) .gt. 0.0001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Externa.en.: ',v(ivExt),vbox(ivExt,ibox),v(ivExt) - vbox(ivExt,ibox)
          end if
          if ( abs(v(ivElect) - vbox(ivElect,ibox)) .gt. 0.000001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Coulomb.en.: ',v(ivElect),vbox(ivElect,ibox),v(ivElect) - vbox(ivElect,ibox)
          end if
          if ( abs(v(ivFlucq) - vbox(ivFlucq,ibox)) .gt. 0.0001) then
             write(io_output,*) '### problem  ###'
             write(io_output,*) ' Fluc Q en.: ',v(ivFlucq),vbox(ivFlucq,ibox),v(ivFlucq) - vbox(ivFlucq,ibox)
          end if
          if ( abs(v(iv3body) - vbox(iv3body,ibox) ) .gt.0.001) then
             write(io_output,*) '### problem ###'
             write(io_output,*) ' 3-body en.: ',v(iv3body),vbox(iv3body,ibox),v(iv3body) - vbox(iv3body,ibox)
          end if
          if ( ldielect ) then
             if ( abs(dipolexo - dipolex(ibox)) .gt. 0.0001) then
                write(io_output,*) '### problem  ###'
                write(io_output,*) ' Dipole X: ',dipolexo,dipolex(ibox),dipolexo - dipolex(ibox)
             end if
          end if
          if (lmipsw) then
             if (abs(vwellipsw-vbox(ivWellIpswb,ibox)).gt.0.001) then
                write(io_output,*) '### problem  ###'
                write(io_output,*) ' well en.: ',vwellipsw,vbox(ivWellIpswb,ibox),vwellipsw-vbox(ivWellIpswb,ibox)
             end if
          end if
       end if
    end do

! -------------------------------------------------------------------
    if (myid.eq.rootid) then
       ! Write out final configurations
       write(io_output,*)
       write(io_output,"(A,"//format_n(nbox,"(1X,F23.10)")//")") ' vstart       =',(vstart(i) ,i=1,nbox)
       write(io_output,"(A,"//format_n(nbox,"(1X,F23.10)")//")") ' vend         =',(vend(i)   ,i=1,nbox)
       write(io_output,"(A,"//format_n(nbox,"(1X,F23.10)")//")") ' vbox         =',(vbox(ivTot,i) ,i=1,nbox)

       ! write out the final configuration for each box, Added by Neeraj 06/26/2006 3M ***
       io_config=get_iounit()
       do ibox = 1,nbox
          string = integer_to_string(ibox)
          file_config = "box"//string(1:len_trim(string))//"config"//run_num(1:len_trim(run_num))//suffix//".xyz"
          !write(file_config,'("box",I1.1,"config",I1.1,A,".xyz")') ibox,run_num,suffix
          open(unit=io_config,access='sequential',action='write',file=file_config,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open box config file '//trim(file_config),jerr)

          nummol = 0
          do i = 1,nchain
             if (nboxi(i).eq.ibox) then
                nummol = nummol + nunit(moltyp(i))
             end if
          end do
          write(io_config,*) nummol
          write(io_config,*)
          do i = 1,nchain
             if(nboxi(i).eq.ibox) then
                imolty = moltyp(i)
                do ii = 1,nunit(imolty)
                   ntii = ntype(imolty,ii)
                   write(io_config,'(a4,5x,3f15.4)') chemid(ntii),rxu(i,ii),ryu(i,ii),rzu(i,ii)
                end do
             end if
          end do
          close(io_config)
       end do

       if (N_add.gt.0) then
          if (lbranch(moltyp2add)) then
             ! If we have an input structure, use it. Setting lrigid to .true. will cause gcmc_setup to call setup_molecule_config and setup_molecule_config will use the stored structure instead of growing a new one.
             lrigid(moltyp2add)=.true.
          else if (.not.lrigid(moltyp2add)) then
             ! otherwise grow it
             call setup_molecule_config(moltyp2add,nchain+1)
             if (N_add.gt.1) then
                do ii = 1, nunit(moltyp2add)
                   rxu(nchain+2:nchain+N_add,ii) = rxu(nchain+1,ii)
                   ryu(nchain+2:nchain+N_add,ii) = ryu(nchain+1,ii)
                   rzu(nchain+2:nchain+N_add,ii) = rzu(nchain+1,ii)
                end do
             end if
          end if

          do j = 1, N_add
             call gcmc_setup(moltyp2add,box2add,i,itemp)
             rxu(i,1) = random(-1) * boxlx(box2add)
             ryu(i,1) = random(-1) * boxly(box2add)
             rzu(i,1) = random(-1) * boxlz(box2add)
             do ii = 2, nunit(moltyp2add)
                rxu(i,ii) = rxu(i,1) + rxu(i,ii)
                ryu(i,ii) = ryu(i,1) + ryu(i,ii)
                rzu(i,ii) = rzu(i,1) + rzu(i,ii)
             end do
          end do
       else if (N_add.lt.0) then
          do j = 1, -N_add
             i = parbox(j,box2add,moltyp2add)
             call gcmc_exchange(i,0)
          end do
       end if

       ! write out the final configuration from the run
       file_config = "config"//run_num(1:len_trim(run_num))//suffix//".dat"
       !write(file_config,'("config",I1.1,A,".dat")') run_num,suffix
       call dump(file_config)

! -------------------------------------------------------------------
       ! Calculate and write out running averages ***
       do ibox=1,nbox
          if ( lpbcz ) then
             do itype = 1, nmolty
                ! number density
                acdens(ibox,itype)=1000.0E0_dp*acdens(ibox,itype)
             end do

             ! sum over all types of molecules
             temacd = sum(acdens(ibox,1:nmolty))

             ! molar volume in mL/mol
             molvol(ibox) = N_Avogadro*1E-21_dp/temacd

             ! specific density
             temspd = 0.0E0_dp
             do itype = 1, nmolty
                temspd = temspd + (acdens(ibox,itype)*masst(itype)/(N_Avogadro*1E-21_dp))
             end do
             speden(ibox) = temspd
          else
             do itype = 1, nmolty
                ! number density
                acdens(ibox,itype)=100.0E0_dp*acdens(ibox,itype)
             end do

             temacd = sum(acdens(ibox,1:nmolty))

             ! molar volume
             molvol(ibox) = 100.0E0_dp / temacd
          end if

          ! energies
          acvsq(:,ibox) = acvsq(:,ibox) - acv(:,ibox) ** 2
          acvkjmol(:,ibox) = acvkjmol(:,ibox)*R_gas/1000_dp

          ! if ( ldielect ) then
          !    flucmom(ibox) = acdipolesq(4,ibox)-acdipole(4,ibox)*acdipole(4,ibox)
          !    flucmom(ibox) = 6.9994685465110493E5_dp*flucmom(ibox)*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
          !    flucmom2(ibox) = acdipolesq(4,ibox)-acdipole(1,ibox)*acdipole(1,ibox)-acdipole(2,ibox)*acdipole(2,ibox)-acdipole(3,ibox)*acdipole(3,ibox)
          !    flucmom2(ibox) = 6.9994685465110493E5_dp*flucmom2(ibox)*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
          ! end if

          ! chemical potential
          ! This expression is correct for the NVT(-Gibbs) and NPT(-Gibbs) ensembles
          ! But the use of the dual cutoff method will introduce some inaccuracies
          do itype = 1, nmolty
             if(bnchem(ibox,itype).gt.0) then
                debroglie = debroglie_factor*sqrt(beta/masst(itype))
                ! determine how many steps it takes to grow molecule
                ! not counting the first inserted bead
                igrow = nugrow(itype)
                if (.not. lrigid(itype)) then
                   call schedule(igrow,itype,steps,1,0,2,0)
                   tmp=real(nchoi1(itype)*nchoih(itype),dp)*(real(nchoi(itype),dp)**steps)*(debroglie**3)
                else
                   call schedule(igrow,itype,steps,1,0,4,0)
                   tmp=real(nchoi1(itype)*nchoir(itype)*nchoih(itype),dp)*(real(nchoi(itype),dp)**steps)*(debroglie**3)
                end if
                acchem(ibox,itype)=(-1.0E0_dp/beta)*log((acchem(ibox,itype)/bnchem(ibox,itype))/tmp)
                itel = 2 + nEnergy + itype
                baver(itel,ibox,:)=(-1.0E0_dp/beta)*log(baver(itel,ibox,:)/tmp)
                if (lucall) call write_deltaG_map(io_output,itype,tmp)
             end if
          end do
       end do

       write(io_output,*)
       write(io_output,"(A,"//format_n(nbox,"(A,I2)")//")") ' Averages and fluctuations                           '&
        ,('       Box ',i,i=1,nbox)
       write(io_output,*)
       write(io_output,"(A,"//format_n(nbox,"(1X,F12.2)")//")") ' pressure                                      [kPa] ='&
        ,(acpres(i),i=1,nbox)
       write(io_output,"(A,"//format_n(nbox,"(1X,F12.6)")//")") ' pressure                         [simulation units] ='&
        ,((acpres(i)*MPa2SimUnits*1E-3_dp),i=1,nbox)
       write(io_output,"(A,"//format_n(nbox,"(1X,F12.4)")//")") ' surface tension                              [mN/m] ='&
        ,(acsurf(i),i=1,nbox)
       write(io_output,"(A,"//format_n(nbox,"(1X,F12.5)")//")") ' compress factor                                     ='&
        ,(accomp(i),i=1,nbox)

       do itype = 1,nmolty
          write(io_output,"(A,I2,A,A10,A,"//format_n(nbox,"(1X,F12.3)")//")") ' chem. potential of type  ',itype,' '&
          ,molecname(itype),'          [K] =',(acchem(i,itype),i=1,nbox)
       end do

       do i=1,3
          write(io_output,"(A,"//format_n(nbox,"(1X,F12.3)")//")") ' boxlength                                       [A] ='&
           ,(acboxl(ibox,i),ibox=1,nbox)
       end do

       if (ANY(lsolid(1:nbox).and..not.lrect(1:nbox))) then
          do i=1,3
             write(io_output,"(A,"//format_n(nbox,"(1X,F12.3)")//")") ' box angle                                     [deg] ='&
              ,(acboxa(ibox,i)*raddeg,ibox=1,nbox)
          end do
       end if

       do itype = 1, nmolty
          write(io_output,"(A,I2,A,A10,A,"//format_n(nbox,"(1X,F12.3)")//")") ' no. of chains of type    ',itype,' '&
           ,molecname(itype),'              =',(acnbox(i,itype),i=1,nbox)
       end do
       if (lpbcz) then
          write(io_output,"(A,"//format_n(nbox,"(1X,F12.3)")//")") ' molar volume                             [cm^3/mol] ='&
           ,(molvol(i),i=1,nbox)
          write(io_output,"(A,"//format_n(nbox,"(1X,F12.6)")//")") ' specific density                           [g/cm^3] ='&
           ,(speden(i),i=1,nbox)
          do itype = 1, nmolty
             write(io_output,"(A,I2,A,A10,A,"//format_n(nbox,"(1X,F12.5)")//")") ' number density of type   ',itype,' '&
              ,molecname(itype),' [chain/nm^3] =',(acdens(i,itype),i=1,nbox)
             if (lexpand(itype)) then
                do itype2=1,numcoeff(itype)
                   write(io_output,"(A,I2,A,I4,A,2F12.5)") ' number density of type   ',itype,' eetype ',itype2,'  ='&
                    ,acdens(itype,itype)*acnbox2(itype,itype,itype2)/(acnbox(itype,itype)*acmove)&
                    ,acnbox2(itype,itype,itype2)/(acnbox(itype,itype)*acmove)
                end do
             end if
          end do
       else
          write(io_output,"(A,"//format_n(nbox,"(1X,F12.4)")//")") ' area per chain                     [A^2/chain] ='&
           ,(molvol(i),i=1,nbox)
          do itype = 1, nmolty
             write(io_output,"(A,I2,A,A10,A,"//format_n(nbox,"(1X,F12.6)")//")") ' number density of type  ',itype,' '&
              ,molecname(itype),' [chain/nm^2] =',(acdens(i,itype),i=1,nbox)
          end do
       end if
       do itype = 1, nmolty
          write(io_output,"(A,I2,A,A10,A,"//format_n(nbox,"(1X,F12.7)")//")") ' molfraction of type      ',itype,' '&
           ,molecname(itype),'              =',(molfra(i,itype),i=1,nbox)
       end do
       do itype = 1, nmolty
          write(io_output,"(A,I2,A,A10,A,"//format_n(nbox,"(1X,F12.3)")//")") ' mean sete length of type ',itype,' '&
           ,molecname(itype),'        [A^2] =',(asetel(i,itype),i=1,nbox)
       end do

       if (lucall) call write_prop_widom(io_output)

       write(io_output,*)
       do j=1,11
          ! only 1 to 11 is the energy information
          write(io_output,"(A14,A,"//format_n(nbox,"(1X,F13.2)")//","//format_n(nbox,"(1X,F10.2)")//")") vname(j)&
           ,'[K per system and kJ/mol per chain] =',acv(j,1:nbox),acvkjmol(j,1:nbox)
       end do

       write(io_output,*)
       write(io_output,"(A,"//format_n(nbox,"(1X,F11.2)")//")") ' fluctuation in <vtot> =',(sqrt(acvsq(1,i)),i=1,nbox)

       write(io_output,*)
       ! Output 2nd virial coefficient data
2000   if (lvirial) then
           nImages = real(nstep/imv,dp)
           write(io_output,*) 'At temperature of ', virtemp
           bvirial = binvir/nImages
           write(io_output,*) 'bvirial ', bvirial, ' [ A^3 / molecule ]'
           write(io_output,*) 'bvirial ', bvirial*N_Avogadro*1E-24_dp, ' [ cm^3 / mole ]'
       end if

       ! solute values
       write(io_output,'(A)') ' type  box       vinter       vintra         vtor        vbend        vtail'
       do itype = 1, nmolty
          do ibox = 1, nbox
             if (solcount(ibox,itype).gt.0) then
                write(io_output,"(2I5,5(1X,F12.5))") itype,ibox,avsolinter(ibox,itype),avsolintra(ibox,itype)&
                 ,avsoltor(ibox,itype),avsolbend(ibox,itype),avsolelc(ibox,itype)
             else
                write(io_output,"(2I5,5(1X,F12.5))") itype,ibox,0.0,0.0,0.0,0.0,0.0
             end if
          end do
       end do

! -------------------------------------------------------------------
       ! Calculate statistical uncertainties from block averages
       if ( nblock .ge. 2 ) then
          ! global averages -
          do i = 1,nprop
             do ibox = 1,nbox
                call calculate_statistics(baver(i,ibox,1:nblock),aver(i,ibox),stdev(i,ibox),errme(i,ibox))
             end do
          end do
          do i = 1,nprop1
             do ibox = 1,nbox-1
                do jbox = ibox+1,nbox
                   call calculate_statistics(baver1(i,ibox,jbox,1:nblock),aver1(i,ibox,jbox),stdev1(i,ibox,jbox)&
                    ,errme1(i,ibox,jbox))
                end do
             end do
          end do

          ! write out the heat of vaporization and solubility parameters
          write(io_output,*)
          do ibox = 1,nbox-1
             do jbox = ibox+1,nbox
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' H_vap      [kJ/mol] btwn box   ',ibox,' and ',jbox,' ='&
                 ,acsolpar(1,ibox,jbox),stdev1(1,ibox,jbox),errme1(1,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' H_vap LJ  [kJ/mol] btwn box    ',ibox,' and ',jbox,' ='&
                 ,acsolpar(2,ibox,jbox),stdev1(2,ibox,jbox),errme1(2,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' H_vap Coul [kJ/mol] btwn box   ',ibox,' and ',jbox,' ='&
                 ,acsolpar(3,ibox,jbox),stdev1(3,ibox,jbox),errme1(3,ibox,jbox)
                ! write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' DeltaU Ext [kJ/mol] btwn box   ',ibox,' and ',jbox,' =',acsolpar(10,ibox,jbox),stdev1(10,ibox,jbox),errme1(10,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' pdV        [kJ/mol] btwn box   ',ibox,' and ',jbox,' ='&
                 ,acsolpar(11,ibox,jbox),stdev1(11,ibox,jbox),errme1(11,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' CED [cal/cc]   btwn box        ',ibox,' and ',jbox,' ='&
                 ,acsolpar(4,ibox,jbox),stdev1(4,ibox,jbox),errme1(4,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' CED_LJ[cal/cc] btwn box        ',ibox,' and ',jbox,' ='&
                 ,acsolpar(5,ibox,jbox),stdev1(5,ibox,jbox),errme1(5,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' CED_Coul[cal/cc] btwn box      ',ibox,' and ',jbox,' ='&
                 ,acsolpar(6,ibox,jbox),stdev1(6,ibox,jbox),errme1(6,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' HSP [(cal/cc)^1/2]  btwn box   ',ibox,' and ',jbox,' ='&
                 ,acsolpar(7,ibox,jbox),stdev1(7,ibox,jbox),errme1(7,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' HSP_LJ[(cal/cc)^1/2] btwn box  ',ibox,' and ',jbox,' ='&
                 ,acsolpar(8,ibox,jbox),stdev1(8,ibox,jbox),errme1(8,ibox,jbox)
                write(io_output,"(2(A,I2),A,3(1X,F14.4))") ' HSP_Cou[(cal/cc)^1/2] btwn box ',ibox,' and ',jbox,' ='&
                 ,acsolpar(9,ibox,jbox),stdev1(9,ibox,jbox),errme1(9,ibox,jbox)
             end do
          end do

          write(io_output,*)
          ! specific density
          do ibox = 1, nbox
             write(io_output,"(A,I2,A,3(1X,E12.5))") ' specific density box ',ibox,' =',aver(1,ibox),stdev(1,ibox),errme(1,ibox)
          end do

          ! system volume
          itel = nEnergy + 4*nmolty + 4
          do ibox = 1, nbox
             write(io_output,"(A,I2,A,3(1X,E12.5))") ' system volume    box ',ibox,' =',aver(itel,ibox),stdev(itel,ibox)&
              ,errme(itel,ibox)
          end do

          ! pressure
          do ibox = 1, nbox
             write(io_output,"(A,I2,A,3(1X,G12.5))") ' pressure         box ',ibox,' =',acpres(ibox),stdev(2,ibox),errme(2,ibox)
          end do

          ! surface tension
          itel = nEnergy + 4*nmolty + 3
          do ibox = 1, nbox
             write(io_output,"(A,I2,A,3(1X,F12.5))") ' surface tension  box ',ibox,' =',acsurf(ibox),stdev(itel,ibox)&
              ,errme(itel,ibox)
          end do

         ! compressibility factor
          itel = nEnergy + 4*nmolty + 7
          do ibox = 1, nbox
             write(io_output,"(A,I2,A,3(1X,F12.5))") ' compressibility  box ',ibox,' =',accomp(ibox),stdev(itel,ibox)&
              ,errme(itel,ibox)
          end do

          write(io_output,*)
          ! energies
          do ibox = 1, nbox
             do j=3,13
                ! only 1 to 10 is the energy information
                write(io_output,"(A17,A,I2,A,3(1X,E12.5))") vname(j-2),' box ',ibox,' =',aver(j,ibox),stdev(j,ibox),errme(j,ibox)
             end do
          end do

          write(io_output,*)
          ! Enthalpy
          do ibox = 1,nbox
             j = nEnergy + 4*nmolty + 5
             write(io_output, "(A,I2,A,3(1X,F12.4))") ' Enthalpy Inst.[kJ/mol] for box ',ibox,' =',acEnthalpy(ibox)&
              ,stdev(j,ibox),errme(j,ibox)
             j = nEnergy + 4*nmolty + 6
             write(io_output,"(A,I2,A,3(1X,F12.4))") ' Enthalpy Ext. [kJ/mol] for box ',ibox,' =',acEnthalpy1(ibox)&
              ,stdev(j,ibox),errme(j,ibox)
          end do

          write(io_output,*)
          ! residual heat capacity, in (J2/mol)
          if (.not.lgibbs) then
             tmp=(enthalpy2-enthalpy*enthalpy)/real(nchain,dp)*R_gas/(temp**2)
             if(lnpt) then
                write(io_output,'(A,G16.9)') ' Cp residual(J/Kmol) = ',tmp
                write(io_output,'(A,G16.9)') '   H2 = ',enthalpy2
                write(io_output,'(A,G16.9)') '   H  = ',enthalpy
             else
                write(io_output,'(A,G16.9)') ' Cv residual(J/Kmol) = ',tmp
                write(io_output,'(A,G16.9)') '   E2 = ',enthalpy2
                write(io_output,'(A,G16.9)') '   E  = ',enthalpy
             end if
          end if

          write(io_output,*)
          ! chemical potential
          do itype = 1, nmolty
             itel = 2 + nEnergy + itype
             do ibox = 1, nbox
                write(io_output,"(2(A,I2),A,3(1X,F12.3))") ' chemical potential  itype ',itype,' box ',ibox,' = '&
                 ,aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
             end do
          end do

          ! square end-to-end length
          do itype = 1, nmolty
             itel = 2 + nEnergy + nmolty + itype
             do ibox = 1, nbox
                write(io_output,"(2(A,I2),A,3(1X,F12.3))") ' mean sete length    itype ',itype,' box ',ibox,' = '&
                 ,aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
             end do
          end do

          ! number density
          do itype = 1, nmolty
             itel = 2 + nEnergy + 2*nmolty + itype
             do ibox = 1, nbox
                if (lpbcz) then
                   write(io_output,"(2(A,I2),A,3(1X,E12.5))") ' number density      itype ',itype,' box ',ibox,' = '&
                    ,1.0E3_dp*aver(itel,ibox),1.0E3_dp*stdev(itel,ibox),1.0E3_dp*errme(itel,ibox)
                else
                   write(io_output,"(2(A,I2),A,3(1X,E12.5))") ' number density      itype ',itype,' box ',ibox,' = '&
                    ,1.0E2_dp*aver(itel,ibox),1.0E2_dp*stdev(itel,ibox),1.0E2_dp*errme(itel,ibox)
                end if
                if (lexpand(itype).and.acnbox(ibox,itype).gt.0.5) then
                   do itype2=1,numcoeff(itype)
                      molfrac = acnbox2(ibox,itype,itype2) /(acmove*acnbox(ibox,itype))
                      write(io_output,"(2(A,I2),A,2(1X,E12.5))") ' number density      itype ',itype,' typ ',itype2,' = '&
                       ,1.0E3_dp*aver(itel,ibox)*molfrac,molfrac
                   end do
                end if
             end do
          end do

          ! molfraction
          do itype = 1, nmolty
             itel = 2 + nEnergy + 3*nmolty + itype
             do ibox = 1, nbox
                write(io_output,"(2(A,I2),A,3(1X,F12.7))") ' mole fraction       itype ',itype,' box ',ibox,' = '&
                 ,aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
             end do
          end do

          if (lgibbs) then
             ! write density results in fitting format ---
             do ibox = 1, nbox-1
                do jbox = ibox+1,nbox
                   if (speden(ibox).lt.speden(jbox)) then
                      ig = ibox
                      il = jbox
                   else
                      ig = jbox
                      il = ibox
                   end if
                   ! write(41,"(2x,f6.1,4(1x,e12.5))") temp,aver(1,ig),stdev(1,ig),aver(1,il),stdev(1,il)
                   ! write ostwald values for each moltyp
                   gconst = R_gas/(1E3_dp*beta)
                   do itype = 1,nmolty
                      itel = 2 + nEnergy + 2 * nmolty + itype
                      ostwald = aver(itel,il)/aver(itel,ig)
                      stdost  = ostwald * sqrt( (stdev(itel,il)/aver(itel,il))**2 + (stdev(itel,ig)/aver(itel,ig))**2 )
                      ! write(42,*) nunit(itype),ostwald,stdost
                      ! write(43,*) nunit(itype),-(gconst*log(ostwald)) + (eta2(ig,itype) - eta2(il,itype)) / 120.27167,gconst*stdost/ostwald
                      write(io_output,"(3(A,I2),A,2(1X,F15.6))") ' Ostwald Coefficient itype ',itype,' between box ',ig&
                       ,' and ',il,' = ',ostwald,stdost
                      write(io_output,"(3(A,I2),A,2(1X,F15.6))") ' Free Enrgy of Trans itype ',itype,' between box ',ig&
                       ,' and ',il,' [kJ/mol] =',-(gconst*log(ostwald))+(eta2(ig,itype)-eta2(il,itype))*R_gas*1E-3_dp&
                       ,gconst*stdost/ostwald
                   end do
                end do
             end do
          end if

          write(io_output,*)

          if (lucall) call write_prop_widom_with_stats(io_output,nblock)

          write(io_output,*)
          write(io_output,*) '-----block averages ------'
          do ibox=1,nbox
             write(io_output,"('  ------------ box: ' ,I2,/, ' block    energy     density    pressure       Z        surf ten   mol fracs')")&
              ibox
             do nbl = 1, nblock
                ! changed so output the same for all ensembles
                ! 06/08/09 KM
                write(io_output,"(1X,I3,"//format_n(nmolty+5,"(1X,E11.4)")//")") nbl,baver(3,ibox,nbl),baver(1,ibox,nbl)&
                 ,baver(2,ibox,nbl),baver(7+nEnergy+4*nmolty,ibox,nbl),baver(nEnergy+4*nmolty+3,ibox,nbl),(baver(2+nEnergy+3*nmolty+zzz,ibox,nbl),zzz=1,nmolty)
             end do
             if (lmipsw) then
                write(io_output,*) 'lambdais', lambdais
                write(io_output,*) 'maginn interphase switch integrand'
                do nbl = 1, nblock
                   write(io_output,*) nbl,baver(nprop,ibox,nbl)
                end do
             end if
          end do
       end if

       ! KM 01/10 remove analysis
       ! if (ianalyze.le.nstep) then
       !    call analysis(2)
       ! end if

       ! ee prob
       IF(lexpee) then
          write(io_output,*)
          write(io_output,*) 'probability of each mstate in ee'
          do nnn = 1, fmstate
             write(io_output,"(i5,1x,i10)") nnn,ee_prob(nnn)
          end do
       end if

       write(io_output,*) 'Program ended at ',time_date_str()
       close(io_output)
    end if

    return
  end subroutine monola

#ifdef w_a
#undef w_a
#endif

#ifdef w_l
#undef w_l
#endif

#ifdef w_i
#undef w_i
#endif

#ifdef w_r
#undef w_r
#endif

#ifdef __TRADITIONAL_CPP__
#define w_a(x) write(io_output,'(1X,A,": ",A)') "x",x
#define w_l(x) write(io_output,'(1X,A,": ",L2)') "x",x
#define w_i(x) write(io_output,'(1X,A,": ",I0)') "x",x
#define w_r(x) write(io_output,'(1X,A,": ",G16.9)') "x",x
#else
#define w_a(x) write(io_output,'(1X,A,": ",A)') #x,x
#define w_l(x) write(io_output,'(1X,A,": ",L2)') #x,x
#define w_i(x) write(io_output,'(1X,A,": ",I0)') #x,x
#define w_r(x) write(io_output,'(1X,A,": ",G16.9)') #x,x
#endif
!> \brief Read input data and initializes the positions
!>
!> reads a starting configuration from unit 7
!> calculates interaction table
!> \remarks It's usually best to print out information right after it is read in and before doing any other operations (which may fail), so that in the case of program error users will know the progress of input processing
  subroutine readdat(file_in)
    use var_type,only:default_string_length
    use const_math,only:raddeg, onepi
    use const_phys,only:MPa2SimUnits,debroglie_factor
    use util_random,only:ranset
    use util_string,only:integer_to_string,real_to_string,uppercase,format_n
    use util_runtime,only:err_exit
    use util_timings,only:time_date_str,time_now
    use util_files,only:get_iounit,readLine
    use util_search,only:indexOf
    use util_memory,only:reallocate
    use util_mp,only:mp_bcast
    use util_kdtree,only:read_kdtree,allocate_kdtree
    use sim_particle,only:allocate_neighbor_list
    use sim_initia,only:get_molecule_config,setup_system_config,setup_cbmc_bend
    ! Q. Paul C. -- add setup_cbmc_bend for tabulated CBMC bending growth in the above line
    use zeolite
    use energy_kspace,only:calp,allocate_kspace,k_max_l,k_max_m,k_max_n,compute_kmax
    use energy_pairwise,only:read_ff,init_ff,type_2body,vdW_nParameter,nonbond_type
    use energy_intramolecular,only:bonds,angles,dihedrals,allocate_energy_bonded
    use energy_3body,only:readThreeBody
    use energy_4body,only:readFourBody
    use energy_sami
    use moves_simple,only:init_moves_simple,averageMaximumDisplacement
    use moves_volume,only:init_moves_volume,allow_cutoff_failure,restore_displ_transl
    use moves_cbmc,only:init_cbmc,allocate_cbmc,llplace
    use moves_ee,only:init_ee,numcoeff,sigm,epsil
    use transfer_shared,only:read_transfer
    use transfer_swap,only:init_swap
    use transfer_swatch,only:init_swatch
    use prop_widom,only:read_prop_widom

    character(LEN=*),intent(in)::file_in

    real,allocatable::ofscale(:),ofscale2(:),qbox(:)
    integer,allocatable::ncarbon(:),inclmol(:),inclbead(:,:),inclsign(:),ainclmol(:),ainclbead(:,:),a15t(:),idummy(:),temphe(:)&
     ,nures(:)
    logical,allocatable::lhere(:)

    character(LEN=default_path_length)::file_input,file_restart,file_struct,file_run,file_movie,file_solute,file_traj&
     ,file_box_movie,file_box_movie_pdb,file_cell,file_cbmc_bend
    ! Q. Paul C. --adding file_cbmc_bend for tabulated CBMC bending growth

    character(LEN=default_string_length)::line_in
    integer::io_input,io_restart,jerr,seed,ij,ii,jj,i,j,k,ncres,nmtres,iensem,inpbc,im,ibox,izz,z,imix
    logical::lprint,L_Ewald_Auto,lmixlb,lmixjo,lmixwh,lmixkong,lsetup,linit,lreadq,lfound,ltmp,tmp_logical
    real::fqtemp,dum,pm,pcumu,debroglie,qtot,min_boxl,rmflucq,rtmp
    !* PARAMETER FOR IONIC SYSTEMS
    logical::lionic !< if LIONIC=.TRUE. System contains charged species, so system may not neutral
    ! variables added (3/24/05) for scaling of 1-4 interactions
    integer::nexclu,inclnum,ainclnum
    ! variables for histograms
    integer::temnc,imol,iutemp,imolty,itype
    ! Variables added (6/30/2006) for fort.4 consistency check
    integer::numvib,numbend,numtor,vib1,bend2,bend3,tor2,tor3,tor4
    integer::vibtype,bendtype,tortype
    ! real::temx(nmax,numax),temy(nmax,numax),temz(nmax,numax)
    ! KM variable added when analysis removed
    integer::nhere
    ! temporary FQ variables
    integer::fqmolty,fqtyp,fqcrosstyp,ntii,ntjj,ntij
    logical::needMFsection

    namelist /io/ file_input,file_restart,file_struct,file_run,file_movie,file_solute,file_traj,io_output&
     ,run_num,suffix,L_movie_xyz,L_movie_pdb,file_cbmc_bend,checkpoint_interval,checkpoint_copies,use_checkpoint,ltraj
    namelist /system/ lnpt,lgibbs,lgrand,lanes,lvirial,lmipsw,lexpee,ldielect,lpbc,lpbcx,lpbcy,lpbcz,lfold,lijall,lchgall,lewald&
     ,lrecip,lcutcm,ltailc,lshift,ldual,L_Coul_CBMC,lneigh&
     ,lexzeo,lslit,lgraphite,lsami,lmuir,lelect_field,lgaro,lionic,L_Ewald_Auto,lmixlb,lmixjo&
     ,lmixwh,lmixkong,L_spline,L_linear,L_vib_table,L_bend_table,L_elect_table,L_cbmc_bend, lkdtree,losmoticnvt
    namelist /mc_shared/ seed,nbox,nmolty,nchain,nmax,nstep,time_limit,lstop,iratio,rmin,softcut&
     ,linit,lreadq,N_add,box2add,moltyp2add
    namelist /analysis/ iprint,imv,iblock,iratp,idiele,iheatcapacity,ianalyze&
     ,nbin,lrdf,lintra,lstretch,lgvst,lbend,lete,lrhoz,bin_width&
     ,lucall,nvirial,starvir,stepvir,ntemp,virtemp
    namelist /mc_flucq/ taflcq,fqtemp,rmflucq,pmflcq,pmfqmt,lflucq,lqtrans,fqegp,nchoiq,nswapq
    namelist /gcmc/ B,nequil,ninstf,ninsth,ndumph
! ===================================================================
    !> read project-wide parameters
    if (myid.eq.rootid) then
       io_input=get_iounit()
       open(unit=io_input,access='sequential',action='read',file=file_in,form='formatted',iostat=jerr,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open input file '//trim(file_in),jerr)
    end if
! -------------------------------------------------------------------
    !> read namelist io
    file_input='fort.4'
    file_restart='fort.77'
    file_struct='input_struc.xyz'
    file_run='run1a.dat'
    file_movie='movie1a.dat'
    file_solute='fort.11'
    file_traj='fort.12'
    file_cbmc_bend='cbmc_bend_table.dat' ! Q. Paul C. -- tabulated CBMC bending growth

    if (myid.eq.rootid) then
       read(UNIT=io_input,NML=io,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: io',jerr)
    end if

    call mp_bcast(file_input,rootid,groupid)
    call mp_bcast(file_restart,rootid,groupid)
    call mp_bcast(file_struct,rootid,groupid)
    call mp_bcast(file_run,rootid,groupid)
    call mp_bcast(file_movie,rootid,groupid)
    call mp_bcast(file_solute,rootid,groupid)
    call mp_bcast(file_traj,rootid,groupid)
    call mp_bcast(io_output,1,rootid,groupid)
    call mp_bcast(run_num,rootid,groupid)
    call mp_bcast(suffix,rootid,groupid)
    call mp_bcast(L_movie_xyz,1,rootid,groupid)
    call mp_bcast(L_movie_pdb,1,rootid,groupid)
    call mp_bcast(ltraj,1,rootid,groupid)
    call mp_bcast(file_cbmc_bend,rootid,groupid) !Q.Paul C. -- tabulated CBMC bending growth
    call mp_bcast(checkpoint_interval,1,rootid,groupid)
    call mp_bcast(checkpoint_copies,1,rootid,groupid)
    call mp_bcast(use_checkpoint,1,rootid,groupid)

    if (myid.eq.rootid.and..not.use_checkpoint) then
       lprint=.true.
    else
       lprint=.false.
    end if

    ! Output unit: if 6 or 0, write to stdout/stderr; otherwise, write to a user designated file
    if (io_output.eq.5) then
       call err_exit(__FILE__,__LINE__,'unit 5 is for standard input',myid+1)
    else if(io_output.ne.6.and.io_output.ne.0.and.myid.eq.rootid) then
       io_output=get_iounit()
       open(unit=io_output,access='stream',action='write',file=file_run,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) then
          call err_exit(__FILE__,__LINE__,'cannot open output file '//trim(file_run),jerr)
       end if
    end if

    if (lprint) then
       write(io_output,'(A,A)') 'Program started at ',time_date_str()
       write(io_output,'(A,I0)') 'Number of processors: ', numprocs
       write(io_output,'(A,I0)') 'Threads per processor: ',thread_num
       call output_version(io_output)
       w_a(run_num)
       w_a(suffix)
       w_l(L_movie_xyz)
       w_l(L_movie_pdb)
    end if
! -------------------------------------------------------------------
    !> read namelist system
    lionic=.false.
    L_Ewald_Auto=.true.
    lmixlb=.true.
    lmixjo=.false.
    lmixwh=.false.
    lmixkong=.false.
    L_cbmc_bend=.false. ! Q. Paul C. -- for tabulated CBMC bending growth
    lkdtree=.false.

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=system,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: system',jerr)
    end if

    call mp_bcast(lnpt,1,rootid,groupid)
    call mp_bcast(lgibbs,1,rootid,groupid)
    call mp_bcast(lgrand,1,rootid,groupid)
    call mp_bcast(lanes,1,rootid,groupid)
    call mp_bcast(lvirial,1,rootid,groupid)
    call mp_bcast(lmipsw,1,rootid,groupid)
    call mp_bcast(lexpee,1,rootid,groupid)
    call mp_bcast(ldielect,1,rootid,groupid)
    call mp_bcast(lpbc,1,rootid,groupid)
    call mp_bcast(lpbcx,1,rootid,groupid)
    call mp_bcast(lpbcy,1,rootid,groupid)
    call mp_bcast(lpbcz,1,rootid,groupid)
    call mp_bcast(lfold,1,rootid,groupid)
    call mp_bcast(lijall,1,rootid,groupid)
    call mp_bcast(lchgall,1,rootid,groupid)
    call mp_bcast(lewald,1,rootid,groupid)
    call mp_bcast(lrecip,1,rootid,groupid)
    call mp_bcast(lcutcm,1,rootid,groupid)
    call mp_bcast(ltailc,1,rootid,groupid)
    call mp_bcast(lshift,1,rootid,groupid)
    call mp_bcast(ldual,1,rootid,groupid)
    call mp_bcast(L_Coul_CBMC,1,rootid,groupid)
    call mp_bcast(lneigh,1,rootid,groupid)
    call mp_bcast(lexzeo,1,rootid,groupid)
    call mp_bcast(lslit,1,rootid,groupid)
    call mp_bcast(lgraphite,1,rootid,groupid)
    call mp_bcast(lsami,1,rootid,groupid)
    call mp_bcast(lmuir,1,rootid,groupid)
    call mp_bcast(lelect_field,1,rootid,groupid)
    call mp_bcast(lgaro,1,rootid,groupid)
    call mp_bcast(lionic,1,rootid,groupid)
    call mp_bcast(L_Ewald_Auto,1,rootid,groupid)
    call mp_bcast(lmixlb,1,rootid,groupid)
    call mp_bcast(lmixjo,1,rootid,groupid)
    call mp_bcast(lmixwh,1,rootid,groupid)
    call mp_bcast(lmixkong,1,rootid,groupid)
    call mp_bcast(L_spline,1,rootid,groupid)
    call mp_bcast(L_linear,1,rootid,groupid)
    call mp_bcast(L_vib_table,1,rootid,groupid)
    call mp_bcast(L_bend_table,1,rootid,groupid)
    call mp_bcast(L_elect_table,1,rootid,groupid)
    call mp_bcast(L_cbmc_bend,1,rootid,groupid) ! Q. Paul C. -- for tabulated CBMC bending growth
    call mp_bcast(lkdtree,1,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A)') '***** PROGRAM  =  THE MAGIC BLACK BOX *****'
       iensem = 0
       if (lgrand) then
          iensem = iensem + 1
          write(io_output,'(A)') 'Grand-canonical ensemble'
       end if
       if (lgibbs) then
          iensem = iensem + 1
          if (lnpt) then
             write(io_output,'(A)') 'NPT Gibbs ensemble'
          else
             if (losmoticnvt) then
                write(io_output,'(A)') 'Osmotic NVT Gibbs ensemble'
             else
                write(io_output,'(A)') 'NVT Gibbs ensemble'
             end if
          end if
       end if
       if (iensem.eq.0) then
          if (lnpt) then
             iensem = iensem + 1
             write(io_output,'(A)') 'Isobaric-isothermal ensemble'
          else
             iensem = iensem + 1
             write(io_output,'(A)') 'Canonical ensemble'
          end if
       end if
       if (iensem.gt.1) call err_exit(__FILE__,__LINE__,'INCONSISTENT ENSEMBLE SPECIFICATION',myid+1)

       if (lanes) write(io_output,'(A)') 'ANES-MC will be performed for polarizable model'
       if (lmipsw) write(io_output,'(A)') 'Thermodynamic integration will be performed'
       if (lexpee) write(io_output,'(A)') 'Expanded ensemble moves enabled'
       if (ldielect) write(io_output,'(A)') 'Computing dielectric constant'

       if ( lpbc ) then
          write(io_output,FMT='(A)',ADVANCE='NO') 'Using periodic boundaries in'
          inpbc = 0
          if (lpbcx) then
             inpbc = inpbc + 1
             write(io_output,'(A)',ADVANCE='NO') ' x'
          end if
          if (lpbcy) then
             inpbc = inpbc + 1
             write(io_output,'(A)',ADVANCE='NO') ' y'
          end if
          if (lpbcz) then
             inpbc = inpbc + 1
             write(io_output,'(A)',ADVANCE='NO') ' z'
          end if
          if (inpbc.eq.0) call err_exit(__FILE__,__LINE__,'INCONSISTENT PBC SPECIFICATION',myid+1)
          write(io_output,'(/,I0,A)') inpbc,'-dimensional periodic box'
       else
          write(io_output,'(A)') 'Cluster mode (no PBC)'
          if (lgibbs.or.lgrand) call err_exit(__FILE__,__LINE__,'INCONSISTENT SPECIFICATION OF LPBC AND ENSEMBLE',myid+1)
       end if

       if (lfold) then
          write(io_output,'(A)') 'Particle coordinates are folded into central box'
          if (.not.lpbc) call err_exit(__FILE__,__LINE__,'INCONSISTENT SPECIFICATION OF LPBC AND LFOLD',myid+1)
       end if

       if (lijall) write(io_output,'(A)') 'All pair interactions are considered (no potential truncation)'
       if (lchgall) write(io_output,'(A)')&
        'All the inter- and intra-molecular Coulombic interactions are considered (no potential truncation)'
       if (lcutcm) then
          write(io_output,'(A)') 'Additional center-of-mass cutoff on computed rcmu'
          if (lijall) call err_exit(__FILE__,__LINE__,'cannot have lijall with lcutcm',myid+1)
          ! if (lchgall) call err_exit(__FILE__,__LINE__,'cannot have lchgall with lcutcm',myid+1)
       end if

       write(io_output,'(A)') 'CBMC simultaneously grows all beads conected to the same bead'
       write(io_output,'(A)') '   with bond lengths/angles generated from Gaussian distribution'
       write(io_output,'(A)') 'Program will call explct() for explicit-hydrogen models'
       if (ldual) write(io_output,'(A)') 'Dual Cutoff Configurational-bias Monte Carlo'
       if (L_Coul_CBMC) then
          write(io_output,'(A)') 'Coulombic interactions will be included in the Rosenbluth weights for CBMC growth'
       else
          write(io_output,'(A)')&
           'Coulombic interactions will NOT be considered during CBMC growth (added only in the final acceptance test)'
       end if

       write(io_output,'(A)') 'Coulombic inter- and intra-molecular interactions will be calculated'
       if (L_elect_table) then
          write(io_output,'(A)') '   using tabulated potential with linear interpolation'
          if (lewald) call err_exit(__FILE__,__LINE__,'L_elect_table and lewald cannot both be true',myid+1)
       else if (lewald) then
          write(io_output,'(A)') '   using Ewald-sum techniques'
       else
          write(io_output,'(A)') '   using (neutral-)group-based cutoff'
       end if

       if (ltailc) write(io_output,'(A)') '   with additional tail corrections'
       if (lshift) then
          write(io_output,'(A)') '   using a shifted potential'
          if (ltailc) then
             call err_exit(__FILE__,__LINE__,'lshift.and.ltailc!',myid+1)
          end if
       end if

       if (lgaro) then
          write(io_output,'(A)') 'Feuston-Garofalini potential -> parameters hard coded in init_garofalini'
       end if
       if (lexzeo) write(io_output,'(A)') 'Tabulated potential enabled for rigid framework zeolites/MOFs'
       if (lslit) write(io_output,'(A)') 'Featureless, planar slit surface(s) created'
       if (lgraphite) write(io_output,'(A)') '"Realistic" graphite surface'
       if (lsami) then
          write(io_output,'(A)') 'External potential for Langmuir films (SAMI)'
          write(io_output,'(A)') 'WARNING: LJ potential hard coded in susami'
          write(io_output,'(A)') 'WARNING: sets potential cut-off to 2.5sigma'
          write(io_output,'(A)') 'WARNING: has build-in tail corrections'
          if (ltailc) call err_exit(__FILE__,__LINE__,'INCONSISTENT SPECIFICATION OF LTAILC AND LSAMI',myid+1)
       end if
       if (lmuir) write(io_output,'(A)') 'External potential for Langmuir monolayers used -> parameters hard coded in sumuir'
       if (lelect_field) write(io_output,'(A)') 'External electric field applied'

       if (lmixlb) then
          write(io_output,'(A)') 'Lorentz-Berthelot combining rules apply'
       else if (lmixjo) then
          write(io_output,'(A)') 'Jorgensen combining rules apply'
       else if (lmixwh) then
          write(io_output,'(A)') 'Waldman-Hagler combining rules apply'
       else if (lmixkong) then
          write(io_output,'(A)') 'Kong combining rules apply'
       else
          call err_exit(__FILE__,__LINE__,'No combining rule is requested',myid+1)
       end if

       w_l(L_spline)
       w_l(L_linear)
       if (L_vib_table) write(io_output,'(A)') 'Stretching energy calculated using tabulated potential with linear interpolation'
       if (L_bend_table) write(io_output,'(A)') 'Bending energy calculated using tabulated potential with linear interpolation'

       write(io_output,'(A)') '*******************************************'
    end if


    imix = 0
    if (lmixlb) imix = imix + 1
    if (lmixjo) imix = imix + 1
    if (lmixwh) imix = imix + 1
    if (lmixkong) imix = imix + 1
    if (imix .gt. 1) call err_exit(__FILE__,__LINE__,'cannot use more than one combining rules!'&
     ,myid+1)

    ! KM for MPI
    if (numprocs.ne.1) then
       if (lneigh) call err_exit(__FILE__,__LINE__,'Cannot run on more than 1 processor with neighbor list!!',myid+1)
       if (lgaro) call err_exit(__FILE__,__LINE__,'Cannot run on more than 1 processor with lgaro = .true.!!',myid+1)
    end if
! -------------------------------------------------------------------
    !> read force field parameters, including intermolecular potentials and intramolecular bonded interactions (stretching, bending, torsion)
    call read_ff(io_input,lmixlb,lmixjo,lmixwh,lmixkong)
! -------------------------------------------------------------------
    if (myid.eq.rootid) close(io_input)
! ===================================================================
    !> read independent-run-specific parameters
    if (myid.eq.rootid) then
       io_input=get_iounit()
       open(unit=io_input,access='sequential',action='read',file=file_input,form='formatted',iostat=jerr,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open main input file '//trim(file_input),myid+1)
    end if
! -------------------------------------------------------------------
    !> read namelist mc_shared
    seed=time_now()
    linit=.false.
    lreadq=.false.
    nmax=0

    if (myid.eq.rootid) then
       read(UNIT=io_input,NML=mc_shared,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_shared',jerr)
    end if

    call mp_bcast(seed,1,rootid,groupid)
    call mp_bcast(nbox,1,rootid,groupid)
    call mp_bcast(nmolty,1,rootid,groupid)
    call mp_bcast(nchain,1,rootid,groupid)
    call mp_bcast(nmax,1,rootid,groupid)
    call mp_bcast(nstep,1,rootid,groupid)
    call mp_bcast(time_limit,1,rootid,groupid)
    call mp_bcast(lstop,1,rootid,groupid)
    call mp_bcast(iratio,1,rootid,groupid)
    call mp_bcast(rmin,1,rootid,groupid)
    call mp_bcast(softcut,1,rootid,groupid)
    call mp_bcast(linit,1,rootid,groupid)
    call mp_bcast(lreadq,1,rootid,groupid)
    call mp_bcast(N_add,1,rootid,groupid)
    call mp_bcast(box2add,1,rootid,groupid)
    call mp_bcast(moltyp2add,1,rootid,groupid)

    nbxmax=nbox+1
    npabmax=nbxmax*(nbxmax+1)/2
    ntmax=nmolty+1
    npamax=ntmax*(ntmax-1)/2
    nprop=nEnergy+(4*ntmax)+7
    if (nmax<nchain+2) nmax=nchain+2
    if (N_add.gt.0) nmax=nmax+N_add
    softlog = 10.0_dp**(-softcut)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_SHARED','------------------------------------------'
       write(io_output,'(A,I0)') 'Random number seed: ',seed
       write(io_output,'(A,I0)') 'number of boxes in the system: ',nbox
       write(io_output,'(A,I0)') 'number of molecule types: ',nmolty
       write(io_output,'(A,I0)') 'number of chains: ',nchain
       if (lstop) then
          write(io_output,'(A,I0)') 'number of steps: ',nstep
       else
          write(io_output,'(A,I0)') 'number of cycles: ',nstep
       end if
       if (time_limit.gt.0) then
          write(io_output,'(A,I0,A)') 'the code will try to run for ',time_limit,' seconds'
       end if
       w_i(iratio)
       write(io_output,'(A,F7.3,A)') 'minimum cutoff (rmin): ',rmin,' [Ang]'
       w_r(softcut)
       write(io_output,'(A,I0,A,I0,A)') 'Write checkpoint file every ',checkpoint_interval,' seconds, and keep the last '&
        ,checkpoint_copies,' copies'
       w_l(linit)
       w_l(lreadq)
    end if

    ! initialize random number generator
    call ranset(seed,numprocs)

    call allocate_system()
    call allocate_sim_cell()
    call allocate_kspace()
    call allocate_kdtree()
! -------------------------------------------------------------------
    !> read name list analysis
    virtemp=0.0_dp

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=analysis,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: analysis',jerr)
    end if

    call mp_bcast(iprint,1,rootid,groupid)
    call mp_bcast(imv,1,rootid,groupid)
    call mp_bcast(iblock,1,rootid,groupid)
    call mp_bcast(iratp,1,rootid,groupid)
    call mp_bcast(idiele,1,rootid,groupid)
    call mp_bcast(iheatcapacity,1,rootid,groupid)
    call mp_bcast(ianalyze,1,rootid,groupid)
    call mp_bcast(nbin,1,rootid,groupid)
    call mp_bcast(lrdf,1,rootid,groupid)
    call mp_bcast(lintra,1,rootid,groupid)
    call mp_bcast(lstretch,1,rootid,groupid)
    call mp_bcast(lgvst,1,rootid,groupid)
    call mp_bcast(lbend,1,rootid,groupid)
    call mp_bcast(lete,1,rootid,groupid)
    call mp_bcast(lrhoz,1,rootid,groupid)
    call mp_bcast(bin_width,1,rootid,groupid)
    call mp_bcast(lucall,1,rootid,groupid)
    call mp_bcast(nvirial,1,rootid,groupid)
    call mp_bcast(starvir,1,rootid,groupid)
    call mp_bcast(stepvir,1,rootid,groupid)
    call mp_bcast(ntemp,1,rootid,groupid)
    call mp_bcast(virtemp,1,rootid,groupid)

    blockm=nstep/iblock
    if (lucall) call read_prop_widom(io_input,lprint,blockm)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST ANALYSIS','------------------------------------------'
       w_i(iprint)
       w_i(imv)
       w_i(iblock)
       w_i(iratp)
       w_i(idiele)
       w_i(iheatcapacity)
       w_i(ianalyze)
       w_i(nbin)
       w_l(lrdf)
       w_l(lintra)
       w_l(lstretch)
       w_l(lgvst)
       w_l(lbend)
       w_l(lete)
       w_l(lrhoz)
       w_r(bin_width)
       w_l(lucall)
       if (lvirial) then
          w_i(nvirial)
          w_r(starvir)
          w_r(stepvir)
          write(io_output,*)&
           '=> Virial coefficient will be calculated at the following temperature:',virtemp
       end if
    end if

    if (lvirial) then
       if (nvirial.gt.maxvir) call err_exit(__FILE__,__LINE__,'nvirial .gt. maxvir',myid+1)
       if (nchain.ne.2) call err_exit(__FILE__,__LINE__,'nchain must equal 2',myid+1)
    end if

    if (allocated(io_box_movie)) deallocate(io_box_movie,stat=jerr)
    if (allocated(io_box_movie_pdb)) deallocate(io_box_movie_pdb,stat=jerr)
    allocate(lhere(nntype),temphe(nntype),io_box_movie(nbxmax),io_box_movie_pdb(nbxmax),ncarbon(ntmax),idummy(ntmax)&
     ,qbox(nbxmax),nures(ntmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readdat: allocating system failed',jerr)
    lhere=.false.
! -------------------------------------------------------------------
    !> read parameters about simulation boxes

    ! Looking for section SIMULATION_BOX
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_CELL:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section SIMULATION_BOX not found',jerr)

          if (UPPERCASE(line_in(1:14)).eq.'SIMULATION_BOX') then
             if (lprint) then
                write(io_output,'(/,A,/,A)') 'SECTION SIMULATION_BOX','------------------------------------------'
             end if
             exit cycle_read_cell
          end if
       END DO CYCLE_READ_CELL
    end if

    do i=1,nbox+1
       if (myid.eq.rootid) then
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SIMULATION_BOX',jerr)
          if (UPPERCASE(line_in(1:18)).eq.'END SIMULATION_BOX') then
             if (i.ne.nbox+1) call err_exit(__FILE__,__LINE__,'Section SIMULATION_BOX not complete!',jerr)
             exit
          else if (i.eq.nbox+1) then
             call err_exit(__FILE__,__LINE__,'Section SIMULATION_BOX has more than nbox records!',jerr)
          end if

          ! boxlx boxly boxlz rcut kalp rcutnn numDimensionIsIstropic lsolid lrect lideal ltwice temperature pressure
          !> provision for different temperatures (just like different pressures) in each box, such as for parallel tempering
          read(line_in,*) boxlx(i),boxly(i),boxlz(i),rcut(i),kalp(i),rcutnn(i),numberDimensionIsIsotropic(i),lsolid(i)&
           ,lrect(i),lideal(i),ltwice(i),rtmp,express(i)

          if ((temp.ge.0).and.(i.gt.1)) then ! temp is initialized in sim_system.F90 to -1.0_dp
             if (temp.ne.rtmp) call err_exit(__FILE__,__LINE__,'Section SIMULATION_BOX: temperature not the same for each box'&
              ,myid+1)
          else
             temp=rtmp
          end if

          if (lideal(i).and.lexpee) call err_exit(__FILE__,__LINE__&
           ,'Cannot have lideal and lexpee both true (if you want this you will have change code)',myid+1)
       else if (i.eq.nbox+1) then
          exit
       end if

       if (i.eq.1 .and. lexzeo) then
          ! load positions of zeolite atoms
          call zeocoord(file_in,lprint)
       end if

       call mp_bcast(boxlx(i),1,rootid,groupid)
       call mp_bcast(boxly(i),1,rootid,groupid)
       call mp_bcast(boxlz(i),1,rootid,groupid)
       call mp_bcast(rcut(i),1,rootid,groupid)
       call mp_bcast(kalp(i),1,rootid,groupid)
       call mp_bcast(rcutnn(i),1,rootid,groupid)
       call mp_bcast(numberDimensionIsIsotropic(i),1,rootid,groupid)
       call mp_bcast(lsolid(i),1,rootid,groupid)
       call mp_bcast(lrect(i),1,rootid,groupid)
       call mp_bcast(lideal(i),1,rootid,groupid)
       call mp_bcast(ltwice(i),1,rootid,groupid)
       call mp_bcast(express(i),1,rootid,groupid)

       if (lprint) then
          write(io_output,'(A,I0,A,2(F8.3," x "),F8.3)') 'Box ',i,': ',boxlx(i),boxly(i),boxlz(i)
          write(io_output,'(2(A,F6.3))') '   rcut: ',rcut(i),' [Ang], kalp: ',kalp(i)
          write(io_output,'(A,F6.3)') '   neighbor list cutoff (rcutnn): ',rcutnn(i)
          write(io_output,'(A,I0)') '   number of dimensions that are isotropic: ',numberDimensionIsIsotropic(i)
          write(io_output,'(4(A,L2))') '   lsolid: ',lsolid(i),', lrect: ',lrect(i),', lideal: ',lideal(i),', ltwice: '&
           ,ltwice(i)
          write(io_output,'(A,F8.3,A)') '   temperature: ',rtmp,' [K]'
          write(io_output,'(A,G16.9,A)') '   external pressure: ',express(i),' [MPa]'
       end if

       if (myid.eq.rootid) then
          ! nchain_1 ... nchain_nmolty ghost_particles
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SIMULATION_BOX',jerr)
          read(line_in,*) (ininch(j,i),j=1,nmolty),ghost_particles(i)

          if (lprint) then
             write(io_output,'(A,'//format_n(nmolty,'(2X,I0)')//')') '   initial number of chains of each type: '&
              ,ininch(1:nmolty,i)
             write(io_output,'(A,I0)') '   Ghost particles: ',ghost_particles(i)
          end if

          ! inix iniy iniz inirot inimix zshift dshift use_linkcell rintramax
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SIMULATION_BOX',jerr)
          read(line_in,*) inix(i),iniy(i),iniz(i),inirot(i),inimix(i),zshift(i),dshift(i),ltmp,rtmp

          if (lprint) then
             write(io_output,'(A,2(I0," x "),I0)') '   initial number of chains in x, y and z directions: ',inix(i),iniy(i)&
              ,iniz(i)
             write(io_output,'(2(A,I0),A,F5.1,A,F7.3)') '   initial rotational displacement: ',inirot(i),', inimix: '&
              ,inimix(i),', zshift: ',zshift(i),', dshift: ',dshift(i)
             if (ltmp) write(io_output,'(A,F8.3)') '   Linked cell structures will be used for this box, rintramax = ',rtmp
          end if

          ! read linkcell information
          if (ltmp) then
             if (licell) call err_exit(__FILE__,__LINE__,'Section SIMULATION_BOX: link cell only allowed for one box',myid+1)
             licell=.true.
             boxlink=i
             rintramax=rtmp
             if (lsolid(boxlink).and.(.not.lrect(boxlink))) call err_exit(__FILE__,__LINE__&
              ,'Linkcell not implemented for nonrectangular boxes',myid+1)
          end if
       end if

       call mp_bcast(ininch(:,i),nmolty,rootid,groupid)
       call mp_bcast(ghost_particles(i),1,rootid,groupid)
       call mp_bcast(inix(i),1,rootid,groupid)
       call mp_bcast(iniy(i),1,rootid,groupid)
       call mp_bcast(iniz(i),1,rootid,groupid)
       call mp_bcast(inirot(i),1,rootid,groupid)
       call mp_bcast(inimix(i),1,rootid,groupid)
       call mp_bcast(zshift(i),1,rootid,groupid)
       call mp_bcast(dshift(i),1,rootid,groupid)
    end do

    call mp_bcast(temp,1,rootid,groupid)
    beta = 1.0_dp / temp
    call mp_bcast(licell,1,rootid,groupid)
    call mp_bcast(boxlink,1,rootid,groupid)
    call mp_bcast(rintramax,1,rootid,groupid)

    temnc = 0
    do i=1,nmolty
       temtyp(i)=sum(ininch(i,1:nbox))
       do j = 1, temtyp(i)
          temnc = temnc + 1
          moltyp(temnc) = i
       end do
    end do

    if (lprint) then
       write(io_output,'(/,A)') 'NUMBER OF MOLECULES OF EACH TYPE'
       write(io_output,'(A,'//format_n(nmolty,'(2X,I0)')//')') ' number of chains of each type: ',temtyp(1:nmolty)
       if (sum(temtyp(1:nmolty)).ne.nchain) then
          do j = 1,nbox
             write(io_output,*) 'ibox:',j,', ininch:',(ininch(i,j),i=1,nmolty)
          end do
          write(io_output,*) 'nchain:',nchain
          call err_exit(__FILE__,__LINE__,'inconsistant number of chains',myid+1)
       end if
    end if

    express = express*MPa2SimUnits

    if (lshift.or.lsami) then
       ! Keep the rcut same for each box
       do ibox = 2,nbox
          if (abs(rcut(1)-rcut(ibox)).gt.1.0E-10_dp) then
             call err_exit(__FILE__,__LINE__,'Keep rcut for each box same',myid+1)
          end if
       end do
    end if
! -------------------------------------------------------------------
    call allocate_neighbor_list()
! -------------------------------------------------------------------
    !> read parameters about molecule types
    !> The division of variables between here and in the
    !> section specific to a particular MC move is somewhat
    !> arbitrary. Variables that are "intrinsic" to the type
    !> of molecules are place here, such as whether it's
    !> rigid (lrigid), is a ring (lring), etc, while
    !> variables that are model- or algorithm-specific are
    !> placed in the corresponding MC move sections, such as
    !> whether the molecule has fluctuating charges
    !> (lflucq), need to treat with expanded ensemble moves
    !> (lexpand), etc.

    ! Looking for section MOLECULE_TYPE
    nugrow=0
    nunit=0
    numax=0
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_MOLTYP:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section MOLECULE_TYPE not found',jerr)

          if (UPPERCASE(line_in(1:13)).eq.'MOLECULE_TYPE') then
             if (lprint) then
                write(io_output,'(/,A,/,A)') 'SECTION MOLECULE_TYPE','------------------------------------------'
             end if
             ntype=0

             do imol=1,nmolty+1
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                if (UPPERCASE(line_in(1:17)).eq.'END MOLECULE_TYPE') then
                   if (imol.ne.nmolty+1) call err_exit(__FILE__,__LINE__,'Section MOLECULE_TYPE not complete!',jerr)
                   exit
                else if (imol.eq.nmolty+1) then
                   call err_exit(__FILE__,__LINE__,'Section MOLECULE_TYPE has more than nmolty records!',jerr)
                end if

                if (scan(line_in(1:10),'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ') >= 1) then
                ! nunit nugrow ncarbon maxcbmc maxgrow iring lelect lring lrigid lbranch lsetup lq14scale qscale iurot isolute
                    read(line_in,*) molecname(imol),nunit(imol),nugrow(imol),ncarbon(imol),nmaxcbmc(imol),maxgrow(imol),iring(imol)&
               ,lelect(imol),lring(imol),lrigid(imol),lbranch(imol),lsetup,lq14scale(imol),qscale(imol),iurot(imol),isolute(imol)
                else
                    read(line_in,*) nunit(imol),nugrow(imol),ncarbon(imol),nmaxcbmc(imol),maxgrow(imol),iring(imol),lelect(imol)&
                    ,lring(imol),lrigid(imol),lbranch(imol),lsetup,lq14scale(imol),qscale(imol),iurot(imol),isolute(imol)
                    molecname(imol)='undefined'
                end if
                if (lprint) then
                   write(io_output,'(A,I2,A,A10)') 'molecule type: ',imol,' ',molecname(imol)
                   write(io_output,'(A,I0)') '   number of units: ',nunit(imol)
                   write(io_output,'(A,I0)') '   number of units for CBMC growth: ', nugrow(imol)
                   write(io_output,'(A,I0)') '   number of carbons for EH alkane: ', ncarbon(imol)
                   write(io_output,'(A,I0)') '   maximum number of units for CBMC: ', nmaxcbmc(imol)
                   write(io_output,'(A,I0)') '   maximum number of interior segments for SAFE-CBMC regrowth: ',maxgrow(imol)
                   write(io_output,'(A,I0)') '   number of atoms in a ring (if lring=.true.): ',iring(imol)
                   write(io_output,'(2(A,I0),6(A,L2),A,F3.1)') '   iurot: ',iurot(imol),', isolute: ',isolute(imol),', lelect: '&
                    ,lelect(imol),', lring: ',lring(imol),', lrigid: ',lrigid(imol),', lbranch: ',lbranch(imol)&
                    ,', lsetup: ',lsetup,', lq14scale: ',lq14scale(imol),', qscale: ',qscale(imol)
                end if

                if (nunit(imol).gt.numax) then
                   numax=nunit(imol)
                   if (numax.gt.ubound(ntype,2)) then
                      call reallocate(ntype,1,ntmax,1,2*numax)
                      call reallocate(riutry,1,ntmax,1,2*numax)
                      call reallocate(leaderq,1,ntmax,1,2*numax)
                      call reallocate(invib,1,ntmax,1,2*numax)
                      call reallocate(itvib,1,ntmax,1,2*numax,1,6)
                      call reallocate(ijvib,1,ntmax,1,2*numax,1,6)
                      call reallocate(inben,1,ntmax,1,2*numax)
                      call reallocate(itben,1,ntmax,1,2*numax,1,12)
                      call reallocate(ijben2,1,ntmax,1,2*numax,1,12)
                      call reallocate(ijben3,1,ntmax,1,2*numax,1,12)
                      call reallocate(intor,1,ntmax,1,2*numax)
                      call reallocate(ittor,1,ntmax,1,2*numax,1,12)
                      call reallocate(ijtor2,1,ntmax,1,2*numax,1,12)
                      call reallocate(ijtor3,1,ntmax,1,2*numax,1,12)
                      call reallocate(ijtor4,1,ntmax,1,2*numax,1,12)
                      call reallocate(irotbd,1,2*numax,1,ntmax)
                      call reallocate(pmrotbd,1,2*numax,1,ntmax)
                   end if
                end if

                if (lrigid(imol)) then
                   ! n, growpoint_1 ... growpoint_n
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                   read(line_in,*) rindex(imol),(riutry(imol,i),i=1,rindex(imol))
                   if (rindex(imol).le.0) then
                      riutry(imol,1) = 1
                   else if (lprint) then
                      write(io_output,'(A,I0,A,'//format_n(rindex(imol),'(1X,I0)')//')') '   The following ',rindex(imol)&
                       ,' beads have flexible side chains:',(riutry(imol,i),i=1,rindex(imol))
                   end if
                end if

                if (lsetup) then
                   call molsetup(io_input,imol,lprint)
                   goto 112
                end if

                do i = 1, nunit(imol)
                   ! linear/branched chain with connectivity table
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                   ! unit ntype leaderq
                   read(line_in,*) j,ntype(imol,i),leaderq(imol,i)
                   if (leaderq(imol,i).gt.j.and..not.lchgall) call err_exit(__FILE__,__LINE__,'group-based cut-off screwed for qq'&
                    ,myid+1)
                   ntype(imol,i)=indexOf(atoms,ntype(imol,i))

                   if (lprint) then
                      write(io_output,'(/,2(A,I0),A,A,A,I0)') '   bead ',i,': bead type ',atoms%list(ntype(imol,i)),' ['&
                       ,trim(chemid(ntype(imol,i))),'], charge leader ',leaderq(imol,i)
                   end if

                   if (ntype(imol,i).eq.0) then
                      call err_exit(__FILE__,__LINE__,'ERROR: atom type undefined!',myid+1)
                   end if

                   ! bond stretching
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                   read(line_in,*) invib(imol,i)
                   if (invib(imol,i).gt.nvib_max) then
                      write(io_output,*) 'imol',imol,'   i',i,'  invib',invib(imol,i)
                      call err_exit(__FILE__,__LINE__,'too many vibrations',myid+1)
                   end if

                   do j = 1, invib(imol,i)
                      call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                      if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                      read(line_in,*) ijvib(imol,i,j),itvib(imol,i,j)
                      if((ijvib(imol,i,j).eq.i).or.(ijvib(imol,i,j).gt.nunit(imol))) then
                         call err_exit(__FILE__,__LINE__,'check vibrations for mol type '//integer_to_string(imol)//' and bead '&
                          //integer_to_string(i),myid+1)
                      end if

                      itvib(imol,i,j)=indexOf(bonds,itvib(imol,i,j))
                      if (itvib(imol,i,j).eq.0) then
                         call err_exit(__FILE__,__LINE__,'ERROR: stretching parameters undefined!',myid+1)
                      else if (brvib(itvib(imol,i,j)).lt.1E-06_dp) then
                         call err_exit(__FILE__,__LINE__,'ERROR: stretching parameters undefined!',myid+1)
                      end if

                      if (lprint) then
                         write(io_output,'(2(A,I0),A,F8.5,A,G16.9)') '      bonded to bead ',ijvib(imol,i,j),', type '&
                          ,bonds%list(itvib(imol,i,j)),', bond length: ',brvib(itvib(imol,i,j)),', k/2: ',brvibk(itvib(imol,i,j))
                      end if
                   end do

                   ! bond bending -
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                   read(line_in,*) inben(imol,i)
                   if (inben(imol,i).gt.nben_max) then
                      call err_exit(__FILE__,__LINE__,'too many bends',myid+1)
                   end if

                   do j = 1, inben(imol,i)
                      call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                      if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                      read(line_in,*) ijben2(imol,i,j),ijben3(imol,i,j),itben(imol,i,j)
                      if ((ijben2(imol,i,j).gt.nunit(imol)).or.(ijben3(imol,i,j).gt.nunit(imol))) then
                         call err_exit(__FILE__,__LINE__,'check bending for molecule type '//integer_to_string(imol)//' bead '&
                          //integer_to_string(i),myid+1)
                      end if
                      if ((ijben2(imol,i,j).eq.i).or.(ijben3(imol,i,j).eq.i).or.(ijben2(imol,i,j).eq.ijben3(imol,i,j))) then
                         call err_exit(__FILE__,__LINE__,'check bending for molecule type '//integer_to_string(imol)//' bead '&
                          //integer_to_string(i),myid+1)
                      end if

                      itben(imol,i,j)=indexOf(angles,itben(imol,i,j))
                      if (itben(imol,i,j).eq.0) then
                         call err_exit(__FILE__,__LINE__,'ERROR: bending parameters undefined!',myid+1)
                      else if (brben(itben(imol,i,j)).lt.1E-06_dp) then
                         call err_exit(__FILE__,__LINE__,'ERROR: bending parameters undefined!',myid+1)
                      end if

                      if (lprint) then
                         write(io_output,'(3(A,I0),A,F8.3,A,G16.9)') '      bending interaction through ',ijben2(imol,i,j)&
                          ,' with bead ',ijben3(imol,i,j),', bending type: ',angles%list(itben(imol,i,j)),', bending angle: '&
                          ,brben(itben(imol,i,j))*raddeg,', k/2: ',brbenk(itben(imol,i,j))
                      end if
                   end do

                   ! bond torsion -
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                   read(line_in,*) intor(imol,i)
                   if (intor(imol,i).gt.ntor_max) then
                      call err_exit(__FILE__,__LINE__,'too many torsions',myid+1)
                   end if

                   do j = 1, intor(imol,i)
                      call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                      if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                      read(line_in,*) ijtor2(imol,i,j),ijtor3(imol,i,j),ijtor4(imol,i,j),ittor(imol,i,j)

                      if (ijtor2(imol,i,j).gt.nunit(imol).or.ijtor3(imol,i,j).gt.nunit(imol).or.ijtor4(imol,i,j).gt.nunit(imol)) then
                         call err_exit(__FILE__,__LINE__,'check torsion for molecule type '//integer_to_string(imol)//' bead '&
                          //integer_to_string(i),myid+1)
                      end if
                      if((ijtor2(imol,i,j).eq.i .or. ijtor3(imol,i,j).eq.i .or. ijtor4(imol,i,j).eq.i)&
                       .or. (ijtor2(imol,i,j).eq.ijtor3(imol,i,j) .or. ijtor2(imol,i,j).eq.(ijtor4(imol,i,j))&
                       .or. (ijtor3(imol,i,j).eq.ijtor4(imol,i,j)))) then
                         call err_exit(__FILE__,__LINE__,'check torsion for molecule type '//integer_to_string(imol)//' bead '&
                          //integer_to_string(i),myid+1)
                      end if

                      ittor(imol,i,j)=indexOf(dihedrals,ittor(imol,i,j))
                      if (ittor(imol,i,j).eq.0) then
                         call err_exit(__FILE__,__LINE__,'ERROR: torsion parameters undefined!',myid+1)
                      end if

                      if (lprint) then
                         write(io_output,'(4(A,I0))') '      torsional interaction through ',ijtor2(imol,i,j),' and '&
                          ,ijtor3(imol,i,j),' with bead ',ijtor4(imol,i,j),', torsional type: ',dihedrals%list(ittor(imol,i,j))
                      end if
                   end do
                end do

   112          continue

                !kea 6/4/09 -- added for multiple rotation centers
                ! To assign multiple rotation centers, set iurot(imol) < 0
                ! Add line after molecule specification, avbmc parameters
                ! First, number of rotation centers
                ! Second, identity of centers (0=COM,integer::> 0 = bead number)
                ! Third, give probability to rotate around different centers
                if(iurot(imol).lt.0) then
                   call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MOLECULE_TYPE',jerr)
                   read(line_in,*) nrotbd(imol),irotbd(1:nrotbd(imol),imol),pmrotbd(1:nrotbd(imol),imol)
                   if (lprint) then
                      write(io_output,'(I0,A)') nrotbd(imol),' rotation centers:'
                      do ii=1,nrotbd(imol)
                         write(io_output,'(A,I0,A,F5.3)') '   around bead: ',irotbd(ii,imol),', probability: ',pmrotbd(ii,imol)
                      end do
                   end if
                end if

                ! starting the self consistency check for the bond vibrations, bending, and torsions
                ! this would help in catching errors in fort.4 connectivity. Starting after continue
                ! so that if we use molsetup subroutine it will provide extra checking. (Neeraj)
                if (.not.lrigid(imol)) then
                   do i = 1,nunit(imol)
                      numvib=invib(imol,i)
                      numbend=inben(imol,i)
                      numtor=intor(imol,i)
                      lfound = .false.
                      do j = 1,numvib
                         vib1 = ijvib(imol,i,j)
                         vibtype  = itvib(imol,i,j)
                         if(invib(imol,vib1).eq.0) then
                            write(io_output,*) 'Check vibration for mol. type:',imol, 'bead',vib1,'with',i
                            call err_exit(__FILE__,__LINE__,'ERROR IN FORT.4 VIBRATIONS',myid+1)
                         end if
                         do k =1,invib(imol,vib1)
                            if(ijvib(imol,vib1,k).eq.i) then
                               lfound = .true.
                               if(vibtype.ne.itvib(imol,vib1,k)) then
                                  write(io_output,*) 'check vibration type of bead',i, 'with',vib1,'molecule type',imol,'vice versa'
                                  call err_exit(__FILE__,__LINE__,'Error in fort.4 vibration specifications',myid+1)
                               end if
                            end if
                         end do
                         if(.not.lfound) then
                            write(io_output,*) 'Check vibration for mol. type:',imol, 'bead ',vib1,'with ',i
                            call err_exit(__FILE__,__LINE__,'Error in fort.4 vibration iformation',myid+1)
                         end if
                      end do
                      lfound= .false.
                      do j = 1,numbend
                         bend2 = ijben2(imol,i,j)
                         bend3 = ijben3(imol,i,j)
                         bendtype = itben(imol,i,j)
                         if(inben(imol,bend3).eq.0) then
                            write(io_output,*) 'Check bending for mol. type:',imol, 'bead ',bend3,'with ',i
                            call err_exit(__FILE__,__LINE__,'ERROR IN FORT.4 BENDING',myid+1)
                         end if
                         do k = 1,inben(imol,bend3)
                            if((ijben2(imol,bend3,k).eq.bend2).and. (ijben3(imol,bend3,k).eq.i)) then
                               lfound = .true.
                               if(itben(imol,bend3,k).ne.bendtype) then
                                  write(io_output,*) 'check bending type of bead',i, 'with',bend3,'mol. typ.',imol,'and vice versa'
                                  call err_exit(__FILE__,__LINE__,'Error in fort.4 bending specifications',myid+1)
                               end if
                            end if
                         end do
                         if(.not.lfound) then
                            write(io_output,*) 'Check bending for mol. type:',imol, 'bead ',bend3,'with ',i
                            call err_exit(__FILE__,__LINE__,'Error in fort.4 bending information',myid+1)
                         end if
                      end do
                      lfound = .false.
                      do j = 1,numtor
                         tor2 = ijtor2(imol,i,j)
                         tor3 = ijtor3(imol,i,j)
                         tor4 = ijtor4(imol,i,j)
                         tortype = ittor(imol,i,j)
                         if(intor(imol,tor4).eq.0) then
                            write(io_output,*) 'Check torsion for mol. type:',imol, 'bead ',tor4,'with ',i,'and vice versa'
                            call err_exit(__FILE__,__LINE__,'ERROR IN FORT.4 TORSION',myid+1)
                         end if
                         do k = 1,intor(imol,tor4)
                            if((ijtor2(imol,tor4,k).eq.tor3).and.(ijtor3(imol,tor4 ,k).eq.tor2).and.(ijtor4(imol,tor4,k).eq.i)) then
                               lfound=.true.
                               if(ittor(imol,tor4,k).ne.tortype) then
                                  write(io_output,*) 'check torsion type of bead',i, 'with',tor4,'mol. typ.',imol,'and vice versa'
                                  call err_exit(__FILE__,__LINE__,'Error in fort.4 torsion specifications',myid+1)
                               end if
                            end if
                         end do
                         if(.not.lfound) then
                            write(io_output,*) 'Check torsion for mol. type:',imol, 'bead ',tor4,'with ',i
                            call err_exit(__FILE__,__LINE__,'Error in fort.4 torsion information',myid+1)
                         end if
                      end do
                   end do
                end if

                ! Neeraj Adding molecule neutrality check
                if (.not.(lgaro.or.lionic.or.lexzeo)) then
                   qtot =0.0E0_dp
                   do j = 1,nunit(imol)
                      qtot = qtot+qelect(ntype(imol,j))
                   end do
                   if(abs(qtot).gt.1E-7_dp) call err_exit(__FILE__,__LINE__,'molecule type '//integer_to_string(imol)&
                    //' not neutral. check charges',myid+1)
                end if
             end do

             if (ALL(.NOT.lelect(1:nmolty))) then
                if (lewald.or.lchgall) then
                   if (lprint) write(io_output,'(A)') 'No charges in the system -> lewald is now turned off'
                end if
             end if

             exit cycle_read_moltyp
          end if
       END DO CYCLE_READ_MOLTYP
    end if

    call mp_bcast(numax,1,rootid,groupid)

    if (numax.gt.0) then
       if (allocated(lplace)) deallocate(lplace,lrigi,stat=jerr)
       allocate(lplace(ntmax,numax),lrigi(ntmax,numax),inclmol(ntmax*numax*numax),inclbead(ntmax*numax*numax,2)&
        ,inclsign(ntmax*numax*numax),ofscale(ntmax*numax*numax),ofscale2(ntmax*numax*numax),ainclmol(ntmax*numax*numax)&
        ,ainclbead(ntmax*numax*numax,2),a15t(ntmax*numax*numax),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readdat: allocating molecule failed',jerr)
       lrigi = .false.

       call reallocate(ntype,1,ntmax,1,numax)
       call reallocate(riutry,1,ntmax,1,numax)
       call reallocate(leaderq,1,ntmax,1,numax)
       call reallocate(invib,1,ntmax,1,numax)
       call reallocate(itvib,1,ntmax,1,numax,1,6)
       call reallocate(ijvib,1,ntmax,1,numax,1,6)
       call reallocate(inben,1,ntmax,1,numax)
       call reallocate(itben,1,ntmax,1,numax,1,12)
       call reallocate(ijben2,1,ntmax,1,numax,1,12)
       call reallocate(ijben3,1,ntmax,1,numax,1,12)
       call reallocate(intor,1,ntmax,1,numax)
       call reallocate(ittor,1,ntmax,1,numax,1,12)
       call reallocate(ijtor2,1,ntmax,1,numax,1,12)
       call reallocate(ijtor3,1,ntmax,1,numax,1,12)
       call reallocate(ijtor4,1,ntmax,1,numax,1,12)
       call reallocate(irotbd,1,numax,1,ntmax)
       call reallocate(pmrotbd,1,numax,1,ntmax)

       call mp_bcast(nunit,nmolty,rootid,groupid)
       call mp_bcast(nugrow,nmolty,rootid,groupid)
       call mp_bcast(ncarbon,nmolty,rootid,groupid)
       call mp_bcast(nmaxcbmc,nmolty,rootid,groupid)
       call mp_bcast(maxgrow,nmolty,rootid,groupid)
       call mp_bcast(iring,nmolty,rootid,groupid)
       call mp_bcast(lelect,nmolty,rootid,groupid)
       call mp_bcast(lring,nmolty,rootid,groupid)
       call mp_bcast(lrigid,nmolty,rootid,groupid)
       call mp_bcast(lbranch,nmolty,rootid,groupid)
       call mp_bcast(lq14scale,nmolty,rootid,groupid)
       call mp_bcast(qscale,nmolty,rootid,groupid)
       call mp_bcast(iurot,nmolty,rootid,groupid)
       call mp_bcast(isolute,nmolty,rootid,groupid)
       call mp_bcast(rindex,nmolty,rootid,groupid)
       call mp_bcast(riutry,ntmax*numax,rootid,groupid)
       call mp_bcast(ntype,ntmax*numax,rootid,groupid)
       call mp_bcast(leaderq,ntmax*numax,rootid,groupid)
       call mp_bcast(invib,ntmax*numax,rootid,groupid)
       call mp_bcast(itvib,ntmax*numax*nvib_max,rootid,groupid)
       call mp_bcast(ijvib,ntmax*numax*nvib_max,rootid,groupid)
       call mp_bcast(inben,ntmax*numax,rootid,groupid)
       call mp_bcast(itben,ntmax*numax*nben_max,rootid,groupid)
       call mp_bcast(ijben2,ntmax*numax*nben_max,rootid,groupid)
       call mp_bcast(ijben3,ntmax*numax*nben_max,rootid,groupid)
       call mp_bcast(intor,ntmax*numax,rootid,groupid)
       call mp_bcast(ittor,ntmax*numax*ntor_max,rootid,groupid)
       call mp_bcast(ijtor2,ntmax*numax*ntor_max,rootid,groupid)
       call mp_bcast(ijtor3,ntmax*numax*ntor_max,rootid,groupid)
       call mp_bcast(ijtor4,ntmax*numax*ntor_max,rootid,groupid)
       call mp_bcast(nrotbd,nmolty,rootid,groupid)
       call mp_bcast(irotbd,numax*ntmax,rootid,groupid)
       call mp_bcast(pmrotbd,numax*ntmax,rootid,groupid)

       do imol = 1,nmolty
          masst(imol) = 0.0E0_dp

          do i = 1,nunit(imol)
             iutemp = ntype(imol,i)

             if (lpl(iutemp)) then
                lplace(imol,i) = .true.
                llplace(imol) = .true.
             else
                lplace(imol,i) = .false.
             end if

             masst(imol)=masst(imol)+mass(iutemp)
             lhere(iutemp) = .true.
          end do
       end do
    end if

    if (lprint) then
       write(io_output,'(/,A,'//format_n(nmolty,'(3X,F10.5)')//')') 'MOLECULAR MASS: ',masst(1:nmolty)
    end if
! -------------------------------------------------------------------
    call allocate_molecule()
    call allocate_energy_bonded()

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'SPECIAL INTERACTION RULES','------------------------------------------'
    end if

    !> Looking for section INTERMOLECULAR_EXCLUSION
    ! read exclusion table for intermolecular interactions
    nexclu=0
    lexclu=.false.
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_EXCLUSION:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_exclusion
          if (UPPERCASE(line_in(1:24)).eq.'INTERMOLECULAR_EXCLUSION') then
             jerr=0
             exit cycle_read_exclusion
          end if
       END DO CYCLE_READ_EXCLUSION
    end if

    call mp_bcast(jerr,1,rootid,groupid)

    if (jerr.eq.0) then
       do
          if (myid.eq.rootid) then
             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section INTERMOLECULAR_EXCLUSION',jerr)
          end if

          call mp_bcast(line_in,rootid,groupid)

          if (UPPERCASE(line_in(1:28)).eq.'END INTERMOLECULAR_EXCLUSION') exit

          nexclu=nexclu+1
          ! mol_1 bead_1 mol_2 bead_2
          read(line_in,*) i,ii,j,jj

          if (lprint) then
             write(io_output,'(4(A,I0))') '      excluding interactions between bead ',ii,' of molecule ',i,' with bead '&
              ,jj,' of molecule ',j
          end if

          lexclu(i,ii,j,jj) = .true.
          lexclu(j,jj,i,ii) = .true.
       end do

       if (lprint) then
          write(io_output,'(A,I0,A,/)') '  Total: ',nexclu,' exclusion rules for intermolecular interactions'
       end if
    end if

! -------------------------------------------------------------------
    !> Looking for section ZEOLITE_EXCLUSION
    ! read exclusion table for intermolecular interactions
    lexclu_zeo=.false.
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_ZEO_EXCLUSION:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit CYCLE_READ_ZEO_EXCLUSION
          if (UPPERCASE(line_in(1:17)).eq.'ZEOLITE_EXCLUSION') then
             jerr=0
             exit CYCLE_READ_ZEO_EXCLUSION
          end if
       END DO CYCLE_READ_ZEO_EXCLUSION
    end if

    call mp_bcast(jerr,1,rootid,groupid)

    if (jerr.eq.0) then
       do
          if (myid.eq.rootid) then
             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section ZEOLITE_EXCLUSION',jerr)
          end if

          call mp_bcast(line_in,rootid,groupid)

          if (UPPERCASE(line_in(1:21)).eq.'END ZEOLITE_EXCLUSION') exit

          ! mol_1 bead_1 mol_2 bead_2
          read(line_in,*) i

          if (lprint) then
             write(io_output,'(4(A,I0))') '      excluding interactions between molecule ',i,' and zeolite'
          end if

          lexclu_zeo(i) = .true.
       end do

    end if

! -------------------------------------------------------------------
    !> Looking for section INTRAMOLECULAR_SPECIAL
    ! read exclusion table for intermolecular interactions
    inclnum=0
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_SPECIAL:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_special

          if (UPPERCASE(line_in(1:22)).eq.'INTRAMOLECULAR_SPECIAL') then
             jerr=0
             exit cycle_read_special
          end if
       END DO CYCLE_READ_SPECIAL
    end if

    call mp_bcast(jerr,1,rootid,groupid)

    if (jerr.eq.0) then
       do
          if (myid.eq.rootid) then
             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section INTRAMOLECULAR_SPECIAL',jerr)
          end if

          call mp_bcast(line_in,rootid,groupid)

          if (UPPERCASE(line_in(1:26)).eq.'END INTRAMOLECULAR_SPECIAL') exit

          inclnum=inclnum+1
          ! inclmol inclbead_1 inclbead_2 inclsign ofscale ofscale2
          read(line_in,*) inclmol(inclnum),inclbead(inclnum,1),inclbead(inclnum,2),inclsign(inclnum),ofscale(inclnum)&
           ,ofscale2(inclnum)

          if (lprint) then
             if (inclsign(inclnum) .eq. 1) then
                write(io_output,'(3(A,I0),2(A,F6.3))') '      including intramolecular interactions for molecule type '&
                 ,inclmol(inclnum), ' between bead ',inclbead(inclnum,1),' and bead ',inclbead(inclnum,2),', ofscale LJ: '&
                 ,ofscale(inclnum),', ofscale Q: ',ofscale2(inclnum)
             else
                write(io_output,'(3(A,I0),2(A,F6.3))') '      excluding intramolecular interactions for molecule type '&
                 ,inclmol(inclnum), ' between bead ',inclbead(inclnum,1),' and bead ',inclbead(inclnum,2),', ofscale LJ: '&
                 ,ofscale(inclnum),', ofscale Q: ',ofscale2(inclnum)
             end if
          end if
       end do

       if (lprint) then
          write(io_output,'(A,I0,A,/)') '  Total: ',inclnum,' inclusion rules for intramolecular interactions'
       end if
    end if
! -------------------------------------------------------------------
    !> Looking for section INTRAMOLECULAR_OH15
    ! read exclusion table for intermolecular interactions
    ainclnum=0
    if (myid.eq.rootid) then
       REWIND(io_input)
       CYCLE_READ_OH15:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_oh15

          if (UPPERCASE(line_in(1:19)).eq.'INTRAMOLECULAR_OH15') then
             jerr=0
             exit cycle_read_oh15
          end if
       END DO CYCLE_READ_OH15
    end if

    call mp_bcast(jerr,1,rootid,groupid)

    if (jerr.eq.0) then
       do
          if (myid.eq.rootid) then
             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section INTRAMOLECULAR_OH15',jerr)
          end if

          call mp_bcast(line_in,rootid,groupid)

          if (UPPERCASE(line_in(1:23)).eq.'END INTRAMOLECULAR_OH15') exit

          ainclnum=ainclnum+1
          ! ainclmol ainclbead_1 ainclbead_2 a15type
          read(line_in,*) ainclmol(ainclnum),ainclbead(ainclnum,1),ainclbead(ainclnum,2) ,a15t(ainclnum)

          if (lprint) then
             write(io_output,'(4(A,I0))') '      repulsive 1-5 OH interaction for molecule type ',ainclmol(ainclnum)&
              ,' between bead ',ainclbead(ainclnum,1),' and bead ',ainclbead(ainclnum,2),' of type ',a15t(ainclnum)
          end if
       end do

       if (lprint) then
          write(io_output,'(A,I0,A)') '  Total: ',ainclnum,' special rules for intramolecular 1-5 OH interactions'
       end if
    end if
!
!> \brief Read biasing potentials from bottom of file INSTEAD of on the 'nunit'
!line for each molecule type
!
!> For simulations requiring biasing potentials for a large number of molecules,
!it is a huge pain to modify them when located
!> at the end of each 'nunit' line in the fort.4. I've restored the original
!location of the biasing potentials to the end of the
!> file, where that are simply listed for each molecule in each box.
  ! initialize biasing potential to be 0 in all boxes
  eta2 = 0.0E0_dp

  ! Looking for section UNIFORM_BIASING_POTENTIALS
    if (myid.eq.rootid) then
       REWIND(io_input)
       UNIFORM_BIASING_POTENTIALS:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit UNIFORM_BIASING_POTENTIALS

          if (UPPERCASE(line_in(1:26)).eq.'UNIFORM_BIASING_POTENTIALS') then
                do imol=1,nmolty+1
                     call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                     if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section UNIFORM_BIASING_POTENTIALS',jerr)
                     if (UPPERCASE(line_in(1:30)).eq.'END UNIFORM_BIASING_POTENTIALS') then
                         if (imol.ne.nmolty+1) call err_exit(__FILE__,__LINE__,'Section UNIFORM_BIASING_POTENTIALS not complete!',jerr)
                         exit
                     else if (imol.eq.nmolty+1) then
                         call err_exit(__FILE__,__LINE__,'Section UNIFORM_BIASING_POTENTIALS has more than nmolty records!',jerr)
                     end if
                     read(line_in,*) (eta2(i,imol),i=1,nbox)

                 end do

                 exit UNIFORM_BIASING_POTENTIALS
          end if
       END DO UNIFORM_BIASING_POTENTIALS
       if (lprint) then
           write(io_output,'(/,A,/,A)') 'SECTION UNIFORM_BIASING_POTENTIALS','------------------------------------------'
           write(io_output,'(A)') 'Molecule type, biasing potential 1 through nbox [K]: '
           do imol=1,nmolty
               write(io_output,'(3(1X,F9.3))') (eta2(i,imol),i=1,nbox)
           end do
       end if
    end if
    call mp_bcast(eta2,nbxmax*ntmax,rootid,groupid)


!> \brief Read in required specification on specific atoms for atom translation
!moves.
!
!> Previously the user could only choose to do atom translations on every atom
!in the system or none at all. Now, users
!> can specify the specific atoms they want to do translations on. This still
!allows the user, if so desired, to do atom
!> translations on all atoms in the system.


  ! Looking for section SPECIFIC_ATOM_TRANSL
    if (myid.eq.rootid) then
        REWIND(io_input)
        READ_ATOM_TRANSL:DO
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section SPECIFIC_ATOM_TRANSL not found',jerr)
            if (UPPERCASE(line_in(1:20)).eq.'SPECIFIC_ATOM_TRANSL') then
                if (lprint) then
                    write(io_output,'(/,A,/,A)') 'SECTION SPECIFIC_ATOM_TRANSL','------------------------------------------'
                end if
                jerr = 0
                exit READ_ATOM_TRANSL
            end if
        end do READ_ATOM_TRANSL
    end if

    call mp_bcast(jerr,1,rootid,groupid)
    if (jerr.eq.0) then
        do
            if (myid.eq.rootid) then
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section SPECIFIC_ATOM_TRANSL',jerr)
            end if

            call mp_bcast(line_in,rootid,groupid)
            if (UPPERCASE(line_in(1:24)).eq.'END SPECIFIC_ATOM_TRANSL') exit

            ! Read in the number of atoms on which to do atom translations
            read(line_in,*) natomtrans_atoms
            allocate(atomtrans_atomlst (natomtrans_atoms))
            allocate(atomtrans_moleclst(natomtrans_atoms))

            ! Read in what those atoms are, and then what molecule they belong to

            if( natomtrans_atoms .gt. 0) then
                if (myid.eq.rootid) call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                call mp_bcast(line_in,rootid,groupid)
                read(line_in,*) (atomtrans_atomlst(j),j=1,natomtrans_atoms)
                if (myid.eq.rootid) call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                call mp_bcast(line_in,rootid,groupid)
                read(line_in,*) (atomtrans_moleclst(j),j=1,natomtrans_atoms)
            end if

            if (lprint) then
                write(io_output,*) 'natomtrans_atoms:   ', natomtrans_atoms
                write(io_output,*) 'atomtrans_atomlst:  ', (atomtrans_atomlst(j),j=1,natomtrans_atoms)
                write(io_output,*) 'atomtrans_moleclst: ', (atomtrans_moleclst(j),j=1,natomtrans_atoms)
            end if
        end do
    end if

    !> set up the inclusion table
    call inclus(inclnum,inclmol,inclbead,inclsign,ncarbon,ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)

! -------------------------------------------------------------------
    ! read kdtree related variables
    if (lkdtree) call read_kdtree(io_input)

! -------------------------------------------------------------------
    !> write out connectivity and bonded interactions
    if (lprint) then
       do imolty=1,nmolty
          write(io_output,'(/,A,I0,/,/,A)') 'Molecule type: ',imolty,'LJ INCLUSION TABLE'
          write(io_output,'(4X,'//format_n(nunit(imolty),'(1X,I3)')//')') (j,j=1,nunit(imolty))
          do j=1,nunit(imolty)
             write(io_output,'(I4,'//format_n(nunit(imolty),'(1X,L3)')//')') j,linclu(imolty,j,1:nunit(imolty))
          end do

          write(io_output,'(/,A)') 'CHARGE INCLUSION TABLE'
          write(io_output,'(4X,'//format_n(nunit(imolty),'(1X,I3)')//')') (j,j=1,nunit(imolty))
          do j=1,nunit(imolty)
             write(io_output,'(I4,'//format_n(nunit(imolty),'(1X,L3)')//')') j,lqinclu(imolty,j,1:nunit(imolty))
          end do

          write(io_output,'(/,A)') '1-4 LJ SCALING FACTORS'
          write(io_output,'(7X,'//format_n(nunit(imolty),'(1X,I6)')//')') (j,j=1,nunit(imolty))
          do j = 1, nunit(imolty)
             write(io_output,'(I7,'//format_n(nunit(imolty),'(1X,F6.3)')//')') j,ljscale(imolty,j,1:nunit(imolty))
          end do

          write(io_output,'(/,A)') '1-4 CHARGE SCALING FACTORS'
          write(io_output,'(7X,'//format_n(nunit(imolty),'(1X,I6)')//')') (j,j=1,nunit(imolty))
          do j = 1, nunit(imolty)
             write(io_output,'(I7,'//format_n(nunit(imolty),'(1X,F6.3)')//')') j,qscale2(imolty,j,1:nunit(imolty))
          end do

          if (nunit(imol).gt.1) then
             write(io_output,'(/,A)') '      i      j type_i type_j   bond length           k/2'
          end if
          do i=1,nunit(imol)
             do j=1,invib(imol,i)
                write(io_output,'(4(1X,I6),1X,F13.4,F14.1)') i,ijvib(imol,i,j),atoms%list(ntype(imol,i))&
                 ,atoms%list(ntype(imol,ijvib(imol,i,j))),brvib(itvib(imol,i,j)),brvibk(itvib(imol,i,j))
             end do
          end do

          if (nunit(imol).gt.2) then
             write(io_output,'(/,A)') '      i      j      k type_i type_j type_k        angle        k/2'
          end if
          do i=1,nunit(imol)
             do j=1,inben(imol,i)
                write(io_output,'(6(1X,I6),2(1X,F12.2))') i,ijben2(imol,i,j),ijben3(imol,i,j),atoms%list(ntype(imol,i))&
                 ,atoms%list(ntype(imol,ijben2(imol,i,j))),atoms%list(ntype(imol,ijben3(imol,i,j)))&
                 ,brben(itben(imol,i,j))*raddeg,brbenk(itben(imol,i,j))
             end do
          end do

          if (nunit(imol).gt.3) then
             write(io_output,'(/,A)') '      i      j      k      l type_i type_j type_k type_l torsion type'
          end if
          do i=1,nunit(imol)
             do j=1,intor(imol,i)
                write(io_output,'(9(1X,I6))') i,ijtor2(imol,i,j),ijtor3(imol,i,j),ijtor4(imol,i,j)&
                 ,atoms%list(ntype(imol,i)),atoms%list(ntype(imol,ijtor2(imol,i,j))),atoms%list(ntype(imol,ijtor3(imol,i,j)))&
                 ,atoms%list(ntype(imol,ijtor4(imol,i,j))),dihedrals%list(ittor(imol,i,j))
             end do
          end do
       end do
    end if
! -------------------------------------------------------------------
    ! write out non-bonded interaction table
    if (lprint) then
       write(io_output,'(/,A)') 'PAIRWISE LJ AND COULOMB INTERACTIONS'
       write(io_output,'(A)') '    i    j         q0(i)         q0(j)     vvdW_1     vvdW_2 ...'
       do i=1,nntype
          do j=1,nntype
             if (lhere(i).and.lhere(j)) then
                ij=type_2body(i,j)
                write(io_output,'(2(1X,I4),2(1X,F13.6),'//format_n(vdW_nParameter(nonbond_type(ij)),'(1X,G12.5)')//')')&
                 atoms%list(i),atoms%list(j),qelect(i),qelect(j),vvdW(1:vdW_nParameter(nonbond_type(ij)),ij)
             end if
          end do
       end do
    end if
! -------------------------------------------------------------------
    pmswat=0.0_dp
    pmflcq=0.0_dp
    pmexpc=0.0_dp
    pmexpc1=0.0_dp
    pm_atom_tra=0.0_dp
    if (lgibbs) then
       pmvol=2.0_dp/nchain
       pmswap=(1.0_dp-pmvol)/4.0_dp+pmvol
       pmcb=(1.0_dp-pmswap)/3.0_dp+pmswap
    else
       pmswap=0.0_dp
       if (lnpt) then
          pmvol=2.0_dp/nchain
       else
          pmvol=0.0_dp
       end if
       pmcb=(1.0_dp-pmvol)/3.0_dp+pmvol
    end if
    pmtra=(1.0_dp-pmcb)/2.0_dp+pmcb

    call read_transfer(io_input,lprint)
    call init_moves_volume(io_input,lprint)
    call init_swatch(io_input,lprint)
    call init_swap(io_input,lprint)
    call init_cbmc(io_input,lprint)
! -------------------------------------------------------------------

    ! Colin Bunner- for osmotic NVT-GEMC simulations, check that at least one molecule type isn't being swapped.
    ! This code check needs to be here because init_swap() reads the swap section, which then populates pmswap.
    if (losmoticnvt) then
       tmp_logical=.true.
       ! Check outright for a zero swap probability
       do imolty = 1, nmolty-1
          if (pmswmt(imolty).lt.1E-6) then
              tmp_logical=.false.
          end if
       end do
       ! Check for a zero swap probability by taking the difference of pmswmt
       ! values, which are cumulative. (abs() necessary for floating point
       ! precision and because we always pad with an additional molecule type)
       do imolty = 1, nmolty-1
          if (abs(pmswmt(imolty+1)-pmswmt(imolty)).lt.1E-6) then
             tmp_logical=.false.
          end if
       end do
       ! If tmp_logical is still true, all molecules are being swapped and the simulation is incorrect.
       if (tmp_logical) call err_exit(__FILE__,__LINE__,'Cannot swap all molecule types in an osmotic NVT-GEMC simulation. &
                                                        & If you are trying to do an osmotic NVT-GEMC simulation, make sure &
                                                        & that swap moves are not specified for the component that your membrane is &
                                                        & impermeable to (look at pmswmt in your fort.4 file and remember these &
                                                        & probabilities are cumulative). If you want a regular NVT-GEMC simulation, &
                                                        & or any other simulation type, make sure losmoticnvt=F in namelist system in &
                                                        & your topmon.inp (or just delete this line because F is the default.',myid+1)

    end if

    !> read namelist mc_flucq
    fqtemp=5.0_dp
    rmflucq=0.1_dp
    lflucq=.false.
    lqtrans=.false.
    fqegp=0.0_dp
    nchoiq=1
    do i=1,nmolty
       pmfqmt(i)=real(i,dp)/nmolty
    end do

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=mc_flucq,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_flucq',jerr)
    end if

    call mp_bcast(taflcq,1,rootid,groupid)
    call mp_bcast(fqtemp,1,rootid,groupid)
    call mp_bcast(rmflucq,1,rootid,groupid)
    call mp_bcast(pmflcq,1,rootid,groupid)
    call mp_bcast(pmfqmt,ntmax,rootid,groupid)
    call mp_bcast(lflucq,ntmax,rootid,groupid)
    call mp_bcast(lqtrans,ntmax,rootid,groupid)
    call mp_bcast(fqegp,ntmax,rootid,groupid)
    call mp_bcast(nchoiq,nbxmax,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_FLUCQ','------------------------------------------'
       write(io_output,'(A,F4.2)') 'target fluctuating charge acceptance ratio (taflcq): ',taflcq
       write(io_output,'(A,F8.3,A)') 'fluctuating charge temperature: ',fqtemp,' [K]'
       write(io_output,'(A,G16.9)') 'initial maximum displacement for fluctuating charge moves: ',rmflucq
       write(io_output,'(A,G16.9)') 'pmflcq: ',pmflcq
       write(io_output,'(A,'//format_n(nbox,'(2X,I0)')//')') '   nchoiq for each box: ',nchoiq(1:nbox)
       write(io_output,'(A,I0)') 'nswapq: ',nswapq
       write(io_output,'(/,A)') 'molecule type:  lflucq lqtrans   pmfqmt            fqegp'
    end if

    if (fqtemp .lt. 0.0_dp) then
       fqbeta = 0.0_dp
       nswapq = 1
    else
       fqbeta = 1.0_dp/fqtemp
    endif
    rmflcq = rmflucq

    !> allow the MC_FLUCQ section to be skipped for non-polarizable models
    needMFsection = .false.

    do i=1,nmolty
       if (lflucq(i)) then
          needMFsection = .true.
       endif

       if (lprint) then
          write(io_output,'(I13,A,2(1X,L7),1X,F8.4,1X,F16.4)') i,':',lflucq(i),lqtrans(i),pmfqmt(i),fqegp(i)
       end if

       if (lflucq(i).and..not.lelect(i)) call err_exit(__FILE__,__LINE__,'lelect must be true if flucq is true',myid+1)
       if (lqtrans(i).and..not.lflucq(i)) call err_exit(__FILE__,__LINE__,'lflucq must be true if intermolecular CT is allowed'&
        ,myid+1)
    end do

    if (ALL(.NOT.lflucq(1:nmolty))) then
       if (lanes) call err_exit(__FILE__,__LINE__,'lanes should be false for nonpolarizable systems!',myid+1)
       if (lfepsi) call err_exit(__FILE__,__LINE__,'lfepsi should be false for nonpolarizable systems!',myid+1)
    end if

    !> Additional MC_FLUCQ section to read in jayq, jayself
    if (needMFsection) then
       REWIND(io_input)
       CYCLE_READ_FQ:DO
          call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section MC_FLUCQ not found',jerr)
          if (UPPERCASE(line_in(1:12)).eq.'END MC_FLUCQ') exit cycle_read_fq
          if (UPPERCASE(line_in(1:8)).eq.'MC_FLUCQ') then

             do i=1,nntype
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_FLUCQ',jerr)

                !> fqmolty = molecule type number for the fluctuating charge molecule
                !> fqtyp = ntii = line number where this bead appears in the ATOMS section
                !> xiq(fqtyp) = hardness parameter for each bead
                !> jayself(fqtype) = coulomb paramter for each bead
                read(line_in,*) fqmolty, fqtyp, xiq(fqtyp), jayself(fqtyp)
                !write(6,*) fqmolty, fqtyp, xiq(fqtyp), jayself(fqtyp)
             end do

             !write(6,*) nntype
             do i=1,nntype*nntype
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_FLUCQ',jerr)

                !> fqcrosstyp = ntij = (ntii - 1)*nntype + ntjj (nntype is total number of lines in ATOMS section)
                !> jayq(fqcrosstyp) = coulomb integral parameter for interacting beads
                read(line_in,*) fqmolty, fqcrosstyp, jayq(fqcrosstyp)
                !write(6,*) fqmolty, fqcrosstyp, jayq(fqcrosstyp)
             end do
          end if
       end do CYCLE_READ_FQ

       if (lprint) then
          do i=1,nmolty
             if (lflucq(i)) then
                write(io_output,'(/,A,I0)') 'FQ parameters for molecule type ',i
                do ii=1,nunit(i)
                   ntii = ntype(i,ii)
                   write(io_output,'(/,A,I0,A,2(1X,F14.4))') "fq self terms, bead ",ii,":", xiq(ntii), jayself(ntii)
                   do jj=1,nunit(i)
                      ntjj = ntype(i,jj)
                      ntij = type_2body(ntii,ntjj)
                      if (ii .ne. jj) write(io_output,'(A,I0,A,I0,A,I0,A,3(1X,F14.4))') "fq cross term btwn ",ii,",",jj," (",ntij,") :", jayq(ntij)
                   end do
                end do
             end if
          end do
       end if
    endif

! -------------------------------------------------------------------
    !> read information for grand-canonical ensemble simulations
    if (lgrand) then
       B=0.0_dp

       if (myid.eq.rootid) then
          rewind(io_input)
          read(UNIT=io_input,NML=gcmc,iostat=jerr)
          if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: gcmc',jerr)
       end if

       call mp_bcast(B,ntmax,rootid,groupid)
       call mp_bcast(nequil,1,rootid,groupid)
       call mp_bcast(ninstf,1,rootid,groupid)
       call mp_bcast(ninsth,1,rootid,groupid)
       call mp_bcast(ndumph,1,rootid,groupid)

       if (lprint) then
          write(io_output,'(/,A,/,A)') 'NAMELIST GCMC','------------------------------------------'
          ! information for histogram output (added 8/30/99)
          w_i(nequil)
          w_i(ninstf)
          w_i(ninsth)
          w_i(ndumph)

          ! chemical potentials (added 8/30/99 by jpotoff)
          do i = 1,nmolty
             write(io_output,'(A,I0,A,G16.9)') 'chemical potential for molecule type ',i,': ',B(i)
          end do
       end if

       ! convert chemical potentials to activities
       ! This B(i) goes in the acceptance rules
       do i=1,nmolty
          debroglie = debroglie_factor*sqrt(beta/masst(i))
          B(i) = exp(B(i)/temp)/(debroglie**3)
       end do

       boxlx(2:)=boxlx(1)
       boxly(2:)=boxly(1)
       boxlz(2:)=boxlz(1)
       lideal(2:)=.true.
    end if
! -------------------------------------------------------------------
    call init_ee(io_input,lprint)
    call init_moves_simple(io_input,lprint)

    ! write out move probabilities, in percentage
    if (lprint) then
       write(io_output,'(/,A)') 'percentage move probabilities:'

       pcumu=0.0_dp
       if (pmvol.gt.pcumu) then
           pm = min(1.0_dp,pmvol)
           pcumu = pmvol
       else
           pm=0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' volume move       :',100.0_dp*pm,' %'

       if (pmswat.gt.pcumu) then
          pm = min(1.0_dp,pmswat - pcumu)
          pcumu = pmswat
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' swatch move       :',100.0_dp*pm,' %'

       if (pmswap.gt.pcumu) then
          pm = min(1.0_dp,pmswap - pcumu)
          pcumu = pmswap
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' swap move         :',100.0_dp*pm,' %'

       if (pmcb.gt.pcumu) then
          pm = min(1.0_dp,pmcb - pcumu)
          pcumu = pmcb
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' CBMC move         :',100.0_dp*pm,' %'

       if (pmflcq.gt.pcumu) then
          pm = min(1.0_dp,pmflcq - pcumu)
          pcumu = pmflcq
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' fluct. charge move:',100.0_dp*pm,' %'

       if (pmexpc.gt.pcumu) then
          pm = min(1.0_dp,pmexpc - pcumu)
          pcumu = pmexpc
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' expanded ens. move:',100.0_dp*pm,' %'

       if (pmexpc1.gt.pcumu) then
          pm = min(1.0_dp,pmexpc1 - pcumu)
          pcumu = pmexpc1
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' new EE move       :',100.0_dp*pm,' %'

       if (pm_atom_tra.gt.pcumu) then
          pm = min(1.0_dp,pm_atom_tra - pcumu)
          pcumu = pm_atom_tra
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' atom trans. move  :',100.0_dp*pm,' %'

       if (pmtra.gt.pcumu) then
          pm = min(1.0_dp,pmtra - pcumu)
          pcumu = pmtra
       else
          pm = 0.0_dp
       end if
       write(io_output,'(A,F8.2,A)') ' translation move  :',100.0_dp*pm,' %'

       pm = max(0.0_dp,1.0_dp - pcumu)
       write(io_output,'(A,F8.2,A)') ' rotation move     :',100.0_dp*pm,' %'
    end if
! -------------------------------------------------------------------
    !> set up force field parameters and read in external potentials
    call init_ff(io_input,lprint)
    if (myid.eq.rootid) close(io_input)
! ===================================================================
    !> Initialize the system or read configuration from the restart file

    ! Q. Paul C. -- for tabulated CBMC bending growth
    if (L_cbmc_bend) then
        call setup_cbmc_bend(file_cbmc_bend)
    end if

    if (linit) then
       do ibox = 1,nbox
          if (lsolid(ibox).and..not.lrect(ibox).and..not.(ibox.eq.1.and.lexzeo)) call err_exit(__FILE__,__LINE__&
           ,'Cannot initialize non-rectangular system',myid+1)
       end do
       call setup_system_config(file_struct)
       nnstep = 0
    else if (myid.eq.rootid) then
       ! begin read restart file with
       ! rootid----------------------------------------------------------------------------
       if (use_checkpoint) file_restart='save-config'

       io_restart=get_iounit()
       ! MJM recl needed for long input records, like nboxi and moltyp for a lot of molecules
       open(unit=io_restart,access='sequential',action='read',file=file_restart,form='formatted',iostat=jerr,recl=4096,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open main restart file '//trim(file_restart),myid+1)

       read(io_restart,*) nnstep
       read(io_restart,*) Armtrax,Armtray,Armtraz
       do im=1,nbox
          do imol = 1,nmolty
             read(io_restart,*) rmtrax(imol,im),rmtray(imol,im),rmtraz(imol,im)
             read(io_restart,*) rmrotx(imol,im),rmroty(imol,im),rmrotz(imol,im)
          end do
          call averageMaximumDisplacement(im,imol)
       end do

       do im=1,nbox
          read(io_restart,*) (rmflcq(i,im),i=1,nmolty)
       end do

       ! changed so fort.77 the same for all ensembles
       ! 06/08/09 KM

       read(io_restart,*) (rmvol(ibox),ibox=1,nbox)

       do ibox=1,nbox
          if (lsolid(ibox).and..not.lrect(ibox)) then
             read(io_restart,*) (rmhmat(ibox,j),j=1,9)
          end if
       end do

       if (lprint) then
          write(io_output,'(/,A,/,A)') 'READING CONFIGURATION FROM RESTART FILE','------------------------------------------'
          write(io_output,'(A)') 'new maximum displacements read from restart-file'
          write(io_output,'(A,3(2X,F10.6))') '   max atom trans. displacement: ',Armtrax,Armtray,Armtraz
          write(io_output,'(A,3(1X,E11.4))') '   max volume displacement: ',(rmvol(ibox), ibox = 1,nbox)
          do ibox=1,nbox
             if (lsolid(ibox).and..not.lrect(ibox)) then
                write(io_output,'(A,I0)') '   rmhmat for box ',ibox
                write(io_output,'(3(4X,G16.9))') rmhmat(ibox,1),rmhmat(ibox,4),rmhmat(ibox,7)
                write(io_output,'(3(4X,G16.9))') rmhmat(ibox,2),rmhmat(ibox,5),rmhmat(ibox,8)
                write(io_output,'(3(4X,G16.9))') rmhmat(ibox,3),rmhmat(ibox,6),rmhmat(ibox,9)
             end if
          end do
          do im=1,nbox
             write(io_output,'(/,A,I0)') 'box      #',im
             do imol = 1,nmolty
                write(io_output,'(A,I0)') '   molecule type ',imol
                write(io_output,'(A,3(1X,F10.6))') '      max trans. displacement:  ',rmtrax(imol,im),rmtray(imol,im)&
                 ,rmtraz(imol,im)
                write(io_output,'(A,3(1X,F10.6))') '      max rot. displacement:    ',rmrotx(imol,im),rmroty(imol,im)&
                 ,rmrotz(imol,im)
                write(io_output,'(A,3(1X,F10.6))') '      max fluc. q displacement: ',rmflcq(imol,im)
             end do
          end do
          write(io_output,'(/,A)') 'reading new box size from restart-file'
       end if

       ! read box dimension info
       do ibox=1,nbox
          if (ibox.eq.1.and.lexzeo) then
             if (lsolid(ibox).and..not.lrect(ibox)) then
                read(io_restart,*) dum,dum,dum,dum,dum,dum,dum,dum,dum
             else
                read(io_restart,*) dum,dum,dum
             end if
          else if (lsolid(ibox) .and. .not. lrect(ibox)) then
             read(io_restart,*) (hmat(ibox,j),j=1,9)
             if (lprint) then
                write(io_output,'(A,I0)') 'ibox:  ', ibox
                write(io_output,'(A)') '   HMAT COORDINATES'
                write(io_output,'(3(4X,G16.9))') hmat(ibox,1),hmat(ibox,4),hmat(ibox,7)
                write(io_output,'(3(4X,G16.9))') hmat(ibox,2),hmat(ibox,5),hmat(ibox,8)
                write(io_output,'(3(4X,G16.9))') hmat(ibox,3),hmat(ibox,6),hmat(ibox,9)
             end if

             call matops(ibox)

             if (lprint) then
                write(io_output,'(A,3(1X,F11.6))') 'min widths:',min_width(ibox,1:3)
                write(io_output,'(A,2X,F12.3)') "cell length |a|:",cell_length(ibox,1)
                write(io_output,'(A,2X,F12.3)') "cell length |b|:",cell_length(ibox,2)
                write(io_output,'(A,2X,F12.3)') "cell length |c|:",cell_length(ibox,3)
                write(io_output,'(A,2X,F12.3)') "cell angle alpha:",cell_ang(ibox,1)*raddeg
                write(io_output,'(A,2X,F12.3)') "cell angle beta: ",cell_ang(ibox,2)*raddeg
                write(io_output,'(A,2X,F12.3)') "cell angle gamma:",cell_ang(ibox,3)*raddeg
             end if

             if ((allow_cutoff_failure.lt.0).and.ANY(rcut(ibox)/min_width(ibox,1:3).gt.0.5_dp)) then
                call err_exit(__FILE__,__LINE__,'rcut > half cell width',myid+1)
             end if
          else
             ! "normal" box info
             read(io_restart,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
             if (lprint) then
                write(io_output,"(' dimension box ',I0,': a = ',F12.6,'  b = ',F12.6,'  c = ',F12.6,' rcut =',F12.6)") ibox,boxlx(ibox)&
                 ,boxly(ibox),boxlz(ibox),rcut(ibox)
             end if
             ! check if rcut too large and address accordingly
             if(((rcut(ibox)/boxlx(ibox).gt.0.5_dp).or.(rcut(ibox)/boxly(ibox).gt.0.5_dp)&
                  .or.(rcut(ibox)/boxlz(ibox).gt.0.5_dp)).and.lcutcm) then
                if (allow_cutoff_failure.lt.0) then
                   ! error exit
                   call err_exit(__FILE__,__LINE__,'rcut > 0.5*boxlx',myid+1)
                else if (allow_cutoff_failure.eq.1) then
                   ! fix rcut to half minimum boxlength
                   rcut(ibox) = min(boxlx(ibox),boxly(ibox))/2.0_dp
                   if (lpbcz) then
                      rcut(ibox)=min(rcut(ibox),boxlz(ibox)/2.0_dp)
                   end if
                   if (rcut(ibox).le.0) then
                      call err_exit(__FILE__,__LINE__,'rcut initialization&
                                                        & error. Maybe a boxlength is negative', myid+1)
                   end if
                   ! restore maximum displacement for translation if needed
                   call restore_displ_transl(ibox)
                end if
             else if (allow_cutoff_failure.eq.1) then
                   call restore_displ_transl(ibox) ! check if need to restore max displ. trans.
                                                   ! and perform if needed
             end if
          end if
       end do

       if (lprint) then
          write(io_output,'(/,A)') 'Finished writing simulation box related info'
       end if

       read(io_restart,*) ncres
       read(io_restart,*) nmtres
       if (lgrand) nchain=ncres

       ! check if the number of particles in fort.4 & fort.77 agree
       if (ncres.ne.nchain.or.nmtres.ne.nmolty) then
          write(io_output,*) 'nchain',nchain,'ncres',ncres
          write(io_output,*) 'nmolty',nmolty,'nmtres',nmtres
          call err_exit(__FILE__,__LINE__,'conflicting information in restart and control files',myid+1)
       end if

       read(io_restart,*) nures(1:nmolty)

       do i=1,nmolty
          if (nures(i).ne.nunit(i)) then
             write(io_output,*) 'unit',i,'nunit',nunit(i),'nures',nures(i)
             call err_exit(__FILE__,__LINE__,'conflicting information in restart and control files',myid+1)
          end if
       end do

       read(io_restart,*) moltyp(1:nchain)
       read(io_restart,*) nboxi(1:nchain)
       if (any(lexpand(1:nmolty))) then
          do i=1,nmolty
             if (lexpand(i)) read(io_restart,*) eetype(i)
          end do
          do i=1,nmolty
             if (lexpand(i)) read(io_restart,*) rmexpc(i)
          end do
       end if

       qbox = 0.0_dp
       do i=1,nchain
          imolty = moltyp(i)
          do j=1,nunit(imolty)
             read(io_restart,*) rxu(i,j),ryu(i,j),rzu(i,j),qqu(i,j)
             if (.not.lreadq) then
                qqu(i,j) = qelect(ntype(imolty,j))
             end if
             qbox(nboxi(i)) = qbox(nboxi(i)) + qqu(i,j)
          end do
       end do

       close(io_restart)

       call check_rigid_structures()


       do i=1,nbox
          if (i.eq.1.and.(lexzeo.or.lionic)) then
             cycle
          else if (abs(qbox(i)).gt.1E-6_dp .and. i.ne.gcbmc_box_num) then
             call err_exit(__FILE__,__LINE__,'box '//integer_to_string(i)//' has a net charge of '//real_to_string(qbox(i)),myid+1)
          end if
       end do
       ! end read restart file with
       ! rootid----------------------------------------------------------------------------
    end if

    ! broadcast variables read in to other processors
    call mp_bcast(nnstep,1,rootid,groupid)
    call mp_bcast(Armtrax,1,rootid,groupid)
    call mp_bcast(Armtray,1,rootid,groupid)
    call mp_bcast(Armtraz,1,rootid,groupid)
    call mp_bcast(rmtrax(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmtray(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmtraz(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmrotx(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmroty(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmrotz(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmflcq(1:nmolty,1:nbox),nmolty*nbox,rootid,groupid)
    call mp_bcast(rmvol(1:nbox),nbox,rootid,groupid)
    call mp_bcast(boxlx(1:nbox),nbox,rootid,groupid)
    call mp_bcast(boxly(1:nbox),nbox,rootid,groupid)
    call mp_bcast(boxlz(1:nbox),nbox,rootid,groupid)
    if (allow_cutoff_failure.eq.1) call mp_bcast(rcut(1:nbox),nbox,rootid,groupid)
    do ibox=1,nbox
       if (lsolid(ibox).and..not.lrect(ibox)) then
          call mp_bcast(rmhmat(ibox,1:9),9,rootid,groupid)
          call mp_bcast(hmat(ibox,1:9),9,rootid,groupid)
          call matops(ibox)
       end if
    end do
    call mp_bcast(nchain,1,rootid,groupid)
    call mp_bcast(moltyp,nchain,rootid,groupid)
    call mp_bcast(nboxi,nchain,rootid,groupid)
    if (any(lexpand(1:nmolty))) then
       call mp_bcast(eetype,nmolty,rootid,groupid)
       call mp_bcast(rmexpc,nmolty,rootid,groupid)
    end if
    call mp_bcast(rxu(1:nchain,1:numax),nchain*numax,rootid,groupid)
    call mp_bcast(ryu(1:nchain,1:numax),nchain*numax,rootid,groupid)
    call mp_bcast(rzu(1:nchain,1:numax),nchain*numax,rootid,groupid)
    call mp_bcast(qqu(1:nchain,1:numax),nchain*numax,rootid,groupid)
! ===================================================================
    !> check that particles are in correct boxes
    !> obtain nchbox, ncmt, parbox, parall
    nchbox = 0
    ncmt = 0
    parbox = 0
    parall = 0
    idummy = 0
    do i=1,nmolty
       if (lexpand(i)) then
          !> \bug problem in expand ensemble
          do j=1,numcoeff(i)
             do ibox=1,2
                ncmt2(ibox,i,j) = 0
             end do
          end do
       end if
    end do

    do i=1,nchain
       ibox = nboxi(i)
       if(ibox.le.nbox) then
          if (ibox.ne.1.and..not.lgibbs) then
             call err_exit(__FILE__,__LINE__,'Particle found outside BOX 1',myid+1)
          end if

          nchbox(ibox) = nchbox(ibox) + 1
          imolty = moltyp(i)
          ncmt(ibox,imolty) = ncmt(ibox,imolty) + 1
          parbox(ncmt(ibox,imolty),ibox,imolty) = i
          idummy(imolty) = idummy(imolty) + 1
          parall(imolty,idummy(imolty)) = i

          if (lexpand(imolty)) then
             if (ibox.gt.2) then
                call err_exit(__FILE__,__LINE__,'put in box 1 and 2 for ee molecules',myid+1)
             end if
             itype = eetype(imolty)
             ncmt2(ibox,imolty,itype) = ncmt2(ibox,imolty,itype) + 1
             do j=1,nunit(imolty)
                sigma_f(imolty,j) = sigm(imolty,j,itype)
                epsilon_f(imolty,j) = epsil(imolty,j,itype)
             end do
          end if
       else
          write(io_output,*) 'i:',i,'nboxi(i)',ibox
          call err_exit(__FILE__,__LINE__,'Particle found in ill-defined box',myid+1)
       end if
    end do

    do ibox=1,nbox
       if (sum(ncmt(ibox,1:nmolty)).ne.nchbox(ibox)) then
          write(io_output,*) 'box ',ibox,', nchbox: ',nchbox(ibox),', ncmt: ',ncmt(ibox,1:nmolty)
          call err_exit(__FILE__,__LINE__,'readdat: nchbox(ibox) .ne. sum(ncmt(ibox,1:nmolty))',myid+1)
       end if
    end do

    ! check that number of particles of each type is consistent
    do i=1,nmolty
       if (lgrand) temtyp(i)=sum(ncmt(1:nbox,i))

       if (sum(ncmt(1:nbox,i)).ne.temtyp(i)) then
          write(io_output,*) 'type ',i,', temtyp: ',temtyp(i),', ncmt: ',ncmt(1:nbox,i)
          call err_exit(__FILE__,__LINE__,'Particle type number inconsistency',myid+1)
       end if
    end do

    if (.not.linit) call get_molecule_config(file_struct,linit=.false.)


! -------------------------------------------------------------------
    !> Set up Ewald parameters, write out both intra- and inter-molecular interaction parameters
    if (lewald) then
       if (lchgall) then
          ! if real space term are summed over all possible pairs in the box
          ! kalp(1) & kalp(2) are fixed while calp(1) & calp(2) change according
          ! to the boxlength, kalp(1) should have a value greater than 5.0
          if (L_Ewald_Auto) kalp = 6.4_dp
          do ibox=1,nbox
             if (lsolid(ibox).and.(.not.lrect(ibox))) then
                min_boxl=minval(min_width(ibox,1:3))
             else
                min_boxl=min(boxlx(ibox),boxly(ibox),boxlz(ibox))
             end if
             calp(ibox) = kalp(ibox)/min_boxl
             if (kalp(ibox).lt.5.6E0_dp.and.lprint) then
                write(io_output,'(A,I0,2(A,G16.9))') 'Warning, kalp too small in box ',ibox,': calp = ',calp(ibox)&
                 ,', min_boxl = ',min_boxl
             end if
          end do
       else
          ! if not lchgall, calp(1) & calp(2) are fixed
          do ibox=1,nbox
             if (L_Ewald_Auto) kalp(ibox) = 3.2_dp/rcut(ibox)
             calp(ibox) = kalp(ibox)
             if (calp(ibox)*rcut(ibox).lt.2.8E0_dp.and.lprint) then
                ! JLR 11-24-09
                ! you may want a smaller kalp, e.g. when comparing to previous work
                ! This does not need to be an error
                write(io_output,'(A,I0,2(A,G16.9))') 'Warning, kalp too small in box ',ibox,': calp = ',calp(ibox)&
                 ,', rcut = ',rcut(ibox)
             end if
          end do
       end if

       if (lprint) then
          write(io_output,'(/,A)') '****Ewald Parameters*****'
          write(io_output,'(A)') 'ibox:      calp  kmaxl  kmaxm  kmaxn         rcut'
          do ibox=1,nbox
             call compute_kmax(ibox)
             write(io_output,'(I4,A,F9.3,3(1X,I6),1X,F12.4)') ibox,': ',calp(ibox),k_max_l(ibox),k_max_m(ibox),k_max_n(ibox)&
              ,rcut(ibox)
          end do
       end if
    else if (lchgall) then
       !kea
       call err_exit(__FILE__,__LINE__,'lewald should be true when lchgall is true',myid+1)
    end if
! -------------------------------------------------------------------
    ! zeolite external potential
    if (lexzeo) call suzeo(lprint)
    call readThreeBody(file_in)
    call readFourBody(file_in)
! ===================================================================
    !> write file headers
    if (myid.eq.rootid) then
       ! write out movie-header
       !*** The very first frame (nnstep=0) is no longer written out because it's usually useless
       if (imv.le.nstep) then
          io_movie=get_iounit()
          open(unit=io_movie,access='stream',action='write',file=file_movie,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open movie file '//trim(file_movie),jerr)
          if (lprint) then
             nhere = 0
             do izz=1,nntype
                if ( lhere(izz) ) then
                   nhere = nhere + 1
                   temphe(nhere) = izz
                end if
             end do

             if (ANY(pmgroup(1:nmolty).gt.0) .and. .not. l_gcbmc_movie) then
                 ! The information about group-CBMC reservoir molecules will not be written in the movie file
                 write(io_movie,*) nstep/imv,nchain-nchbox(gcbmc_box_num),nmolty,nbox-1,nhere
                 write(io_movie,*) (rcut(ibox),ibox=1,nbox-1) !< assume the last box is the reservoir box!!
             else
                 write(io_movie,*) nstep/imv,nchain,nmolty,nbox,nhere
                 write(io_movie,*) (rcut(ibox),ibox=1,nbox)
             end if

             write(io_movie,*) (temphe(izz),izz=1,nhere)

             do imolty = 1,nmolty
                write(io_movie,*) nunit(imolty)
                ! output bond connectivity information
                do ii=1,nunit(imolty)
                   write(io_movie,*) invib(imolty,ii),(ijvib(imolty,ii,z),z=1,invib(imolty,ii))
                end do

                ! output torsional connectivity information
                do j = 1,nunit(imolty)
                   write(io_movie,*) intor(imolty,j),(ijtor2(imolty,j,ii),ijtor3(imolty,j,ii),ijtor4(imolty,j,ii)&
                    ,ii=1,intor(imolty,j))
                end do
             end do
          end if
       else
          io_movie=-1
       end if

       ! open xyz movie file for individual simulation box
       if (L_movie_xyz) then
          do ibox = 1,nbox
             string = integer_to_string(ibox)
             file_box_movie = "box"//string(1:len_trim(string))//"movie"//run_num(1:len_trim(run_num))//suffix//".xyz"
             !write(file_box_movie,'("box",I1.1,"movie",I1.1,A,".xyz")') ibox ,run_num,suffix
             io_box_movie(ibox)=get_iounit()
             open(unit=io_box_movie(ibox),access='stream',action='write',file=file_box_movie,form='formatted',iostat=jerr&
              ,status='unknown')
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open box movie file '//trim(file_box_movie),jerr)
          end do
       else
          io_box_movie=-1
       end if

       ! open pdb movie file for individual simulation box
       if (L_movie_pdb) then
          do ibox = 1,nbox
             string = integer_to_string(ibox)
             file_box_movie_pdb = "box"//string(1:len_trim(string))//"movie"//run_num(1:len_trim(run_num))//suffix//".pdb"
             !write(file_box_movie,'("box",I1.1,"movie",I1.1,A,".xyz")') ibox
             !,run_num,suffix
             io_box_movie_pdb(ibox)=get_iounit()
             open(unit=io_box_movie_pdb(ibox),access='stream',action='write',file=file_box_movie_pdb,form='formatted'&
              ,iostat=jerr,status='unknown')
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open box movie file '//trim(file_box_movie_pdb),jerr)
          end do
       else
          io_box_movie_pdb=-1
       end if

       ! write out isolute movie header
       if (ANY(isolute(1:nmolty).le.nstep)) then
          io_solute=get_iounit()
          open(unit=io_solute,access='stream',action='write',file=file_solute,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open solute movie file '//trim(file_solute),jerr)
          if (lprint) then
             write(io_solute,*) nmolty
             do imol = 1, nmolty
                write(io_solute,*) imol,nunit(imol),(nstep/isolute(imol))*temtyp(imol)
             end do
          end if
       else
          io_solute=-1
       end if

       ! cell parameters using crystallographic convention
       if (ANY(lsolid.and..not.lrect)) then
          file_cell = "cell_param"//run_num(1:len_trim(run_num))//suffix//".dat"
          !write(file_cell,'("cell_param",I1.1,A,".dat")') run_num,suffix
          io_cell=get_iounit()
          open(unit=io_cell,access='stream',action='write',file=file_cell,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open cell file '//trim(file_cell),jerr)
       else
          io_cell=-1
       end if

       ! set up info at beginning of fort.12 for analysis
       if (ltraj) then
          io_traj=get_iounit()
          open(unit=io_traj,access='stream',action='write',file=file_traj,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open traj file '//trim(file_traj),jerr)
          if (lprint) write(io_traj,*) nstep,iratp,nbox,nmolty,(masst(i),i=1,nmolty)
       end if
    end if
! -------------------------------------------------------------------
    ! KM 01/10 remove analysis
    ! if (ianalyze.lt.nstep) then
    !    nhere = 0
    !    do izz=1,nntype
    !       if ( lhere(izz) ) then
    !          nhere = nhere + 1
    !          temphe(nhere) = izz
    !          beadtyp(nhere)=izz
    !       end if
    !    end do
    !    do izz = 1,nhere
    !       atemp = temphe(izz)
    !       decode(atemp) = izz
    !    end do
    ! end if
! ===================================================================

    deallocate(ncarbon,idummy,qbox,nures,lhere,temphe,inclmol,inclbead,inclsign,ofscale,ofscale2&
     ,ainclmol,ainclbead,a15t,stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readdat: deallocation failed',jerr)

    return
  end subroutine readdat

!> \brief Perform periodic updates and block averages
  subroutine monper()
    use const_math,only:raddeg
    use const_phys,only:N_Avogadro,R_gas,MPa2SimUnits
    use util_math,only:update_average,store_block_average
    use util_string,only:format_n
    use util_kdtree,only:update_tree_height,construct_kdtree, check_tree_coord
    use energy_intramolecular,only:U_bonded
    use energy_pairwise,only:energy,coru
    use moves_simple,only:update_translation_rotation_max_displacement
    use moves_volume,only:update_volume_max_displacement
    use moves_cbmc,only:opt_safecbmc
    use transfer_shared,only:opt_bias,lopt_bias,freq_opt_bias
    use transfer_swap,only:acchem,bnchem
    use prop_pressure,only:pressure
    use parser_pdb,only:writePDBmovie
    use prop_widom,only:blk_avg_prop_widom

    logical::lfq,ovrlap
    integer::im,i,ibox,jbox,Temp_nmol,m,mm,imolty,nummol,ii,ntii,intg,ilunit,k,j,itype,itel,zzz,jmolty
    real::ratflcq,press1,surf,comp,temvol,Temp_Mol_Vol,Temp_Energy,Heat_vapor_T,Heat_vapor_LJ,Heat_vapor_COUL,CED_T,CED_LJ&
     ,CED_COUL,HSP_T,HSP_LJ,HSP_COUL,pdV,setx,sety,setz,setel,v(nEnergy),vol,rho,temmass,dpr,dpp

#ifdef __DEBUG__
    write(io_output,*) 'begin MONPER in ',myid
#endif
! -------------------------------------------------------------------

    ! re-construct mol_tree if the height exceeds two times the optimal height
    ! NOT used when volume move is used, because volume move reconstructs the whole tree
    if (lkdtree) then
        do ibox = 1, nbox
            if (lkdtree_box(ibox)) then
                call update_tree_height(mol_tree(ibox)%tree)
                !< reconstruct the tree if the height gets too large
                if (mol_tree(ibox)%tree%height .ge. (2.0*tree_height(ibox))) call construct_kdtree(ibox, ibox, .false.)
            end if
        end do
    end if

    ! Optimize and output MC move parameters
    if (ANY(lopt_bias).and.mod(nnn,freq_opt_bias).eq.0) then
       call opt_bias()
    end if

    if (mod(nnn,iratio).eq.0) then
       call update_translation_rotation_max_displacement(io_output)

       ! adjust maximum charge displacement for fluc Q
       lfq = .false.
       do im = 1,nbox
          do i = 1,nmolty
             if ( bnflcq(i,im) .gt. 0.5E0_dp ) then
                lfq = .true.
                ratflcq = bsflcq(i,im)/(bnflcq(i,im)*taflcq)
                if ( ratflcq .lt. 0.1E0_dp ) then
                   rmflcq(i,im) = rmflcq(i,im) * 0.1E0_dp
                else
                   rmflcq(i,im) = rmflcq(i,im) * ratflcq
                end if
             end if
             ! accumulate flcq info for final output
             bsflcq2(i,im) = bsflcq2(i,im) + bsflcq(i,im)
             bnflcq2(i,im) = bnflcq2(i,im) + bnflcq(i,im)
             ! rezero flcq
             bsflcq(i,im) = 0.0E0_dp
             bnflcq(i,im) = 0.0E0_dp
          end do
       end do
       if ( lfq.and.myid.eq.rootid ) then
          ! write out information about fluctuating charge success
          write(io_output,*) 'Box:   rmflcq for moltyps'
          do im =1,nbox
             write(io_output,*) im,(rmflcq(i,im),i=1,nmolty)
          end do
       end if
    end if

    if (mod(nnn,iupdatefix).eq.0) then
       call opt_safecbmc()
    end if

    if ((lgibbs.or.lnpt).and.(mod(nnn,iratv).eq.0).and.(pmvol.gt.0.0E0_dp)) then
       call update_volume_max_displacement(io_output)
    end if
! -------------------------------------------------------------------
    ! Calculate pressure and other related properties
    if (mod(nnn,iratp).eq.0) then
       ! calculate pressure ***
       acnp = acnp + 1
       do ibox = 1, nbox
          call pressure( press1, surf, comp, ibox )
          pres(ibox) = press1
          compress(ibox) = comp
          call update_average(acpres(ibox),press1,acnp)
          call update_average(acsurf(ibox),surf,acnp)
          call update_average(accomp(ibox),comp,acnp)

          if (lsolid(ibox) .and. .not. lrect(ibox)) then
             temvol = cell_vol(ibox)
          else
             if ( lpbcz ) then
                temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
             else
                temvol = boxlx(ibox)*boxly(ibox)
             end if
          end if

          ! Enthalpy calculation
          Temp_nmol = sum(ncmt(ibox,1:nmolty))
          ! molar volume in m3/mol, energies in kJ/mol
          Temp_Mol_Vol = temvol/Temp_nmol*N_Avogadro*1E-30_dp
          Temp_Energy  = vbox(ivTot,ibox)/Temp_nmol*R_gas/1000_dp
          call update_average(acEnthalpy(ibox),Temp_Energy+pres(ibox)*Temp_Mol_Vol,acnp)
          call update_average(acEnthalpy1(ibox),Temp_Energy+(express(ibox)*1E3_dp/MPa2SimUnits)*Temp_Mol_Vol,acnp)
       end do

       ! cannot calculate a heat of vaporization for only one box,
       ! and some compilers choke because Heat_vapor_T will not be
       ! defined if nbox == 1
       if (lgibbs) then
          do ibox = 1,nbox-1
             do jbox = ibox+1,nbox
                ! WRITE(io_output,*) 'ieouwfe ',ibox,jbox
                call calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ,Heat_vapor_COUL,pdV,CED_T,CED_LJ,CED_COUL,HSP_T,HSP_LJ&
                 ,HSP_COUL,ibox,jbox)

                ! Heat of vaporization
                call update_average(acsolpar(1,ibox,jbox),Heat_vapor_T,acnp)
                call update_average(acsolpar(2,ibox,jbox),Heat_vapor_LJ,acnp)
                call update_average(acsolpar(3,ibox,jbox),Heat_vapor_COUL,acnp)
                call update_average(acsolpar(4,ibox,jbox),CED_T,acnp)
                call update_average(acsolpar(5,ibox,jbox),CED_LJ,acnp)
                call update_average(acsolpar(6,ibox,jbox),CED_COUL,acnp)
                call update_average(acsolpar(7,ibox,jbox),HSP_T,acnp)
                call update_average(acsolpar(8,ibox,jbox),HSP_LJ,acnp)
                call update_average(acsolpar(9,ibox,jbox),HSP_COUL,acnp)
                !call update_average(acsolpar(10,ibox,jbox),DeltaU_Ext,acnp)
                call update_average(acsolpar(11,ibox,jbox),pdV,acnp)
             end do
          end do
       end if
    end if
! -------------------------------------------------------------------
    ! Print out summary current simulation status
    if ((mod(nnn,iprint).eq.0).and.(myid.eq.rootid)) then
       ! write out runtime information ***
       write(io_output,FMT='(i6,i8,e12.4,f10.3,f12.1,f12.2,'//format_n(nmolty,"i12")//')') nnn,tmcc,vbox(ivTot,1),boxlx(1),pres(1)&
        ,compress(1),(ncmt(1,imolty),imolty=1,nmolty)
       if ( lgibbs ) then
          do ibox = 2, nbox
             write(io_output,FMT='(14x,e12.4,f10.3,f12.1,f12.2,'//format_n(nmolty,"i12")//')') vbox(ivTot,ibox),boxlx(ibox),pres(ibox)&
              ,compress(ibox),(ncmt(ibox,imolty),imolty=1,nmolty)
          end do
       end if
    end if

    if((lgibbs.or.lnpt.or.lgrand)&
      .and.(.not.lvirial).and.(myid.eq.rootid).and.ltraj) then
       do ibox = 1,nbox
          if ( lpbcz ) then
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                write(io_cell,'(i8,6f12.4)') tmcc,cell_length(ibox,1),cell_length(ibox,2),cell_length(ibox,3)&
                 ,cell_ang(ibox,1)*raddeg,cell_ang(ibox,2)*raddeg,cell_ang(ibox,3)*raddeg
                write(io_traj,FMT='(8E13.5,'//format_n(nmolty,"i5")//')') hmat(ibox,1),hmat(ibox,4),hmat(ibox,5),hmat(ibox,7)&
                 ,hmat(ibox,8),hmat(ibox,9),vbox(ivTot,ibox),pres(ibox),(ncmt(ibox,itype),itype=1,nmolty)
             else
                write(io_traj,'(5E13.5,'//format_n(nmolty,"i5")//')') boxlx(ibox),boxly(ibox),boxlz(ibox),vbox(ivTot,ibox)&
                 ,pres(ibox),(ncmt(ibox,itype),itype=1,nmolty)
             end if
          else
             write(io_traj,'(3E12.5,'//format_n(nmolty,"i4")//')') boxlx(ibox)*boxly(ibox),vbox(ivTot,ibox),pres(ibox)&
              ,(ncmt(ibox,itype),itype=1,nmolty)
          end if
       end do
    end if
! -------------------------------------------------------------------
    ! Write out movie frames
    if (mod(nnn,imv).eq.0) then
       if ( lvirial ) then
          call virial(binvir,binvir2)
       else if (myid.eq.rootid) then
          ! write out the movie configurations ***
          write(io_movie,*) nnn
          do ibox = 1, nbox
             if (ibox.ne.gcbmc_box_num .or. l_gcbmc_movie) then
                 !< write out box information unless this is the reservoir box and l_gcbmc_movie is false
                 write(io_movie,*) (ncmt(ibox,zzz),zzz=1,nmolty)

                 if (lsolid(ibox) .and. .not. lrect(ibox)) then
                    write(io_movie,*) (hmat(ibox,zzz),zzz=1,9)
                 else
                    write(io_movie,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
                 end if
             end if
          end do
          do m = 1, nchain
             imolty = moltyp(m)
             if (l_gcbmc_movie .or. nboxi(m).ne.gcbmc_box_num) then
                 !< only write if the molecule is not in the reservoir box, or l_gcbmc_movie is true
                 write(io_movie,'(4(1x,i5),3(1x,f16.6))') m,imolty,nunit(imolty),nboxi(m),xcm(m),ycm(m),zcm(m)
                 do mm = 1, nunit(imolty)
                     write(io_movie,'(4(1x,f14.6),i5)') rxu(m,mm),ryu(m,mm),rzu(m,mm),qqu(m,mm),atoms%list(ntype(imolty,mm))
                 end do
             end if
          end do

          ! KM for MPI
          if (L_movie_xyz) then
             do ibox = 1,nbox
                if (l_gcbmc_movie .or. ibox.ne.gcbmc_box_num) then
                    nummol = 0
                    do i = 1,nchain
                       if (nboxi(i).eq.ibox) then
                          nummol = nummol + nunit(moltyp(i))
                       end if
                    end do
                    write(io_box_movie(ibox),*) nummol
                    write(io_box_movie(ibox),*)
                    do i = 1,nchain
                       if(nboxi(i).eq.ibox) then
                          imolty = moltyp(i)
                          do ii = 1,nunit(imolty)
                             ntii = ntype(imolty,ii)
                             write(io_box_movie(ibox),'(a4,5x,3f15.4)') chemid(ntii), rxu(i,ii), ryu(i,ii), rzu(i,ii)
                          end do
                       end if
                    end do
                 end if
             end do
          end if

          if (L_movie_pdb) then
              do ibox = 1,nbox
                  call writePDBmovie(io_box_movie_pdb(ibox),ibox)
              end do
          end if

       end if
    end if
! -------------------------------------------------------------------
    ! calculate the integrand of thermosynamic integration
    if (lmipsw.and.(mod(nnn,iratipsw).eq.0)) then
       acipsw = acipsw+1
       call deriv(1)
       call update_average(acdvdl,dvdl,acipsw)
    end if

    ! JLR 11-11-09
    ! do not call analysis if ianalyze is greater than number of cycles
    ! KM 01/10 remove analysis
    !if(mod(nnn,ianalyze).eq.0) then
    !call analysis(1)
    ! end if
    ! END JLR 11-11-09

    do intg = 1, nchain
       ibox = nboxi(intg)
       imolty = moltyp(intg)
       ! accumulate m-n-box and m-s-e-t-e-l ***
       ! only count the main chain - not the hydrogens
       ilunit = nugrow(imolty)
       setx = rxu(intg,1) - rxu(intg,ilunit)
       sety = ryu(intg,1) - ryu(intg,ilunit)
       setz = rzu(intg,1) - rzu(intg,ilunit)
       setel = setx*setx + sety*sety + setz*setz
       mnbox(ibox,imolty) = mnbox(ibox,imolty) + 1
       call update_average(asetel(ibox,imolty),setel,mnbox(ibox,imolty))
    end do

    do imolty = 1, nmolty
       ! KM for MPI
       ! only do anything for lsolute if monper is not called from readdat (nnn.ne.0)
       ! calculate energy and write out movie for lsolute
       if ((nnn.ne.0).and.mod(nnn,isolute(imolty)).eq.0) then
          do ibox = 1, nbox
             do k = 1, ncmt(ibox,imolty)
                i = parbox(k,ibox,imolty)
                ! set coords for energy and write out conformations
                if (myid.eq.rootid) write(io_solute,*) imolty,ibox,nunit(imolty)
                do j = 1, nunit(imolty)
                   rxuion(j,1) = rxu(i,j)
                   ryuion(j,1) = ryu(i,j)
                   rzuion(j,1) = rzu(i,j)
                   if (myid.eq.rootid) write(io_solute,*) ntype(imolty,j),rxuion(j,1),ryuion(j,1),rzuion(j,1),qqu(j,1)
                end do

                call energy(i,imolty,v,1,ibox,1,nunit(imolty),.true.,ovrlap,.false.,.false.,.false.,.false.)
                if (ovrlap) write(io_output,*)  '*** DISASTER, OVERLAP IN MONPER'

                if (ltailc) then
                   ! tail corrections
                   if (lsolid(ibox).and..not.lrect(ibox)) then
                      vol = cell_vol(ibox)
                   else
                      vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
                   end if
                   v(ivTail) = 0.0E0_dp
                   do jmolty = 1, nmolty
                      rho = ncmt(ibox,jmolty) / vol
                      v(ivTail) = v(ivTail) + ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
                   end do
                end if

                call U_bonded(i,imolty,v(ivStretching),v(ivBending),v(ivTorsion))

                solcount(ibox,imolty) = solcount(ibox,imolty) + 1
                call update_average(avsolinter(ibox,imolty),v(ivInterLJ)/2.0_dp+v(ivTail),solcount(ibox,imolty))
                call update_average(avsolintra(ibox,imolty),v(ivIntraLJ),solcount(ibox,imolty))
                call update_average(avsoltor(ibox,imolty),v(ivTorsion),solcount(ibox,imolty))
                call update_average(avsolbend(ibox,imolty),v(ivBending),solcount(ibox,imolty))
                call update_average(avsolelc(ibox,imolty),v(ivElect)+v(ivEwald),solcount(ibox,imolty))
             end do
          end do
       end if
    end do
! -------------------------------------------------------------------
    ! calculation of block averages
    if (mod(nnn,iblock).eq.0)  then
       nblock = nblock + 1
       do ibox=1,nbox
          ! specific density
          temmass = 0.0E0_dp
          if ( lpbcz ) then
             do itype = 1, nmolty
                temmass = temmass + masst(itype)*acdens(ibox,itype)
             end do
             dpr = temmass*1E24_dp/N_Avogadro
          else
             do itype = 1, nmolty
                temmass = temmass + acdens(ibox,itype)
             end do
             dpr = temmass
          end if
          call store_block_average(baver(1,ibox,nblock),dpr,acmove,bccold(1,ibox),nccold(1,ibox))

          ! pressure
          call store_block_average(baver(2,ibox,nblock),acpres(ibox),acnp,bccold(2,ibox),nccold(2,ibox))

          ! surface tension
          itel = nEnergy + 4*nmolty + 3
          call store_block_average(baver(itel,ibox,nblock),acsurf(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))

          ! compressibility factor
          itel = nEnergy + 4*nmolty + 7
          call store_block_average(baver(itel,ibox,nblock),accomp(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))

          ! box volume
          itel = nEnergy + 4*nmolty + 4
          call store_block_average(baver(itel,ibox,nblock),acvolume(ibox),acmove,bccold(itel,ibox),nccold(itel,ibox))

          ! energies
          do j=1,nEnergy
             itel = 2 + j
             call store_block_average(baver(itel,ibox,nblock),acv(j,ibox),acmove,bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! chemical potential
          do itype = 1, nmolty
             itel = 2 + nEnergy + itype
             dpp = acchem(ibox,itype)/bnchem(ibox,itype)
             call store_block_average(baver(itel,ibox,nblock),dpp,bnchem(ibox,itype),bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! square end-to-end length
          do itype = 1, nmolty
             itel = 2 + nEnergy + nmolty + itype
             call store_block_average(baver(itel,ibox,nblock),asetel(ibox,itype),mnbox(ibox,itype),bccold(itel,ibox)&
              ,nccold(itel,ibox))
          end do

          ! number density
          do itype = 1, nmolty
             itel = 2 + nEnergy + 2*nmolty + itype
             call store_block_average(baver(itel,ibox,nblock),acdens(ibox,itype),acmove,bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! mol fraction
          do itype = 1, nmolty
             itel = 2 + nEnergy + 3*nmolty + itype
             call store_block_average(baver(itel,ibox,nblock),molfra(ibox,itype),acmove,bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! Enthalpy
          itel = nEnergy + 4*nmolty + 5
          call store_block_average(baver(itel,ibox,nblock),acEnthalpy(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))
          itel = nEnergy + 4*nmolty + 6
          call store_block_average(baver(itel,ibox,nblock),acEnthalpy1(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))
       end do
       do ibox = 1,nbox-1
          do jbox = ibox+1,nbox
             do i = 1,nprop1
                call store_block_average(baver1(i,ibox,jbox,nblock),acsolpar(i,ibox,jbox),acnp,bccold1(i,ibox,jbox)&
                 ,nccold1(i,ibox,jbox))
             end do
          end do
       end do

       if (lucall) call blk_avg_prop_widom(nblock)
    end if

    if (ldielect.and.(mod(nnn,idiele).eq.0).and.myid.eq.rootid) then
       ! use fort.27 to calculate dielectric constant
       !---------not currently used---------!
       ! If you really want this quantity comment should be taken out**
       ! dielect = 6.9994685465110493E5_dp*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
       ! write(14,*) nnn,dielect*acdipolesq(4,ibox)/acmove
       ! write(15,*) nnn,dielect*(acdipolesq(4,ibox)/acmove-(acdipole(1,ibox)/acmove)**2-(acdipole(2,ibox)/acmove)**2 - (acdipole(3,ibox)/acmove)**2)
       ! write(16,*) nnn,acdipole(1,ibox)/acmove,acdipole(2,ibox)/acmove,acdipole(3,ibox)/acmove
       ! write(25,*) nnn+nnstep,vbox(1,1)
       !---------not currently used---------!

       ! output the fluctuation information
       write(14,*) nnn,acvol(ibox)
       write(15,*) nnn,acvolsq(ibox)
       write(16,*) nnn,acv(1,ibox)
       write(17,*) nnn,acvsq(1,ibox)
       write(18,*) nnn,acvol(ibox)*acv(1,ibox)
       write(19,*) nnn,acv(14,ibox)-acvol(ibox)*acv(1,ibox)
       do ibox = 1,nbox
          write(27,*) dipolex(ibox),dipoley(ibox),dipolez(ibox)
       end do
    end if

    ! ibox = 1
    ! imolty = 1
    ! do i = 1,nchain
    !    if ( nboxi(i) .eq. ibox ) then
    !       if ( moltyp(i) .eq. imolty ) then
    !          bin = aint(zcm(i)/binstep) + 1
    !          temvol = boxlx(ibox)*boxly(ibox)*binstep
    !          profile(bin) = profile(bin)+1.0E0_dp/temvol
    !       end if
    !    end if
    ! end do

    ! Output running averages of residual heat capacity, in (J2/mol)
    if (mod(nnn,iheatcapacity).eq.0.and..not.lgibbs.and.myid.eq.rootid) then
       write(56,'(I12,F18.6,F18.2,F18.6)') nnn,(enthalpy2-enthalpy*enthalpy)/real(nchain,dp)*R_gas/(temp**2),enthalpy2,enthalpy
    end if

    ! check rigid structures
    call check_rigid_structures()

#ifdef __DEBUG__
    write(io_output,*) 'end MONPER in ',myid
#endif
    return
  end subroutine monper

  subroutine read_checkpoint_main(file_chkpt)
    use util_mp,only:mp_bcast
    use moves_simple,only:read_checkpoint_simple
    use moves_cbmc,only:read_checkpoint_cbmc
    use moves_volume,only:read_checkpoint_volume
    use transfer_swap,only:read_checkpoint_swap
    use transfer_swatch,only:read_checkpoint_swatch
    use transfer_shared,only:read_checkpoint_transfer_shared
    use moves_ee,only:read_checkpoint_ee
    use prop_widom,only:read_checkpoint_prop_widom
    character(LEN=*),intent(in)::file_chkpt
    integer::io_chkpt,jerr,i,pos_output,pos_movie,pos_box_movie(nbox),pos_solute,pos_cell,pos_traj

    if (myid.eq.rootid) then
       io_chkpt=get_iounit()
       open(unit=io_chkpt,access='stream',action='read',file=file_chkpt,form='unformatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot read checkpoint file '//trim(file_chkpt),jerr)

       read(io_chkpt) pos_output,pos_movie,pos_box_movie,pos_solute,pos_cell,pos_traj,nnstep,nnn,acmove,vbox,vstart,acdens,molfra&
        ,acnbox,acnbox2,acv,acvsq,acvkjmol,acdipole,acdipolesq,acboxa,acboxl,acvol,acvolsq,acvolume,enthalpy,enthalpy2,acnp,acpres&
        ,acsurf,accomp,acEnthalpy,acEnthalpy1,solcount,acsolpar,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc,mnbox,asetel,acipsw&
        ,acdvdl,nblock,nccold1,bccold1,baver1,nccold,bccold,baver
    end if

    if (lucall) call read_checkpoint_prop_widom(io_chkpt,blockm)
    call read_checkpoint_simple(io_chkpt)
    call read_checkpoint_cbmc(io_chkpt)
    call read_checkpoint_volume(io_chkpt)
    call read_checkpoint_swap(io_chkpt)
    call read_checkpoint_swatch(io_chkpt)
    call read_checkpoint_transfer_shared(io_chkpt)

    ! fluctuating charges
    if (myid.eq.rootid) read(io_chkpt) bnflcq,bsflcq,bnflcq2,bsflcq2

    call read_checkpoint_ee(io_chkpt)

    call mp_bcast(nnstep,1,rootid,groupid)
    call mp_bcast(nnn,1,rootid,groupid)
    call mp_bcast(acmove,1,rootid,groupid)
    call mp_bcast(vbox,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(vstart,nbxmax,rootid,groupid)
    call mp_bcast(acdens,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(molfra,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acnbox,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acnbox2,nbxmax*ntmax*20,rootid,groupid)
    call mp_bcast(acv,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(acvsq,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(acvkjmol,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(acdipole,4*nbxmax,rootid,groupid)
    call mp_bcast(acdipolesq,4*nbxmax,rootid,groupid)
    call mp_bcast(acboxa,nbxmax*3,rootid,groupid)
    call mp_bcast(acboxl,nbxmax*3,rootid,groupid)
    call mp_bcast(acvol,nbxmax,rootid,groupid)
    call mp_bcast(acvolsq,nbxmax,rootid,groupid)
    call mp_bcast(acvolume,nbxmax,rootid,groupid)
    call mp_bcast(enthalpy,1,rootid,groupid)
    call mp_bcast(enthalpy2,1,rootid,groupid)
    call mp_bcast(acnp,1,rootid,groupid)
    call mp_bcast(acpres,nbxmax,rootid,groupid)
    call mp_bcast(acsurf,nbxmax,rootid,groupid)
    call mp_bcast(accomp,nbxmax,rootid,groupid)
    call mp_bcast(acEnthalpy,nbxmax,rootid,groupid)
    call mp_bcast(acEnthalpy1,nbxmax,rootid,groupid)
    call mp_bcast(solcount,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acsolpar,nprop1*nbxmax*nbxmax,rootid,groupid)
    call mp_bcast(avsolinter,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsolintra,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsolbend,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsoltor,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsolelc,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(mnbox,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(asetel,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acipsw,1,rootid,groupid)
    call mp_bcast(acdvdl,1,rootid,groupid)
    call mp_bcast(nblock,1,rootid,groupid)
    call mp_bcast(nccold1,nprop1*nbxmax*nbxmax,rootid,groupid)
    call mp_bcast(bccold1,nprop1*nbxmax*nbxmax,rootid,groupid)
    call mp_bcast(baver1,nprop1*nbxmax*nbxmax*blockm,rootid,groupid)
    call mp_bcast(nccold,nprop*nbxmax,rootid,groupid)
    call mp_bcast(bccold,nprop*nbxmax,rootid,groupid)
    call mp_bcast(baver,nprop*nbxmax*blockm,rootid,groupid)
    call mp_bcast(bnflcq,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsflcq,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bnflcq2,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsflcq2,ntmax*nbxmax,rootid,groupid)

    if (myid.eq.rootid) then
       close(io_chkpt)
       write(UNIT=io_output,FMT="()",ADVANCE='NO',POS=pos_output)
       if (io_movie.ge.0) write(UNIT=io_movie,FMT="()",ADVANCE='NO',POS=pos_movie)
       do i=1,nbox
          if (io_box_movie(i).ge.0) write(UNIT=io_box_movie(i),FMT="()",ADVANCE='NO',POS=pos_box_movie(i))
       end do
       if (io_solute.ge.0) write(UNIT=io_solute,FMT="()",ADVANCE='NO',POS=pos_solute)
       if (io_cell.ge.0) write(UNIT=io_cell,FMT="()",ADVANCE='NO',POS=pos_cell)
       if (ltraj) write(UNIT=io_traj,FMT="()",ADVANCE='NO',POS=pos_traj)
    end if
  end subroutine read_checkpoint_main

  subroutine write_checkpoint_main(file_chkpt)
    use util_files,only:flush_force
    use moves_simple,only:write_checkpoint_simple
    use moves_cbmc,only:write_checkpoint_cbmc
    use moves_volume,only:write_checkpoint_volume
    use transfer_swap,only:write_checkpoint_swap
    use transfer_swatch,only:write_checkpoint_swatch
    use transfer_shared,only:write_checkpoint_transfer_shared
    use moves_ee,only:write_checkpoint_ee
    use prop_widom,only:write_checkpoint_prop_widom
    character(LEN=*),intent(in)::file_chkpt
    integer::io_chkpt,jerr,i,pos_output,pos_movie,pos_box_movie(nbox),pos_solute,pos_cell,pos_traj

    pos_output=flush_force(io_output)
    pos_movie=flush_force(io_movie)
    do i=1,nbox
       pos_box_movie(i)=flush_force(io_box_movie(i))
    end do
    pos_solute=flush_force(io_solute)
    pos_cell=flush_force(io_cell)
    pos_traj=flush_force(io_traj)

    io_chkpt=get_iounit()
    open(unit=io_chkpt,access='stream',action='write',file=file_chkpt,form='unformatted',iostat=jerr,status='unknown')
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open checkpoint file '//trim(file_chkpt),jerr)

    write(io_chkpt) pos_output,pos_movie,pos_box_movie,pos_solute,pos_cell,pos_traj,nnstep,nnn,acmove,vbox,vstart,acdens,molfra&
     ,acnbox,acnbox2,acv,acvsq,acvkjmol,acdipole,acdipolesq,acboxa,acboxl,acvol,acvolsq,acvolume,enthalpy,enthalpy2,acnp,acpres&
     ,acsurf,accomp,acEnthalpy,acEnthalpy1,solcount,acsolpar,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc,mnbox,asetel,acipsw,acdvdl&
     ,nblock,nccold1,bccold1,baver1,nccold,bccold,baver

    if (lucall) call write_checkpoint_prop_widom(io_chkpt)
    call write_checkpoint_simple(io_chkpt)
    call write_checkpoint_cbmc(io_chkpt)
    call write_checkpoint_volume(io_chkpt)
    call write_checkpoint_swap(io_chkpt)
    call write_checkpoint_swatch(io_chkpt)
    call write_checkpoint_transfer_shared(io_chkpt)

    ! fluctuating charges
    write(io_chkpt) bnflcq,bsflcq,bnflcq2,bsflcq2

    call write_checkpoint_ee(io_chkpt)

    close(io_chkpt)
  end subroutine write_checkpoint_main

#ifdef w_nl
#undef w_nl
#endif

#ifdef wa_nl
#undef wa_nl
#endif

#ifdef __TRADITIONAL_CPP__
#define w_nl(x) write(io_unit,*) '    ',"x",'=',x
#define wa_nl(x,n) if (n.gt.0) write(io_unit,*) '    ',"x",'=',x(1:n)
#else
#define w_nl(x) write(io_unit,*) '    ',#x,'=',x
#define wa_nl(x,n) if (n.gt.0) write(io_unit,*) '    ',#x,'=',x(1:n)
#endif
!> \brief Generate standard input files containing all parameters
!>
!> Put the following calling statement right after having finished reading
!> fort.4 in readdat
!> call generate_standard_input(lionic,L_Ewald_Auto,lmixlb,lmixjo,seed,linit,lreadq,fqtemp,inclnum,inclmol,inclbead,inclsign,ncarbon,ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)
  subroutine generate_standard_input(lionic,L_Ewald_Auto,lmixlb,lmixjo,lmixwh,lmixkong,seed,linit,lreadq,fqtemp,inclnum,inclmol,inclbead&
   ,inclsign,ncarbon,ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)
    use const_phys,only:MPa2SimUnits
    use util_string,only:format_n
    use energy_intramolecular,only:bonds,angles,dihedrals
    integer,parameter::io_unit=512
    logical,intent(in)::lionic,L_Ewald_Auto,lmixlb,lmixjo,lmixwh,lmixkong,linit,lreadq
    integer,intent(in)::seed,inclnum,ainclnum,inclmol(:),inclbead(:,:),inclsign(:),ncarbon(:),ainclmol(:),ainclbead(:,:),a15t(:)
    real,intent(in)::fqtemp,ofscale(:),ofscale2(:)
    integer::i,imol,j,avbmc_version(nmolty)
    real::rmvolume,rmtra,rmrot,armtra,rmflucq
! ===================================================================
    open(unit=io_unit,access='sequential',action='write',file='topmon.inp.std',form='formatted',status='new')
! -------------------------------------------------------------------
    write(io_unit,'(/," &io")')
    w_nl(io_output)
    w_nl(run_num)
    w_nl(suffix)
    w_nl(L_movie_xyz)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &system")')
    w_nl(lnpt)
    w_nl(lgibbs)
    w_nl(lgrand)
    w_nl(lanes)
    w_nl(lvirial)
    w_nl(lmipsw)
    w_nl(lexpee)
    w_nl(ldielect)
    w_nl(lpbc)
    w_nl(lpbcx)
    w_nl(lpbcy)
    w_nl(lpbcz)
    w_nl(lfold)
    w_nl(lijall)
    w_nl(lchgall)
    w_nl(lewald)
    w_nl(lcutcm)
    w_nl(ltailc)
    w_nl(lshift)
    w_nl(ldual)
    w_nl(L_Coul_CBMC)
    w_nl(lneigh)
    w_nl(lexzeo)
    w_nl(lslit)
    w_nl(lgraphite)
    w_nl(lsami)
    w_nl(lmuir)
    w_nl(lelect_field)
    w_nl(lgaro)
    w_nl(lionic)
    w_nl(L_Ewald_Auto)
    w_nl(lmixlb)
    w_nl(lmixjo)
    w_nl(lmixwh)
    w_nl(lmixkong)
    w_nl(L_spline)
    w_nl(L_linear)
    w_nl(L_vib_table)
    w_nl(L_bend_table)
    w_nl(L_elect_table)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    close(io_unit)
! ===================================================================
    open(unit=io_unit,access='sequential',action='write',file='fort.4.std',form='formatted',status='new')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_shared")')
    w_nl(seed)
    w_nl(nbox)
    w_nl(nmolty)
    w_nl(nchain)
    w_nl(nstep)
    w_nl(lstop)
    w_nl(iratio)
    w_nl(rmin)
    w_nl(softcut)
    w_nl(linit)
    w_nl(lreadq)
    w_nl(N_add)
    w_nl(box2add)
    w_nl(moltyp2add)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &analysis")')
    w_nl(iprint)
    w_nl(imv)
    w_nl(iblock)
    w_nl(iratp)
    w_nl(idiele)
    w_nl(iheatcapacity)
    w_nl(ianalyze)
    w_nl(nbin)
    w_nl(lrdf)
    w_nl(lintra)
    w_nl(lstretch)
    w_nl(lgvst)
    w_nl(lbend)
    w_nl(lete)
    w_nl(lrhoz)
    w_nl(bin_width)
    w_nl(lucall)
    w_nl(nvirial)
    w_nl(starvir)
    w_nl(stepvir)
    w_nl(ntemp)
    w_nl(virtemp)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &external_field")')
    !wa_nl(Elect_field,nbox)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_volume")')
    w_nl(tavol)
    w_nl(iratv)
    wa_nl(pmvlmt,nbox)
    w_nl(nvolb)
    wa_nl(pmvolb,nvolb)
    wa_nl(box5,nvolb)
    wa_nl(box6,nvolb)
    w_nl(pmvol)
    w_nl(pmvolx)
    w_nl(pmvoly)
    rmvolume=rmvol(1)
    w_nl(rmvolume)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_swatch")')
    w_nl(pmswat)
    w_nl(nswaty)
    wa_nl(pmsatc,nswaty)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_swap")')
    w_nl(pmswap)
    wa_nl(pmswmt,nmolty)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_cbmc")')
    w_nl(rcutin)
    w_nl(pmcb)
    wa_nl(pmcbmt,nmolty)
    wa_nl(pmall,nmolty)
    wa_nl(nchoi1,nmolty)
    wa_nl(nchoi,nmolty)
    wa_nl(nchoir,nmolty)
    wa_nl(nchoih,nmolty)
    wa_nl(nchtor,nmolty)
    wa_nl(nchbna,nmolty)
    wa_nl(nchbnb,nmolty)
    wa_nl(icbdir,nmolty)
    wa_nl(icbsta,nmolty)
    w_nl(rbsmax)
    w_nl(rbsmin)
    do i=1,nmolty
       if (lavbmc1(i)) then
          avbmc_version(i)=1
       else if (lavbmc2(i)) then
          avbmc_version(i)=2
       else if (lavbmc3(i)) then
          avbmc_version(i)=3
       else
          avbmc_version(i)=0
       end if
    end do
    wa_nl(avbmc_version,nmolty)
    wa_nl(pmbias,nmolty)
    wa_nl(pmbsmt,nmolty)
    wa_nl(pmbias2,nmolty)
    wa_nl(pmfix,nmolty)
    wa_nl(lrig,nmolty)
    w_nl(lpresim)
    w_nl(iupdatefix)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_ee")')
    w_nl(pmexpc)
    wa_nl(pmeemt,nmolty)
    w_nl(pmexpc1)
    wa_nl(lexpand,nmolty)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_simple")')
    armtra=(Armtrax+Armtray+Armtraz)/3.0_dp
    w_nl(armtra)
    rmtra=(rmtrax(1,1)+rmtray(1,1)+rmtraz(1,1))/3.0_dp
    w_nl(rmtra)
    rmrot=(rmrotx(1,1)+rmroty(1,1)+rmrotz(1,1))/3.0_dp
    w_nl(rmrot)
    w_nl(tatra)
    w_nl(tarot)
    w_nl(pmtra)
    wa_nl(pmtrmt,nmolty)
    wa_nl(pmromt,nmolty)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/," &mc_flucq")')
    w_nl(taflcq)
    w_nl(fqtemp)
    rmflucq=rmflcq(1,1)
    w_nl(rmflucq)
    w_nl(pmflcq)
    wa_nl(pmfqmt,nmolty)
    wa_nl(lflucq,nmolty)
    wa_nl(lqtrans,nmolty)
    wa_nl(fqegp,nmolty)
    wa_nl(nchoiq,nbox)
    write(io_unit,'(" /",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"SIMULATION_BOX")')
    do i=1,nbox
       write(io_unit,'("! boxlx boxly boxlz rcut kalp rcutnn numDimensionIsIstropic lsolid lrect lideal ltwice temperature &
        &pressure",/,3(F8.3,1X),3(F6.3,1X),I0,1X,4(L1,1X),F7.2,1X,F7.1)') boxlx(i),boxly(i),boxlz(i),rcut(i),kalp(i),rcutnn(i)&
        ,numberDimensionIsIsotropic(i),lsolid(i),lrect(i),lideal(i),ltwice(i),temp,express(i)/MPa2SimUnits
       write(io_unit,'("! nchain_1 ... nchain_nmolty ghost_particles",/,'//format_n(nmolty,'(I0,1X)')//',I0)')&
        (ininch(j,i),j=1,nmolty),ghost_particles(i)
       ! use_linkcell is always .false., and rintramax is 0.0
       write(io_unit,'("! inix iniy iniz inirot inimix zshift dshift use_linkcell rintramax",/,5(I0,1X),F4.1,1X,F6.3,L1,1X,F3.1)')&
        inix(i),iniy(i),iniz(i),inirot(i),inimix(i),zshift(i),dshift(i),.false.,0.0_dp
    end do
    write(io_unit,'("END SIMULATION_BOX",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"MOLECULE_TYPE")')
    do imol=1,nmolty
       ! lsetup is always .false.
       write(io_unit,'("! nunit nugrow ncarbon maxcbmc maxgrow iring lelect lring lrigid lbranch lsetup lq14scale qscale &
        &iurot isolute ",/,6(I0,1X),7(L1,1X),F3.1,2(1X,I0),'//format_n(nbox,'(1X,F7.1)')//')') nunit(imol),nugrow(imol)&
        ,ncarbon(imol),nmaxcbmc(imol),maxgrow(imol),iring(imol),lelect(imol),lring(imol),lrigid(imol),lbranch(imol)&
        ,.false.,lq14scale(imol),qscale(imol),iurot(imol),isolute(imol)
       if (lrigid(imol)) write(io_unit,'("! n, growpoint_1 ... growpoint_n",/,I0,'//format_n(rindex(imol),'(1X,I0)')//')')&
        rindex(imol),(riutry(imol,i),i=1,rindex(imol))
       do i=1,nunit(imol)
          write(io_unit,'("! unit ntype leaderq",/,2(I0,1X),I0,/,"! stretching",/,I0)') i,atoms%list(ntype(imol,i))&
           ,leaderq(imol,i),invib(imol,i)
          do j=1,invib(imol,i)
             write(io_unit,'(I0,1X,I0)') ijvib(imol,i,j),bonds%list(itvib(imol,i,j))
          end do
          write(io_unit,'("! bending",/,I0)') inben(imol,i)
          do j=1,inben(imol,i)
             write(io_unit,'(2(I0,1X),I0)') ijben2(imol,i,j),ijben3(imol,i,j),angles%list(itben(imol,i,j))
          end do
          write(io_unit,'("! torsion",/,I0)') intor(imol,i)
          do j=1,intor(imol,i)
             write(io_unit,'(3(I0,1X),I0)') ijtor2(imol,i,j),ijtor3(imol,i,j),ijtor4(imol,i,j),dihedrals%list(ittor(imol,i,j))
          end do
       end do
    end do
    write(io_unit,'("END MOLECULE_TYPE",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"MC_SWATCH")')
    do i=1,nswaty
       write(io_unit,'("! moltyp1<->moltyp2 nsampos 2xncut",/,4(I0,1X),I0)') nswatb(i,1:2),nsampos(i),ncut(i,1:2)
       write(io_unit,'("! gswatc 2x(ifrom, iprev)",/,'//format_n(2*(ncut(i,1)+ncut(i,2))-1,'(I0,1X)')//',I0,/,"! splist")')&
        (gswatc(i,j,1:2*ncut(i,j)),j=1,2)
       do j=1,nsampos(i)
          write(io_unit,'(I0,1X,I0)') splist(i,j,1:2)
       end do
       write(io_unit,'("! nswtcb pmswtcb",/,I0,'//format_n(nswtcb(i),'(F6.2,1X)')//',/,"! box numbers")') nswtcb(i)&
        ,pmswtcb(i,1:nswtcb(i))
       do j=1,nswtcb(i)
          write(io_unit,'(I0,1X,I0)') box3(i,j),box4(i,j)
       end do
    end do
    write(io_unit,'("END MC_SWATCH",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"MC_SWAP")')
    do i=1,nmolty
       write(io_unit,'("! nswapb pmswapb",/,I0,'//format_n(nswapb(i),'(F6.2,1X)')//',/,"! box1 box2")') nswapb(i)&
        ,(pmswapb(i,j),j=1,nswapb(i))
       do j=1,nswapb(i)
          write(io_unit,'(I0,1X,I0)') box1(i,j),box2(i,j)
       end do
    end do
    write(io_unit,'("END MC_SWAP",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"SAFE_CBMC")')
    do imol=1,nmolty
       if (lrig(imol)) then
          write(io_unit,'(I0)') nrig(imol)
          if (nrig(imol).gt.0) then
             do i=1,nrig(imol)
                write(io_unit,'(I0,1X,I0)') irig(imol,i),frig(imol,i)
             end do
          else
             write(io_unit,'(I0,1X,I0)') nrigmin(imol),nrigmax(imol)
          end if
       end if
    end do
    write(io_unit,'("END SAFE_CBMC",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"INTERMOLECULAR_EXCLUSION",/,"! mol_1 bead_1 mol_2 bead_2")')
    write(io_unit,'("END INTERMOLECULAR_EXCLUSION",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"INTRAMOLECULAR_SPECIAL",/,"! inclmol inclbead_1 inclbead_2 inclsign ofscale ofscale2")')
    do i=1,inclnum
       write(io_unit,'(4(I0,1X),F5.1,1X,F5.1)') inclmol(i),inclbead(i,1),inclbead(i,2),inclsign(i),ofscale(i),ofscale2(i)
    end do
    write(io_unit,'("END INTRAMOLECULAR_SPECIAL",/)')
! -------------------------------------------------------------------
    write(io_unit,'(/,"INTRAMOLECULAR_OH15",/,"! ainclmol ainclbead_1 ainclbead_2 a15type")')
    do i=1,ainclnum
       write(io_unit,'(3(I0,1X),I0)') ainclmol(i),ainclbead(i,1),ainclbead(i,2) ,a15t(i)
    end do
    write(io_unit,'("END INTRAMOLECULAR_OH15",/)')
! -------------------------------------------------------------------
    close(io_unit)
! ===================================================================
  end subroutine generate_standard_input

  subroutine check_rigid_structures()
    use transfer_swap,only:compute_beg,beg
    integer::iunit,i,j,imolty
    real::initial_value,my_dx,my_dy,my_dz,my_dist
    ! test if all rigid molecules have same configuration
    ! initialize all vectors
    initial_value = 1000.0_dp
    rigid_intra_dist = initial_value
    do i=1,nchain
       imolty = moltyp(i)
       iunit = nunit(imolty)
       if (lrigid(imolty)) then
          call compute_beg(imolty)
          do j = beg+1,iunit ! beg+1, not all beads will be considered if growpoint
              my_dx = rxu(i,beg) - rxu(i,j)
              my_dy = ryu(i,beg) - ryu(i,j)
              my_dz = rzu(i,beg) - rzu(i,j)
              my_dist = my_dx*my_dx + my_dy*my_dy + my_dz*my_dz
              if (abs(rigid_intra_dist(j,imolty)-initial_value) < distance_tolerance) then
                  ! store distance
                  rigid_intra_dist(j,imolty) =  my_dist
              else if (abs(rigid_intra_dist(j,imolty) - my_dist) > distance_tolerance) then
                  ! error exit
                  write(io_output,*) 'stored distance: ',rigid_intra_dist(j,imolty), 'new vector: ', my_dx
                  write(io_output,*) 'imolty: ',imolty, 'bead: ',j
                  write(io_output,*) 'ichain: ',i
                  call err_exit(__FILE__,__LINE__,'rigid struc problem(x): not all structs the same ',myid+1)
              end if
          end do
       end if
    end do
  end subroutine check_rigid_structures
END MODULE topmon_main
