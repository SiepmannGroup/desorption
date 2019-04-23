module sim_system
  use var_type,only:dp,default_path_length,atom_symbol_length,long_int
  use util_runtime,only:err_exit
  use util_string,only:integer_to_string
  use util_files,only:get_iounit
  use util_search,only:LookupTable
  implicit none
  save

!==========================================================!
!=============== General Programming Advice ===============!
!==========================================================!
! Ideally, global variables should be kept to a MINIMUM.
! Variables below should go into their separate modules
! (e.g., MPI-related variables to util_mp,
! volume-moves-related variables to moves_volume). Exposing
! too many variables where they are not needed increases the
! chance that they are modified inadvertently.

! In many cases, things that you feel are "required"
! globally (for example when outputing swap statistics, you
! may want counters defined in transfer_swap available in
! the main program) can be delegated (calling the subroutine
! output_swap_stats, defined also in transfer_swap,
! eliminates the need to expose internal variables of
! transfer_swap)

! The ideas of software engineering are broad. Not every
! principle there can be conveniently implemented using
! Fortran. However, you should always try to look for a
! general solution than an ad-hoc fix. When problems you are
! dealing with grow in complexity (they inevitably will),
! you will find solutions that are general in nature and
! clear in logic are easier to extend and maintain.

! With computer programming, if you find yourself doing the
! same thing twice, odds are that there is a better way to
! it. Spend 10 minutes (even hours) to learn a technique
! that can be implemented in 1 minute, rather than spend 1
! minute to come up with a technique that takes 10 minutes.
! Even if the learning phase takes longer, it will usually
! pay off in the long run: you may likely apply the better
! method many times in the future and you have sharpened
! your programming skills.
!==========================================================!

  !=== Information about the system ===

  !*** PARAMETERS FOR ENSEMBLE ***
  logical::lnpt=.false.& !< if LNPT=.TRUE. then a NPT volume move is used to equilibrate with a pressure bath (implies cubic simulation boxes) else an NVT simulation is performed
   ,lgibbs=.false.& !< if LGIBBS=.TRUE. then a Gibbs-ensemble simulation is performed (implies cubic simulation boxes)
   ,lgrand=.false.& !< if LGRAND=.TRUE. then simulation is performed in the grand-canonical ensemble
   ,lanes=.false.& !< if LANES=.TRUE. then simulation is performed in the adiabatic nuclear and electronic sampling technique for polarizable force fields
   ,lvirial=.false.& !<! if LVIRIAL=.TRUE. then one chain will be simulated in each box independently and the second virial coefficient will be calculated for their interactions at a series of distances along the x-axis
   ,lmipsw=.false.& !< if lmipsw is true, then thermodynamic integration is performed for the phases, fort.35 must be supplied if so
   ,lexpee=.false.& !< if lexpee is true, then expanded esnemble is performed for the phases, fort.44 must be supplied if so
   ,ldielect=.false.& !< if LDIELECT=.TRUE. then dielectric constant will be calculated and LEWALD must be .TRUE. Correct only in NVT ensemble
   ,losmoticnvt=.false. !< True for osmotic NVT Gibbs simulations, which add a delta{P}*delta{V} term to the energy term in the Boltzmann factor

  !*** PARAMETERS FOR BOUNDARY CONDITIONS ***
  logical::lpbc=.true.& !< if LPBC = .TRUE. then periodic boundaries are used
   ,lpbcx=.true.& !< if LPBCX = .TRUE. then periodic boundary in x-directions is used
   ,lpbcy=.true.& !< if LPBCY = .TRUE. then periodic boundary in y-directions is used
   ,lpbcz=.true.& !< if LPBCZ = .TRUE. then periodic boundary in z-directions is used
   ,lfold=.true. !< if LFOLD = .TRUE. then coordinates are always folded into central box

  !*** PARAMETERS FOR INTERACTION CUT-OFFS ***
  logical::lijall=.false.& !< if LIJALL = .TRUE. then all i-j interactions are considered (no potential cut-off). Must set lcutcm and ltailc to false for lijall = true *** check top of sumup.f if lijall = true and lchgall is false
   ,lchgall=.false.& !< if LCHGALL= .TRUE. then all the electrostatic interaction are considered
   ,lewald=.false.& !< if LEWALD=.TRUE. then ewald-sum will be used to calculate the electrostatic interactions.
   ,lrecip=.false.& !< if LRECIP=.TRUE. then sin and cos in reciprocal space will be calculated recursively.
   ,lcutcm=.true.& !< if LCUTCM=.TRUE. then a cutoff of the centers of mass will be used with a value of rcmu as calculated in ctrmas
   ,L_Ewald_Auto=.true.

  !*** PARAMETERS FOR TAIL CORRECTIONS ***
  logical::ltailc=.true.& !< if LTAILC=.TRUE. tail corrections are added (WARNING: .lsami. in external.inc switches an intrinsic tail correction on)
   ,lshift=.false. !< truncated and shifted potentials

  !*** PARAMETERS FOR CBMC-PSEUDO-POTENTIAL ***
  logical::ldual=.true.& !< if LDUAL=.TRUE. then the external potential during a CBMC growth will only go out to a radius of rcutin and then will be corrected to the full rcut at the end. This is Dual Cutoff Configurational-bias Monte Carlo (DC-CBMC)
   ,L_Coul_CBMC=.true. !< if L_Coul_CBMC=.TRUE. then the electrostatic interaction is computed during CBMC/SWAP

  !*** PARAMETERS OF NEIGHBOR LIST ***
  logical::lneigh=.false.& !< if LNIEGH=.TRUE. the nearest neighbor list will be used with a value of rcutnn specified by fort.4
   ,lneighbor

  !*** External surface ***
  logical::lexzeo=.false.& !< implicit rigid framework for zeolites and metal-organic frameworks. Tabulated potential will be used for the interactions with the rigid framework. Parameters are read in subroutine SUZEO
   ,lslit=.false.& !< featureless, planar slit surface(s)
   ,lgraphite=.false.& !< x,y-dependent graphite surface
   ,lsami=.false.& !< SAMI
   ,lmuir=.false.& !< Langmuir monolayer
   ,lelect_field=.false. !< external electric field

  !*** Thermodynamic integration ***
  real,allocatable::rxwell(:,:),rywell(:,:),rzwell(:,:),sxwell(:,:),sywell(:,:),szwell(:,:)&
   ,vwellipswot(:),vwellipswnt(:),vipswnt(:),vipswot(:),awell(:,:,:)
  integer,allocatable::nwell(:)
  logical,allocatable::lwell(:)
  integer,parameter::nw=4000
  integer::iratipsw=100000000
  real::dvdl,vipsw,pipsw,vwellipsw,pwellipsw,etais,lambdais,bwell,vipswo,vipswn,vwellipswo,vwellipswn&
   ,lena,lenc,pwellips(3,3),pips(3,3),dhmat(3,3)
  logical::lstagea,lstageb,lstagec

  !=== Analysis ===
  integer,parameter::nEnergy=14,ivTot=1,ivInterLJ=2,ivTail=3,ivIntraLJ=4,ivStretching=5,ivBending=6,ivTorsion=7,ivElect=8,ivExt=9&
   ,iv3body=10,ivFlucq=11,ivIpswb=12,ivWellIpswb=13,ivEwald=14
  real,allocatable::vbox(:,:) !< (j,ibox): energies of ibox;
                              !< j = 1: total energy;
                              !< 2: intermolecular LJ;
                              !< 3: tail correction;
                              !< 4: intramolecular non-bonded LJ;
                              !< 5: stretching; 6: bending; 7: torsion;
                              !< 8: electrostatic;
                              !< 9: external field;
                              !< 10: 3-body Feuston-Garofalini;
                              !< 11: fluctuating charge (vflucqb);
                              !< 12: vipswb; 13: vwellipswb;
                              !< 14: Ewald reciprocal-space electrostatic
  character(LEN=1)::suffix='a'
  character(LEN=5)::run_num='1'
  integer::io_output=6,iprint=10000000,imv=10000000,iblock=10000000,iratp=500,idiele=10000000,iheatcapacity=10000000,nprop,ianalyze=10000000,nbin=1
  real::bin_width=0.2_dp
  logical::L_movie_xyz=.false.,lrdf=.false.,lintra=.false.,lstretch=.false.,lgvst=.false.,lbend=.false.,lete=.false.&
   ,lrhoz=.false.,lucall=.false.,L_movie_pdb=.false.,ltraj=.true.

  !*** Histograms for grand-canonical ensemble ***
  integer::nequil=0,ninstf=0,ninsth=0,ndumph=0

  !*** 2nd virial coefficient ***
  integer,parameter::maxvir=1,maxntemp=1 !< maximum number of bins for the 2nd virial coefficient
  integer::ntemp=0,nvirial=0
  real::virtemp,starvir=0.0_dp,stepvir=0.0_dp

  !*** Dielectric constant ***
  real,allocatable::dipolex(:),dipoley(:),dipolez(:)
  real::dipolexo,dipoleyo,dipolezo

  !*** RPLC ***
  logical,allocatable::lrplc(:)& !< if lrplc=.true. there are some special rules in CBMC for how to grow chains
   ,ltwice(:) !< if ltwice=.true. then mimage is applied twice

  !=== Force field parameters ===
  integer,parameter::nvib_max=6,nben_max=12,ntor_max=12,max_num_rigid_beads=24
  integer::nntype !< number of types of beads
  logical::L_spline=.false.,L_linear=.false.,L_vib_table=.false.,L_bend_table=.false.,L_elect_table=.false.
  real,allocatable::mass(:),vvdW_b(:,:),qelect(:),vvdW(:,:),ecut(:)&
   ,brvib(:),brvibk(:),brben(:),brbenk(:),maxRegrowVib(:),minRegrowVib(:)&
   ,ljscale(:,:,:),qscale2(:,:,:)
  character(len=atom_symbol_length),allocatable::chemid(:)
  integer,allocatable::a15type(:,:,:)
  logical,allocatable::lij(:),lqchg(:),lexclu(:,:,:,:),linclu(:,:,:),lqinclu(:,:,:),lainclu(:,:,:),lpl(:),lexclu_zeo(:)

  !*** Feuston-Garofalini force field ***
  logical::lgaro !< if LGARO=.TRUE. Feuston-Garofalini potential will be used

  !=== Information about simulation box ===
  integer::nbox=1,nbxmax& !< maximum number of boxes
   ,npabmax& !< maximum number of box pairs (for swatch and swap)
   ,boxlink
  real::beta,temp=-1.0_dp,rintramax=0.0_dp
  logical::licell=.false.
  real,allocatable,target::boxlx(:),boxly(:),boxlz(:),rcut(:),rcutnn(:),kalp(:)
  logical,allocatable,target::lsolid(:),lrect(:)
  real,allocatable::express(:),zshift(:),dshift(:)
  integer,allocatable::numberDimensionIsIsotropic(:),ininch(:,:),inix(:),iniy(:),iniz(:),inirot(:),inimix(:)
  integer(long_int),allocatable::ghost_particles(:) ! compatibility with very large vapor boxes
  logical,allocatable::lideal(:) !< if lideal=.true. then intermolecular interactions are not computed

  !=== Information about molecule types ===
  integer::nmolty=1,ntmax& !< maximum number of types of chains
   ,npamax& !< maximum number of pairs to switch or swatch
   ,numax !< maximum number of units
  integer,allocatable::nunit(:),isolute(:),rindex(:),riutry(:,:)
  logical,allocatable::lelect(:),lrigid(:),lq14scale(:),lbranch(:)
  real,allocatable::eta2(:,:),qscale(:),B(:) !< chemical potential
  integer,allocatable::ntype(:,:),leaderq(:,:)&
   ,invib(:,:),itvib(:,:,:),ijvib(:,:,:)&
   ,inben(:,:),itben(:,:,:),ijben2(:,:,:),ijben3(:,:,:)&
   ,intor(:,:),ittor(:,:,:),ijtor2(:,:,:),ijtor3(:,:,:),ijtor4(:,:,:)
  real,allocatable::rigid_intra_dist(:,:)
  character(LEN=10),allocatable::molecname(:)

  !=== tolerance values ===
  real::distance_tolerance=1E-6_dp

  !=== Information about particles in the system ===
  integer::nchain=4,nmax !< maximum number of chains + 2
  integer,allocatable::moltyp(:),nboxi(:)
  real,allocatable::masst(:),rcmu(:),xcm(:),ycm(:),zcm(:)& !< center-of-mass coordinates of each chain
   ,sxcm(:),sycm(:),szcm(:)& !< center-of-mass coordinates of each chain in scaled units, for use in non-orthorhombic simulation cells
   ,rxu(:,:),ryu(:,:),rzu(:,:),qqu(:,:)&
   ,rxu_update(:),ryu_update(:),rzu_update(:)

  !=== Counters ===
  integer,allocatable::ncmt(:,:)& !< (ibox,itype): number of molecules of itype in ibox
   ,ncmt2(:,:,:)& !< (ibox,itype,itype2): number of molecules of itype in ibox, in stage itype2 as in expanded ensemble
   ,nchbox(:)& !< number of molecules (of any type) in each box
   ,parall(:,:)& !< (itype,j): index of the j-th molecule of itype
   ,parbox(:,:,:)& !< (j,ibox,itype): index of j-th molecule of itype in ibox
   ,temtyp(:) !< number of molecules of each molecule type

  !=== Information about Monte Carlo moves ===
  integer::tmcc,iratio=500
  real::rmin=1.2_dp,softcut=100.0_dp,softlog

  !*** Volume moves ***
  real::pmvol,pmvolx=0.333_dp,pmvoly=0.666_dp,tavol=0.5_dp
  real,allocatable::pmvlmt(:),pmvolb(:),rmvol(:),rmhmat(:,:)
  integer::nvolb,iratv=500
  integer,allocatable::box5(:),box6(:)
  logical::l_bilayer
  real::pm_consv, pmvol_xy

  !*** CBMC particle identity switch (swatch) moves ***
  real::pmswat
  integer::nswaty
  real,allocatable::pmsatc(:),pmswtcb(:,:)
  integer,allocatable::nswatb(:,:),nsampos(:)& !< number of beads that remain in the same position
   ,ncut(:,:),splist(:,:,:),gswatc(:,:,:),nswtcb(:),box3(:,:),box4(:,:)
  integer,allocatable::ncutsafe(:,:), gswatcsafe(:,:,:) !< Paul -- safe-swatch variables
  !> \bug liswatch, other, liswinc not initialized
  ! logical,allocatable::liswinc(:,:)
  ! logical::liswatch !< prevents non-grown beads from being included in the new growth in boltz
  ! integer::other

  !*** CBMC particle transfer (swap) moves ***
  real::pmswap
  real,allocatable::pmswmt(:),pmswapb(:,:)
  integer,allocatable::nswapb(:),box1(:,:),box2(:,:)

  !** AVBMC moves **
  logical,allocatable::lavbmc1(:),lavbmc2(:),lavbmc3(:),lbias(:)
  real,allocatable::pmbias(:),pmbsmt(:),pmbias2(:),favor(:),favor2(:)
  real::rbsmax=3.5_dp,rbsmin=2.5_dp,vol_eff

  !*** Regular CBMC moves ***
  real::pmcb
  real,allocatable::pmcbmt(:),pmall(:)

  !** SAFE CBMC **
  real,allocatable::pmfix(:)
  integer,allocatable::iring(:),nrig(:),irig(:,:)& !< the site rigid sites will be grown from
   ,frig(:,:)& !< the previous site (not kept rigid)
   ,nrigmin(:),nrigmax(:) !< the minimum and maximum amounts of the chain to keep rigid
  logical,allocatable::lring(:),lrig(:),lrigi(:,:)
  integer,parameter::maxbin=201 !< SAFECBMC max number of bins
  integer::iupdatefix=100
  logical::lpresim=.false.

  !*** CBMC_bend_table (Bin's tabulated angle CBMC method)
  logical::L_cbmc_bend=.false.  ! whether to use tabulated bending table for CBMC move
  real,allocatable::lin_bend_type(:),lin_bend_table(:,:),lin_bend_prob(:,:) ! Tabulated bending table for linear CBMC bead growth
  real,allocatable::br_bend_type(:),br_bend_theta1(:,:),br_bend_theta2(:,:,:),br_bend_phi12(:,:,:,:),br_bend_prob(:,:) ! Tabulated bending table for one-branch CBMC bead growth
  integer,allocatable::lin_bend_dim(:),br_bend_dim1(:),br_bend_dim2(:),br_bend_dim3(:) ! The dimensions for the above bending tables

  !*** Group_CBMC
  real,allocatable::pmgroup(:) !< probability of performing group-CBMC
  integer::gcbmc_box_num=0 !< box number that contains the repeat units
  integer::gcbmc_mol_num !< number of molecules using group CBMC
  integer,allocatable::gcbmc_mol_list(:) !< corresponding list of moltype with gcbmc_moltype
  integer,allocatable::gcbmc_unit_num(:) !< number of repeat unit for moltype with gcbmc_moltype
  integer,allocatable::gcbmc_unit_moltype(:,:) !< molecule type of repeat units
  integer,allocatable::gcbmc_unit_list(:,:,:) !< corresponding list of beads
  logical::l_gcbmc_movie=.false. !< whether to write gcbmc reservoir box information to the movie file

  !*** CBMC shared variables ***
  integer,allocatable::nugrow(:),nmaxcbmc(:),iurot(:),maxgrow(:)&
   ,nchoi1(:),nchoi(:),nchoir(:),nchoih(:),nchtor(:),nchoig(:),nchbna(:),nchbnb(:)& !< number of candidates during CBMC regrowth for the first bead, subsequent beads (flexible and rigid), explicit-hydrogen, torsion, and bendings
   ,nrotbd(:),irotbd(:,:),icbdir(:),icbsta(:)&
   ,growfrom(:)& !< (index): the bead from which the new beads are to be grown at the index-th step
   ,growprev(:)& !< (index): the bead that exists and is connected to growfrom(index)
   ,grownum(:)& !< (index): the number of beads to be grown from growfrom(index)
   ,growlist(:,:) !< (index,count): the count-th one of the new beads to be grown from growfrom(index), grownum(index) entries
  real,allocatable::pmrotbd(:,:),rxp(:,:),ryp(:,:),rzp(:,:)& !< (iunit,itrial) coordinates of iunit of the selected molecule for the itrial-th candidate configuration, for use in rosenbluth
   ,bfac(:),vtr(:,:),vtrorient(:),vtrelect_intra(:),vtrelect_inter(:)& !< Boltzmann factors and energies for the configuration in r{x,y,z}p
   ,rxuion(:,:),ryuion(:,:),rzuion(:,:),qquion(:,:)& !< (iunit,flag): coordinates of iunit of the selected molecule in its old (flag=1) or new (flag=2) state, for subroutine energy (one-particle energy)
   ,rxnew(:),rynew(:),rznew(:)& !< coordinates of iunit of the selected molecule in its new configuration, for use in config, swap and swatch
   ,xvec(:,:),yvec(:,:),zvec(:,:)
  logical,allocatable::lovr(:),lexist(:),lplace(:,:)
  integer::nchmax& !< maximum number of choices for trial sites in CBMC growth
   ,nchtor_max& !< maximum number of choices for torsion CBMC growth
   ,nchbn_max& !< maxium number of choices for bending CBMC growth
   ,moltion(2)
  real::rcutin=6.0_dp,weight,weiold,vold(nEnergy),vnew(nEnergy),vneworient,voldorient


  !*** Fluctuating charge moves ***
  real,allocatable::xiq(:),jayself(:),jayq(:)
  logical,parameter::lfepsi=.false. !< If lfepsi is true, the fluctuation of epsilon is used instead of the fluctuation of the sigma.
  logical,allocatable::lflucq(:),lqtrans(:)
  real::pmflcq,taflcq=0.95_dp,fqbeta
  real,allocatable::bnflcq(:,:),bsflcq(:,:),bnflcq2(:,:),bsflcq2(:,:)& !< accumulators for statistics of fluctuating charge moves;
   !< bnflcq and bsflcq record the number of all and successful attempts between updates of maximum displacements, while bnflcq2
   !< and bsflcq record the values for the entire run
   ,fqegp(:),pmfqmt(:),rmflcq(:,:)
  integer,allocatable::nchoiq(:)
  integer::nswapq

  !*** Expanded ensemble moves ***
  integer,parameter::smax=35 !< max no. of mstates
  logical::leemove,lmstate,leeacc
  integer::fmstate,sstate1,sstate2,nstate,box_state(smax),eepointp,eeirem,boxrem1,boxins1,ee_prob(smax)
  real::wee_ratio,psi(smax),um_markov(smax,smax),eeratio
  integer::mstate,ee_moltyp(smax),nmolty1
  logical,allocatable::lexpand(:)
  real::pmexpc
  real,allocatable::pmeemt(:),sigma_f(:,:),epsilon_f(:,:),ee_qqu(:,:),rminee(:)
  integer,allocatable::rmexpc(:),eetype(:)

  !*** New expanded ensemble moves ***
  real::pmexpc1

  !*** Atom translation ***
  real::pm_atom_tra,Armtrax,Armtray,Armtraz

  !*** Specific atom translation ***
  integer natomtrans_atoms
  integer,allocatable::atomtrans_atomlst(:),atomtrans_moleclst(:)

  !*** Translation ***
  real::pmtra,tatra=0.5_dp
  real,allocatable::pmtrmt(:),rmtrax(:,:),rmtray(:,:),rmtraz(:,:)

  !*** Rotation ***
  real::tarot=0.5_dp
  real,allocatable::pmromt(:),rmrotx(:,:),rmroty(:,:),rmrotz(:,:)

  !=== MPI-related ===
  integer::myid=0,numprocs=1,groupid=0
  integer,parameter::rootid=0

  !=== OPENMP-related ===
  integer::thread_id=0,thread_num=1,thread_num_max=1,thread_num_proc=1

  !** Parameters for 2nd virial coefficient calculation ***
  ! real,parameter::a0 = 0.2003E0_dp,b0 = 1.3946E0_dp,aslope = 8.85E5_dp,bslope = 158.25E0_dp,ashift = 7.227E5_dp,bshift = 505.97E0_dp ! slope=0.3 a=3.05
  ! real,parameter::a0 = 0.15561E0_dp,b0 = 1.2960E0_dp,aslope = 5.5475E5_dp,bslope = 131.23E0_dp,ashift = 4.121E5_dp,bshift = 381.96E0_dp ! slope=0.3 a=2.90
  ! real,parameter::a0 = 0.16833E0_dp,b0 = 1.3299E0_dp,aslope = 6.5125E5_dp,bslope = 139.75E0_dp,ashift = 4.9975E5_dp,bshift = 419.81E0_dp ! slope=0.3 a=2.95
  ! real,parameter::a0 = 0.13150E0_dp,b0 = 1.2574E0_dp,aslope = 4.2863E5_dp,bslope = 117.5E0_dp,ashift = 3.0208E5_dp,bshift = 323.6E0_dp ! slope=0.3 a=2.82
  ! real,parameter::a0 = 0.11037E0_dp,b0 = 1.1942E0_dp,aslope = 3.4013E5_dp,bslope = 108E0_dp,ashift = 2.2868E5_dp,bshift = 285.04E0_dp ! slope=0.3 a=2.75
  ! real,parameter::a0 = 0.14035E0_dp,b0 = 1.2596E0_dp,aslope = 4.7263E5_dp,bslope = 123.25E0_dp,ashift = 3.3978E5_dp,bshift = 347.63E0_dp ! slope=0.3 a=2.85
  ! real,parameter::a0 = 0.14820E0_dp,b0 = 1.2689E0_dp,aslope = 5.2125E5_dp,bslope = 128.75E0_dp,ashift = 3.8220E5_dp,bshift = 371.29E0_dp ! slope=0.3 a=2.88
  real,parameter::a0 = -0.22818E0_dp,b0 = 0.44662E0_dp,aslope = 14.0738E5_dp,bslope = 351.25E0_dp,ashift = 3.5852E5_dp,bshift = 358.94E0_dp ! slope=0.5 a=2.71
  ! real,parameter::a0 = -0.23379E0_dp,b0 = 0.43630E0_dp,aslope = 12.78125E5_dp,bslope = 337.5E0_dp,ashift = 3.1904E5_dp,bshift = 337.86E0_dp ! slope=0.5 a=2.68
  ! real,parameter::a0 = -0.23808E0_dp,b0 = 0.425507E0_dp,aslope = 11.7775E5_dp,bslope = 322.04E0_dp,ashift = 2.8902E5_dp,bshift = 326.875E0_dp ! slope=0.5 a=2.655
  ! real,parameter::a0 = -0.23641E0_dp,b0 = 0.43409E0_dp,aslope = 12.171E5_dp,bslope = 330E0_dp,ashift = 3.0075E5_dp,bshift = 326.52E0_dp ! slope=0.5 a=2.665
  ! real,parameter::a0 = -0.23903E0_dp,b0 = 0.42231E0_dp,aslope = 11.5875E5_dp,bslope = 325E0_dp,ashift = 2.8339E5_dp,bshift = 319.44E0_dp ! slope=0.5 a=2.65
  ! real,parameter::a0 = -0.24796E0_dp,b0 = 0.40345E0_dp,aslope = 9.8238E5_dp,bslope = 304E0_dp,ashift = 2.3206E5_dp,bshift = 288.72E0_dp ! slope=0.5 a=2.60
  ! real,parameter::a0 = -0.25703E0_dp,b0 = 0.38299E0_dp,aslope = 8.3075E5_dp,bslope = 284.38E0_dp,ashift = 1.8945E5_dp,bshift = 261.09E0_dp ! slope=0.5 a=2.55
  ! real,parameter::a0 = -0.39106E0_dp,b0 = 0.08431E0_dp,aslope = 25.22875E5_dp,bslope = 662.75E0_dp,ashift = 3.0689E5_dp,bshift = 335.43E0_dp ! slope = 0.7 a=2.52
  ! real,parameter::a0 = -0.39233E0_dp,b0 = 0.08203E0_dp,aslope = 24.425E5_dp,bslope = 653.75E0_dp,ashift = 2.9499E5_dp,bshift = 328.60E0_dp ! slope = 0.7 a=2.51
  ! real,parameter::a0 = -0.39357E0_dp,b0 = 0.07946E0_dp,aslope = 23.641E5_dp,bslope = 645E0_dp,ashift = 2.835E5_dp,bshift = 322.13E0_dp ! slope = 0.7 a=2.50
  ! real,parameter::a0 = -0.394906E0_dp,b0 = 0.075678E0_dp,aslope = 22.88625E5_dp,bslope = 637.735E0_dp,ashift = 2.72469E5_dp,bshift = 316.22E0_dp ! slope = 0.7 a=2.49
  ! real,parameter::a0 = -0.396175E0_dp,b0 = 0.072973E0_dp,aslope = 22.149125E5_dp,bslope = 629E0_dp,ashift = 2.61798E5_dp,bshift = 309.95E0_dp ! slope = 0.7 a=2.48
  ! real,parameter::a0 = -0.39745E0_dp,b0 = 0.070157E0_dp,aslope = 21.4335E5_dp,bslope = 620.75E0_dp,ashift = 2.5151E5_dp,bshift = 306.89E0_dp ! slope = 0.7 a=2.47
  ! real,parameter::a0 = -0.399983E0_dp,b0 = 0.064159E0_dp,aslope = 20.06425E5_dp,bslope = 604.75E0_dp,ashift = 2.32037E5_dp,bshift = 292.09E0_dp ! slope = 0.7 a=2.45
  ! real,parameter::a0 = -0.40629E0_dp,b0 = 0.050088E0_dp,aslope = 16.9775E5_dp,bslope = 565.50E0_dp,ashift = 1.89213E5_dp,bshift = 263.93E0_dp ! slope = 0.7 a=2.40
  ! real,parameter::a0 = -0.48332E0_dp,b0 = -0.12334E0_dp,aslope = 34.64125E5_dp,bslope = 1014.25E0_dp,ashift = 2.2815E5_dp,bshift = 294.26E0_dp ! slope = 0.9 a=2.30
  ! real,parameter::a0 = -0.48137E0_dp,b0 = -0.11887E0_dp,aslope = 36.99125E5_dp,bslope = 1041.25E0_dp,ashift = 2.4744E5_dp,bshift = 306.23E0_dp ! slope = 0.9 a=2.32
  ! real,parameter::a0 = 0.0E0_dp,b0 = 0.0E0_dp,aslope = 3.0E5_dp,bslope = 0.0E0_dp,ashift = 8.0E5_dp,bshift = 1200.0E0_dp ! slope = 0.3 a=2.85

  type(LookupTable)::atoms

  ! kdtree related variables
  ! Data types
  type :: interval
      real :: lower, upper
  end type interval

  type :: tree_node
      ! internal tree node
      type(tree_node), pointer :: left_node !< left node
      type(tree_node), pointer :: right_node !< right node
      type(tree_node), pointer :: parent_node !< parent node
      real, dimension(3) :: coord !< coordinates
      integer :: ichain !< chain no.
      integer :: ibead  !< bead no.
      integer :: ix     !< with a value of 0, -1, 1, indicating which periodic image the particle is
      integer :: iy
      integer :: iz
      type(interval), pointer :: cube !< the bounding box of the node
      logical :: l_cube_updated !< if the cube has been updated
      integer :: cut_dim !< the cutting dimension of the node
      integer :: height !< height of the current node
  end type tree_node

  type :: tree
     type(tree_node), pointer :: tree_root !< root of this tree
     integer :: height !< height of the tree
     integer :: node_num   !< number of nodes in this tree
     type(interval), pointer :: cube(:) !< the bounding cube of the tree
     type(interval), pointer :: bound(:) !< the min and max of the coordinates in the center box
     type(interval), pointer :: bound_all(:) !< the min and max of the coordinates in the periodic images
     integer :: box     !< which box the tree represents
  end type tree

  type :: tree_ptr
      type(tree), pointer :: tree
  end type tree_ptr

  type(tree_ptr), allocatable :: mol_tree(:)
  logical :: lkdtree=.false.             !< whether kdtree is used in general
  logical, allocatable :: lkdtree_box(:) !< whether this box uses kdtree
  real, allocatable :: kdtree_buffer_len(:)  !< the length of the buffer region used in kdtree
  integer, allocatable :: tree_height(:)  !< the optimal height of each tree

CONTAINS
  subroutine allocate_system()
    integer,parameter::initial_size=15
    integer::jerr

    if (allocated(boxlx)) deallocate(boxlx,boxly,boxlz,rcut,rcutnn,kalp,lsolid,lrect,express,ghost_particles&
        ,numberDimensionIsIsotropic,inix,iniy,iniz,inirot,inimix,nchoiq,box5,box6,zshift,dshift,rmvol,pmvlmt&
        ,pmvolb,lideal,ltwice,rmhmat,dipolex,dipoley,dipolez,nchbox,vbox,stat=jerr)
    allocate(boxlx(nbxmax),boxly(nbxmax),boxlz(nbxmax),rcut(nbxmax),rcutnn(nbxmax),kalp(nbxmax),lsolid(nbxmax),lrect(nbxmax)&
     ,express(nbxmax),ghost_particles(nbxmax),numberDimensionIsIsotropic(nbxmax),inix(nbxmax),iniy(nbxmax)&
     ,iniz(nbxmax),inirot(nbxmax),inimix(nbxmax),nchoiq(nbxmax),box5(npabmax),box6(npabmax),zshift(nbxmax),dshift(nbxmax)&
     ,rmvol(nbxmax),pmvlmt(nbxmax),pmvolb(npabmax),lideal(nbxmax),ltwice(nbxmax),rmhmat(nbxmax,9),dipolex(nbxmax)&
     ,dipoley(nbxmax),dipolez(nbxmax),nchbox(nbxmax),vbox(nEnergy,nbxmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'allocate_system.1: allocation failed',jerr)

    if (allocated(nrotbd)) deallocate(nrotbd,xcm,ycm,zcm,pmsatc,pmswtcb,nswatb,nsampos,ncut,gswatc,nswtcb,box3,box4,ncutsafe,gswatcsafe,temtyp,B,molecname,nunit,nugrow,nmaxcbmc,iurot,maxgrow,isolute,iring,nrig,irig,frig,nrigmin,nrigmax,rindex,riutry,lelect,lflucq,lqtrans,lexpand,lavbmc1,lavbmc2,lavbmc3,lbias,lring,lrigid,lrig,lq14scale,fqegp,eta2,qscale,pmbias,pmbsmt,pmbias2,rmtrax,rmtray,rmtraz,rmrotx,rmroty,rmrotz,lbranch,ininch,rmflcq,pmswmt,pmswapb,pmcbmt,pmall,pmfix,pmfqmt,pmeemt,pmtrmt,pmromt,pmgroup,nswapb,box1,box2,nchoi1,nchoi,nchoir,nchoih,nchoig,nchtor,nchbna,nchbnb,icbdir,icbsta,lrplc,masst,rmexpc,eetype,ncmt,ncmt2,parall,parbox,bnflcq,bsflcq,bnflcq2,bsflcq2,rxwell,rywell,rzwell,sxwell,sywell,szwell,nwell,lwell,moltyp,rcmu,sxcm,sycm,szcm,nboxi,favor,favor2,ntype,leaderq,invib,itvib,ijvib,inben,itben,ijben2,ijben3,intor,ittor,ijtor2,ijtor3,ijtor4,irotbd,pmrotbd,rigid_intra_dist,stat=jerr)
    allocate(nrotbd(ntmax),xcm(nmax),ycm(nmax),zcm(nmax),pmsatc(npamax),pmswtcb(npamax,npabmax),nswatb(npamax,2)&
     ,nsampos(npamax),ncut(npamax,2),gswatc(npamax,2,2*npamax),nswtcb(npamax),box3(npamax,npabmax),box4(npamax,npabmax)&
     ,ncutsafe(npamax,2),gswatcsafe(npamax,3,2*npamax)&
     ,temtyp(ntmax),B(ntmax),molecname(ntmax),nunit(ntmax),nugrow(ntmax),nmaxcbmc(ntmax),iurot(ntmax),maxgrow(ntmax),isolute(ntmax),iring(ntmax)&
     ,nrig(ntmax),irig(ntmax,6),frig(ntmax,6),nrigmin(ntmax),nrigmax(ntmax),rindex(ntmax),riutry(ntmax,initial_size),lelect(ntmax)&
     ,lflucq(ntmax),lqtrans(ntmax),lexpand(ntmax),lavbmc1(ntmax),lavbmc2(ntmax)&
     ,lavbmc3(ntmax),lbias(ntmax),lring(ntmax),lrigid(ntmax),lrig(ntmax),lq14scale(ntmax),fqegp(ntmax),eta2(nbxmax,ntmax)&
     ,qscale(ntmax),pmbias(ntmax),pmbsmt(ntmax),pmbias2(ntmax),rmtrax(ntmax,nbxmax),rmtray(ntmax,nbxmax),rmtraz(ntmax,nbxmax)&
     ,rmrotx(ntmax,nbxmax),rmroty(ntmax,nbxmax),rmrotz(ntmax,nbxmax),lbranch(ntmax),ininch(ntmax,nbxmax),rmflcq(ntmax,nbxmax)&
     ,pmswmt(ntmax),pmswapb(ntmax,npabmax),pmcbmt(ntmax),pmall(ntmax),pmfix(ntmax),pmfqmt(ntmax),pmeemt(ntmax),pmtrmt(ntmax)&
     ,pmromt(ntmax),pmgroup(ntmax),nswapb(ntmax),box1(ntmax,npabmax),box2(ntmax,npabmax),nchoig(ntmax),nchoi1(ntmax),nchoi(ntmax),nchoir(ntmax)&
     ,nchoih(ntmax),nchtor(ntmax),nchbna(ntmax),nchbnb(ntmax),icbdir(ntmax),icbsta(ntmax),lrplc(ntmax),masst(ntmax)&
     ,rmexpc(ntmax),eetype(ntmax),ncmt(nbxmax,ntmax),ncmt2(nbxmax,ntmax,20),parall(ntmax,nmax)&
     ,parbox(nmax,nbxmax,ntmax),bnflcq(ntmax,nbxmax),bsflcq(ntmax,nbxmax),bnflcq2(ntmax,nbxmax),bsflcq2(ntmax,nbxmax)&
     ,rxwell(nw,ntmax),rywell(nw,ntmax),rzwell(nw,ntmax),sxwell(nw,ntmax),sywell(nw,ntmax),szwell(nw,ntmax),nwell(ntmax)&
     ,lwell(ntmax),moltyp(nmax),rcmu(nmax),sxcm(nmax),sycm(nmax),szcm(nmax),nboxi(nmax)&
     ,favor(nmax),favor2(nmax),ntype(ntmax,initial_size),leaderq(ntmax,initial_size)&
     ,invib(ntmax,initial_size),itvib(ntmax,initial_size,nvib_max),ijvib(ntmax,initial_size,nvib_max)&
     ,inben(ntmax,initial_size),itben(ntmax,initial_size,nben_max),ijben2(ntmax,initial_size,nben_max)&
     ,ijben3(ntmax,initial_size,nben_max),intor(ntmax,initial_size),ittor(ntmax,initial_size,ntor_max)&
     ,ijtor2(ntmax,initial_size,ntor_max),ijtor3(ntmax,initial_size,ntor_max),ijtor4(ntmax,initial_size,ntor_max)&
     ,irotbd(initial_size,ntmax),pmrotbd(initial_size,ntmax),rigid_intra_dist(max_num_rigid_beads,ntmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'allocate_system.2: allocation failed',jerr)
  end subroutine allocate_system

  subroutine allocate_molecule()
    integer::jerr
    if (allocated(splist)) deallocate(splist,lexist,lexclu,a15type,epsilon_f,sigma_f,ljscale,qscale2,ee_qqu&
     ,rxnew,rynew,rznew,rxu,ryu,rzu,qqu,rxuion,ryuion,rzuion,qquion,rxu_update,ryu_update,rzu_update,linclu&
     ,lqinclu,lainclu,xvec,yvec,zvec,growfrom,growprev,grownum,growlist,awell,lexclu_zeo,stat=jerr)
    allocate(splist(npamax,numax,2),lexist(numax),lexclu(ntmax,numax,ntmax,numax),a15type(ntmax,numax,numax)&
     ,epsilon_f(2,numax),sigma_f(2,numax),ljscale(ntmax,numax,numax),qscale2(ntmax,numax,numax),ee_qqu(numax,smax)&
     ,rxnew(numax),rynew(numax),rznew(numax),rxu(nmax,numax),ryu(nmax,numax),rzu(nmax,numax),qqu(nmax,numax)&
     ,rxuion(numax,2),ryuion(numax,2),rzuion(numax,2),qquion(numax,2),rxu_update(numax),ryu_update(numax),rzu_update(numax)&
     ,linclu(ntmax,numax,numax),lqinclu(ntmax,numax,numax),lainclu(ntmax,numax,numax)&
     ,xvec(numax,numax),yvec(numax,numax),zvec(numax,numax),growfrom(numax)&
     ,growprev(numax),grownum(numax),growlist(numax,numax),awell(numax,numax,ntmax)&
     ,lexclu_zeo(ntmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_molecule: allocation failed',jerr)
    end if
  end subroutine allocate_molecule

  subroutine checkAtom()
    if (.not.allocated(atoms%list)) call err_exit(__FILE__,__LINE__,": ATOMS section has not been defined!",myid+1)
  end subroutine checkAtom

  subroutine setup_mpi(nrank,tid,gid)
    integer,intent(in)::nrank,tid,gid

    numprocs=nrank
    myid=tid
    groupid=gid

  end subroutine setup_mpi
end module sim_system
