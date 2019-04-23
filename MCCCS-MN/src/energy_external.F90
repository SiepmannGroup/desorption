MODULE energy_external
  use var_type,only:dp
  use const_math,only:onepi,twopi,fourpi
  use util_math,only:mbessel
  use sim_system,only:vvdW,qqu
  implicit none
  private
  save
  public::U_ext,init_energy_external

  integer::surface_type=1& !< 1: Steele 10-4-3 potential, with type 1 (Lennard-Jones 12-6) for input parameters
   !< 2: Steele 9-3 potential, with type 1 (Lennard-Jones 12-6) for input parameters
   !< 3: Hautman 12-3 potential, with type 7 (Lennard-Jones 12-6-8) nonbond interactions for input parameters
   ,ntsubst=190 !< bead type of substrate, used to input potential parameters
  real::rsol=0.114E0_dp& !< unit 1/A^3
   ,delta=3.40E0_dp& !< delta and a1 are in Angstroms
   ,a1=2.460E0_dp
  logical::double_surface=.true.
  real,allocatable::Elect_field(:)

contains
!> \brief Calculates the energy of a bead with a featureless, planar slit surface
!> \par Literature
!> 1. William A Steele, "The interaction of gases with solid surfaces", Pergamon Press, 1974
!> 2. Joseph Hautman and Michael L Klein, "Simulation of a monolayer of alkyl thiol chains", J Chem Phys 1989,91:4994
  function slitpore(z,ntij)
    real::slitpore
    real,intent(in)::z
    integer,intent(in)::ntij
    real::sig3,dz3,dz12

    ! Note that vvdW(1,ntij) is 4pi

    if (surface_type.eq.1) then
       ! Steele 10-4-3 potential
       slitpore = onepi/2.0_dp*rsol*vvdW(1,ntij)*vvdW(3,ntij)*delta*(0.4_dp*(vvdW(2,ntij)/z)**10-(vvdW(2,ntij)/z)**4&
        -vvdW(3,ntij)**2/(3.0_dp*delta*(z+0.61_dp*delta)**3))
    else if (surface_type.eq.2) then
       ! Steele 9-3 potential
       sig3 = vvdW(2,ntij)**3
       dz3 = sig3/(z**3)
       slitpore = onepi/6.0_dp*rsol*vvdW(1,ntij)*sig3*dz3*(2.0_dp/15.0_dp*dz3*dz3-1.0_dp)
    else if (surface_type.eq.3) then
       ! Hautman 12-3 potential
       dz3 = (z-vvdW(3,ntij))**3
       dz12 = dz3**3
       slitpore = vvdW(1,ntij)/dz12 - vvdW(2,ntij)/dz3
    end if
  end function slitpore

!> \brief Calculates the energy of a bead with a graphite surface with x,y-dependent potentials
  function exgrph(x,y,z,ntij)
    real::exgrph
    real,intent(in)::x,y,z
    integer,intent(in)::ntij
    real::aa2,a1sq,e1,fxy,bb,cc,dd,k2,k5,zzz

    exgrph = slitpore(z,ntij)

    e1 = 0.0_dp
    fxy = 0.0_dp
    a1sq = a1**2
    aa2 = (vvdW(3,ntij)/a1sq)**3
    bb = aa2*onepi*vvdW(1,ntij)/sqrt(3.0_dp)
    ! bb = onepi*vvdW(1,ntij)*vvdW(3,ntij)**3/(sqrt(3.0_dp)*a1**6)
    cc = aa2/(30.0_dp*(twopi/sqrt(3.0_dp))**5)
    ! cc = vvdW(3,ntij)**6/(30.0_dp*a1**6*(twopi/sqrt(3.0_dp)**5))
    dd = 2.0_dp*(twopi/sqrt(3.0_dp))**2
    zzz = fourpi*z/(sqrt(3.0_dp)*a1)
    k2 = mbessel(zzz,2.0_dp)
    k5 = mbessel(zzz,5.0_dp)
    e1 = bb*(cc * k5 * (a1/z)**5 - dd * k2 * (a1/z)**2)
    fxy = -2.0_dp*(cos(twopi*(x/a1 + y/sqrt(3.0_dp)/a1)) + cos(twopi*(x/a1 - y/sqrt(3.0_dp)/a1)) + cos(fourpi*y/sqrt(3.0_dp)/a1))
    exgrph = exgrph + e1*fxy
  end function exgrph

  function U_ext(ibox,i,j,ntj)
    use const_phys,only:eXV_to_K
    use sim_system,only:rxu,ryu,rzu,nntype,lelect_field,lexzeo,lslit,lgraphite,lsami,lmuir,io_output,nchain,boxlz
    use energy_sami,only:exsami,exmuir
    use zeolite
    real::U_ext
    integer,intent(in)::ibox,i,j,ntj
    integer::ntij
    real::vtmp

    U_ext=0.0_dp

    if (lelect_field) then
       !> Calculates interaction of molecule i with an external field E
       !> \par Units
       !> E in V/A, q in e, rz in A \n
       !> E*q*rz = V*e \n
       !> 1 V*e = 11600 K
       !> \author 06/24/07 by KM
       U_ext = U_ext - Elect_field(ibox)*rzu(i,j)*qqu(i,j)*eXV_to_K
    end if

    if (ibox.ne.1) return

    if (lexzeo) then
       vtmp=exzeo(rxu(i,j),ryu(i,j),rzu(i,j),ntj,ignoreTable=.false.)
       if (i.le.nchain.and.abs(vtmp).gt.upper_limit_zeo) then
          write(io_output,*) '###problem: energy for molnum = ',i,'bead =',j,'is larger than'
          write(io_output,*) 'upper_limit_zeo. You will want to increase upperLimit in topmon.inp and/or'
          write(io_output,*) 'look into this. rx, ry, rz =',rxu(i,j),ryu(i,j),rzu(i,j),'vtmp = ',vtmp
       end if
       U_ext = U_ext + vtmp
    end if

    if (lslit) then
       ! carbon slitpore
       ntij = (ntj-1)*nntype + ntsubst
       ! calculate interaction with surfaces at the bottom and the top of the box
       U_ext = U_ext + slitpore(rzu(i,j),ntij)
       if (double_surface) U_ext = U_ext + slitpore(boxlz(ibox)-rzu(i,j),ntij)
    end if

    if (lgraphite) then
       ntij = (ntj-1)*nntype + ntsubst
       U_ext = U_ext + exgrph(rxu(i,j),ryu(i,j),rzu(i,j),ntij)
    end if

    if (lsami)  U_ext = U_ext + exsami(rzu(i,j),ntj)
    if (lmuir)  U_ext = U_ext + exmuir(rzu(i,j),ntj)

  end function U_ext

  subroutine init_energy_external(io_input,lprint)
    use util_runtime,only:err_exit
    use util_string,only:format_n
    use util_mp,only:mp_bcast
    use sim_system,only:nbox,nbxmax,io_output,myid,rootid,groupid
    INTEGER,INTENT(IN)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::jerr
    namelist /external_field/ surface_type,ntsubst,rsol,delta,a1,Elect_field,double_surface

    !> read namelist external_field
    if (allocated(Elect_field)) deallocate(Elect_field,stat=jerr)
    allocate(Elect_field(nbxmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_energy_external: allocation failed',jerr)
    Elect_field=0.0_dp

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=external_field,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: external_field',jerr)
    end if

    call mp_bcast(surface_type,1,rootid,groupid)
    call mp_bcast(ntsubst,1,rootid,groupid)
    call mp_bcast(rsol,1,rootid,groupid)
    call mp_bcast(delta,1,rootid,groupid)
    call mp_bcast(a1,1,rootid,groupid)
    call mp_bcast(Elect_field,nbox,rootid,groupid)
    call mp_bcast(double_surface,1,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST EXTERNAL_FIELD','------------------------------------------'
       if (surface_type.eq.1) then
          write(io_output,'(A)') 'Steele 10-4-3 slit pore'
       else if (surface_type.eq.2) then
          write(io_output,'(A)') 'Steele 9-3 slit pore'
       else if (surface_type.eq.3) then
          write(io_output,'(A)') 'Hautman 12-3 monolayer'
       end if
       if (double_surface) then
          write(io_output,'(A)') 'On both sides of the simulation box (z = 0 & z = boxlz)'
       else
          write(io_output,'(A)') 'On one side of the simulation box (z = 0)'
       end if
       write(io_output,'(A,I0)') 'Surface material bead type: ',ntsubst
       write(io_output,'(A,F8.5,A)') 'Surface atom density: ',rsol,' [Ang^-3]'
       write(io_output,'(A,F8.5,A)') 'Surface layer spacing: ',delta,' [Ang]'
       write(io_output,'(A,F8.5,A)') 'a1: ',a1,' [Ang]'
       write(io_output,'(A,'//format_n(nbox,'(3X,G16.9)')//',A)') 'Electric field in z direction:',Elect_field(1:nbox),' [V/A]'
    end if
  end subroutine init_energy_external
end MODULE energy_external
