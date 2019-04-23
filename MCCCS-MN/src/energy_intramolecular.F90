MODULE energy_intramolecular
  use var_type,only:dp,default_string_length
  use const_math,only:onepi,twopi,raddeg
  use util_runtime,only:err_exit
  use util_memory,only:reallocate
  use util_string,only:uppercase,integer_to_string,stripComments
  use util_files,only:readLine
  use util_search,only:LookupTable,initiateTable,destroyTable,addToTable
  use sim_system,only:io_output,brvib,brvibk,brben,brbenk,maxRegrowVib,minRegrowVib,myid,L_spline,L_linear
  implicit none
  private
  save
  public::read_ff_bonded,vtorso,lininter_vib,lininter_bend,U_torsion,U_bonded,allocate_energy_bonded,bonds,angles,dihedrals

  integer,parameter::torsion_nParameter(-1:8)=(/0,0,4,10,3,2,5,10,4,5/)
  integer,allocatable::vib_type(:),ben_type(:),torsion_type(:) !< type -1: tabulated potential
  !< type 0: dummy torsion type for setting up interaction table
  !< type 1: OPLS potential (three terms)
  !< V(phi) = vtt0 + vtt1*[1+cos(phi)] + vtt2*[1-cos(2*phi)] + vtt3*[1+cos(3*phi)]
  !< type 2: Ryckaert-Bellemans potential
  !< V(psi) = vtt0 + vtt1*cos(psi) + vtt2*cos(psi)^2 + ... + vtt9*cos(psi)^9
  !< type 3: periodic type
  !< V(phi) = vtt0*[1+cos(vtt1*phi-vtt2)]
  !< type 4: harmonic type
  !< V(psi) = vtt0*(psi-vtt1)^2
  !< type 5: OPLS potential (four terms)
  !< V(phi) = vtt0 + vtt1*[1+cos(phi)] + vtt2*[1-cos(2*phi)] + vtt3*[1+cos(3*phi)] + vtt4*[1-cos(4*phi)]
  !< type 6: nine-term Fourier cosine series
  !< V(phi) = vtt0 + vtt1*cos(phi) + vtt2*cos(2*phi) + ... + vtt9*cos(9*phi)
  !< type 7: Ryckaert-Bellemans potential (three terms)
  !< V(psi) = vtt0 + vtt1*cos(psi) + vtt2*cos(psi)^2 + vtt3*cos(psi)^3
  !< type 8: Ryckaert-Bellemans potential (four terms)
  !< V(psi) = vtt0 + vtt1*cos(psi) + vtt2*cos(psi)^2 + vtt3*cos(psi)^4 + vtt3*cos(psi)^4
  !< \note psi uses polymer convention (trans is 0 deg), while phi uses protein convention (trans is 180 deg)
  real,allocatable::vtt(:,:)

  integer,allocatable::vibsplits(:),bendsplits(:),splpnts(:)
  real,allocatable::vib(:,:),tabvib(:,:),bend(:,:),tabbend(:,:),deg(:,:),tabtorso(:,:),torderiv2(:,:)
  integer::ntabvib,ntabbend,nttor

  type(LookupTable)::bonds,angles,dihedrals

  real,allocatable::rxvec(:,:),ryvec(:,:),rzvec(:,:),distij2(:,:),distanceij(:,:)
contains
  subroutine read_ff_bonded(io_ff)
    use util_search,only:tightenTable
    use util_mp,only:mp_bcast
    use sim_system,only:rootid,groupid
    integer,intent(in)::io_ff
    integer,parameter::initial_size=20
    integer::i,n,jerr,readstat
    character(LEN=default_string_length)::line_in

    !> Looking for section BONDS
    if (allocated(bonds%list)) then
       call destroyTable(bonds)
       deallocate(vib_type,brvib,brvibk,minRegrowVib,maxRegrowVib,stat=jerr)
    end if
    n=0
    if (myid.eq.rootid) then
       REWIND(io_ff)
       CYCLE_READ_BONDS:DO
          call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_bonds

          if (UPPERCASE(line_in(1:5)).eq.'BONDS') then
             call initiateTable(bonds,initial_size)
             allocate(vib_type(1:initial_size),brvib(1:initial_size),brvibk(1:initial_size),minRegrowVib(1:initial_size),maxRegrowVib(1:initial_size),stat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_intramolecular: bonds allocation failed',myid)
             brvib=0.0E0_dp
             brvibk=0.0E0_dp
             maxRegrowVib=2.0E0_dp
             minRegrowVib=0.0E0_dp
             do
                call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section BONDS',jerr)
                if (UPPERCASE(line_in(1:9)).eq.'END BONDS') exit
                n=n+1
                read(line_in,*) i
                i=addToTable(bonds,i,expand=.true.)
                if (i.gt.ubound(vib_type,1)) then
                   call reallocate(vib_type,1,2*ubound(vib_type,1))
                   call reallocate(brvib,1,2*ubound(brvib,1))
                   call reallocate(brvibk,1,2*ubound(brvibk,1))
                   call reallocate(minRegrowVib,1,2*ubound(minRegrowVib,1))
                   call reallocate(maxRegrowVib,1,2*ubound(maxRegrowVib,1))
                end if
                readstat=0
                call stripComments(line_in)
                read(line_in,*,iostat=readstat) jerr,vib_type(i),brvib(i),brvibk(i),minRegrowVib(i),maxRegrowVib(i)
                ! if the max and min regrows are not specified, you will get an
                ! error in iostat. The following is a dirty way to keep
                ! backwards compability if the max and min regrows are not
                ! specified
                if ((readstat .eq. -1 ) .or. (readstat .eq. -2)) then ! if the max and min regrows are not specified
                   read(line_in,*) jerr,vib_type(i),brvib(i),brvibk(i)
                   minRegrowVib(i)=0.0E0_dp
                   maxRegrowVib(i)=2.0E0_dp
                else if(.not.(readstat.eq.0)) then
                    call err_exit(__FILE__,__LINE__,'error in reading bond values',readstat)
                end if
             end do
             exit cycle_read_bonds
          end if
       END DO CYCLE_READ_BONDS
    end if

    call mp_bcast(n,1,rootid,groupid)

    if (n.gt.0) then
       if (myid.eq.rootid) then
          call tightenTable(bonds)
          call reallocate(vib_type,1,n)
          call reallocate(brvib,1,n)
          call reallocate(brvibk,1,n)
          call reallocate(minRegrowVib,1,n)
          call reallocate(maxRegrowVib,1,n)
       else
          call initiateTable(bonds,n)
          allocate(vib_type(1:n),brvib(1:n),brvibk(1:n),minRegrowVib(1:n),maxRegrowVib(1:n),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_intramolecular: bonds allocation failed',myid)
       end if

       call mp_bcast(bonds%size,1,rootid,groupid)
       call mp_bcast(bonds%list,bonds%size,rootid,groupid)
       call mp_bcast(vib_type,n,rootid,groupid)
       call mp_bcast(brvib,n,rootid,groupid)
       call mp_bcast(brvibk,n,rootid,groupid)
       call mp_bcast(minRegrowVib,n,rootid,groupid)
       call mp_bcast(maxRegrowVib,n,rootid,groupid)
    end if

    !> Looking for section ANGLES
    if (allocated(angles%list)) then
       call destroyTable(angles)
       deallocate(ben_type,brben,brbenk,stat=jerr)
    end if
    n=0
    if (myid.eq.rootid) then
       REWIND(io_ff)
       CYCLE_READ_ANGLES:DO
          call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_angles

          if (UPPERCASE(line_in(1:6)).eq.'ANGLES') then
             call initiateTable(angles,initial_size)
             allocate(ben_type(1:initial_size),brben(1:initial_size),brbenk(1:initial_size),stat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_intramolecular: angles allocation failed',myid)
             brben=0.0E0_dp
             brbenk=0.0E0_dp
             do
                call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section ANGLES',jerr)
                if (UPPERCASE(line_in(1:10)).eq.'END ANGLES') exit
                n=n+1
                read(line_in,*) i
                i=addToTable(angles,i,expand=.true.)
                if (i.gt.ubound(ben_type,1)) then
                   call reallocate(ben_type,1,2*ubound(ben_type,1))
                   call reallocate(brben,1,2*ubound(brben,1))
                   call reallocate(brbenk,1,2*ubound(brbenk,1))
                end if
                read(line_in,*) jerr,ben_type(i),brben(i),brbenk(i)
                brben(i) = brben(i) * onepi / 180.0E0_dp
             end do
             exit cycle_read_angles
          end if
       END DO CYCLE_READ_ANGLES
    end if

    call mp_bcast(n,1,rootid,groupid)

    if (n.gt.0) then
       if (myid.eq.rootid) then
          call tightenTable(angles)
          call reallocate(ben_type,1,n)
          call reallocate(brben,1,n)
          call reallocate(brbenk,1,n)
       else
          call initiateTable(angles,n)
          allocate(ben_type(1:n),brben(1:n),brbenk(1:n),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_intramolecular: angles allocation failed',myid)
       end if

       call mp_bcast(angles%size,1,rootid,groupid)
       call mp_bcast(angles%list,angles%size,rootid,groupid)
       call mp_bcast(ben_type,n,rootid,groupid)
       call mp_bcast(brben,n,rootid,groupid)
       call mp_bcast(brbenk,n,rootid,groupid)
    end if

    !> Looking for section DIHEDRALS
    if (allocated(dihedrals%list)) then
       call destroyTable(dihedrals)
       deallocate(vtt,torsion_type,stat=jerr)
    end if
    n=0
    if (myid.eq.rootid) then
       REWIND(io_ff)
       CYCLE_READ_DIHEDRALS:DO
          call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
          if (jerr.ne.0) exit cycle_read_dihedrals

          if (UPPERCASE(line_in(1:9)).eq.'DIHEDRALS') then
             call initiateTable(dihedrals,initial_size)
             allocate(vtt(0:9,1:initial_size),torsion_type(1:initial_size),stat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_intramolecular: dihedrals allocation failed',myid)
             vtt=0.0E0_dp
             do
                call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section DIHEDRALS',jerr)
                if (UPPERCASE(line_in(1:13)).eq.'END DIHEDRALS') exit
                n=n+1
                read(line_in,*) i
                i=addToTable(dihedrals,i,expand=.true.)
                if (i.gt.UBOUND(torsion_type,1)) then
                   call reallocate(torsion_type,1,2*UBOUND(torsion_type,1))
                   call reallocate(vtt,0,9,1,2*UBOUND(vtt,2))
                end if
                read(line_in,*) jerr,torsion_type(i),vtt(0:torsion_nParameter(torsion_type(i))-1,i)
                if (torsion_type(i).eq.3) then
                   vtt(2,i)=vtt(2,i)*onepi/180.0E0_dp
                else if (torsion_type(i).eq.1.or.torsion_type(i).eq.5) then
                   ! convert OPLS to Ryckaert-Bellemans because the latter is more efficient
                   vtt(0,i)=vtt(0,i)+vtt(1,i)+2.0E0_dp*vtt(2,i)+vtt(3,i)
                   vtt(1,i)=-vtt(1,i)+3.0E0_dp*vtt(3,i)
                   vtt(2,i)=-2.0E0_dp*vtt(2,i)+8.0E0_dp*vtt(4,i)
                   vtt(3,i)=-4.0E0_dp*vtt(3,i)
                   vtt(4,i)=-8.0E0_dp*vtt(4,i)
                   if (torsion_type(i).eq.1) then
                      torsion_type(i)=7
                   else if (torsion_type(i).eq.5) then
                      torsion_type(i)=8
                   end if
                end if
             end do
             exit cycle_read_dihedrals
          end if
       END DO CYCLE_READ_DIHEDRALS
    end if

    call mp_bcast(n,1,rootid,groupid)

    if (n.gt.0) then
       if (myid.eq.rootid) then
          call tightenTable(dihedrals)
          call reallocate(torsion_type,1,n)
          call reallocate(vtt,0,9,1,n)
       else
          call initiateTable(dihedrals,n)
          allocate(vtt(0:9,1:n),torsion_type(1:n),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_intramolecular: dihedrals allocation failed',myid)
       end if

       call mp_bcast(dihedrals%size,1,rootid,groupid)
       call mp_bcast(dihedrals%list,dihedrals%size,rootid,groupid)
       call mp_bcast(torsion_type,n,rootid,groupid)
       call mp_bcast(vtt,10*n,rootid,groupid)
    end if

    call read_tabulated_ff_bonded()

    return
  end subroutine read_ff_bonded

!DEC$ ATTRIBUTES FORCEINLINE :: vtorso
  function vtorso(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,itype)
    real::vtorso
    real,intent(in)::xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3
    integer,intent(in)::itype
    real::thetac,theta,tac2,tac3,tac4,tac5,tac6,tac7,tac8,tac9

    if (torsion_type(itype).eq.0) then
       ! type 0: dummy torsion type for setting up interaction table
       vtorso=0.0E0_dp
    else if (torsion_type(itype).eq.-1) then
       call dihedral_angle(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,thetac,theta,.true.)
       vtorso=inter_tor(theta,itype)
    else
       call dihedral_angle(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,thetac,theta,.false.)

       if (torsion_type(itype).eq.7) then
          !< type 7: Ryckaert-Bellemans potential (three terms), angle in polymer convention (trans is 0 deg)
          tac2 = thetac*thetac
          tac3 = tac2*thetac
          vtorso = vtt(0,itype)+vtt(1,itype)*thetac+vtt(2,itype)*tac2+vtt(3,itype)*tac3
       else if (torsion_type(itype).eq.8) then
          !< type 8: Ryckaert-Bellemans potential (four terms), angle in polymer convention (trans is 0 deg)
          tac2 = thetac*thetac
          tac3 = tac2*thetac
          tac4 = tac3*thetac
          vtorso = vtt(0,itype)+vtt(1,itype)*thetac+vtt(2,itype)*tac2+vtt(3,itype)*tac3+vtt(4,itype)*tac4
       else if (torsion_type(itype).eq.2) then
          ! type 2: Ryckaert-Bellemans potential, angle in polymer convention (trans is 0 deg)
          tac2 = thetac*thetac
          tac3 = tac2*thetac
          tac4 = tac3*thetac
          tac5 = tac4*thetac
          tac6 = tac5*thetac
          tac7 = tac6*thetac
          tac8 = tac7*thetac
          tac9 = tac8*thetac
          vtorso = vtt(0,itype)+vtt(1,itype)*thetac+vtt(2,itype)*tac2+vtt(3,itype)*tac3+vtt(4,itype)*tac4+vtt(5,itype)*tac5&
           +vtt(6,itype)*tac6+vtt(7,itype)*tac7+vtt(8,itype)*tac8+vtt(9,itype)*tac9
       else if (torsion_type(itype).eq.3) then
          ! type 3: periodic type, angle in protein convention (trans is 180 deg)
          theta=theta+onepi
          vtorso=vtt(0,itype)*(1+cos(vtt(1,itype)*theta-vtt(2,itype)))
       else if (torsion_type(itype).eq.4) then
          ! type 4: harmonic type, angle in polymer convention (trans is 0 deg)
          tac2=theta-vtt(1,itype)
          vtorso=vtt(0,itype)*tac2*tac2
       else if (torsion_type(itype).eq.1) then
          ! type 1: OPLS potential (three terms), angle in protein convention (trans is 180 deg)
          theta=theta+onepi
          ! remember: 1 + cos( theta+onepi ) = 1 - cos( theta )
          vtorso = vtt(0,itype) + vtt(1,itype)*(1.0E0_dp-thetac) + vtt(2,itype)*(1.E0_dp-cos(2.E0_dp*theta))&
           + vtt(3,itype)*(1.E0_dp+cos(3.E0_dp*theta))
       else if (torsion_type(itype).eq.5) then
          ! type 5: OPLS potential (four terms), angle in protein convention (trans is 180 deg)
          theta=theta+onepi
          vtorso = vtt(0,itype) + vtt(1,itype)*(1.0E0_dp-thetac) + vtt(2,itype)*(1.E0_dp-cos(2.E0_dp*theta))&
           + vtt(3,itype)*(1.E0_dp+cos(3.E0_dp*theta)) + vtt(4,itype)*(1.E0_dp-cos(4.E0_dp*theta))
       else if (torsion_type(itype).eq.6) then
          ! type 6: nine-term Fourier cosine series, angle in protein convention (trans is 180 deg)
          theta=theta+onepi
          vtorso=vtt(0,itype)-vtt(1,itype)*thetac+vtt(2,itype)*cos(2.0E0_dp*theta)+vtt(3,itype)*cos(3.0E0_dp*theta)&
           +vtt(4,itype)*cos(4.0E0_dp*theta)+vtt(5,itype)*cos(5.0E0_dp*theta)+vtt(6,itype)*cos(6.0E0_dp*theta)&
           +vtt(7,itype)*cos(7.0E0_dp*theta)+vtt(8,itype)*cos(8.0E0_dp*theta)+vtt(9,itype)*cos(9.0E0_dp*theta)
       else
          call err_exit(__FILE__,__LINE__,'vtorso: undefined torsional type',myid+1)
       end if
    end if

    return
  end function vtorso

!> \brief Calculate the dihedral angle and its cosine
!>
!> The dihedral is formed between vectors (1,2), (2,3), and (3,4) using polymer convention (trans is 0 degree)
!DEC$ ATTRIBUTES FORCEINLINE :: dihedral_angle
  subroutine dihedral_angle(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,thetac,theta,extended)
    real,intent(in)::xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3
    real,intent(out)::thetac,theta
    logical,intent(in)::extended !< whether to extend the dihedral angle to the range of -180 -- +180 degree

    real::x12,y12,z12,x23,y23,z23,d12,d23,dot,tcc,xcc,ycc,zcc

! calculate cross products d_a x d_a-1
    x12 = yvec1 * zvec2 - zvec1 * yvec2
    y12 = zvec1 * xvec2 - xvec1 * zvec2
    z12 = xvec1 * yvec2 - yvec1 * xvec2

! calculate cross products d_a-1 x d_a-2
    x23 = yvec2 * zvec3 - zvec2 * yvec3
    y23 = zvec2 * xvec3 - xvec2 * zvec3
    z23 = xvec2 * yvec3 - yvec2 * xvec3

! calculate lengths of cross products ***
    d12 = sqrt ( x12*x12 + y12*y12 + z12*z12 )
    d23 = sqrt ( x23*x23 + y23*y23 + z23*z23 )

! Addition for table look up for Torsion potential
! calculate dot product of cross products ***
    dot = x12*x23 + y12*y23 + z12*z23
    thetac = - (dot / ( d12 * d23 ))
    if (thetac.gt.1.0E0_dp) thetac=1.0E0_dp
    if (thetac.lt.-1.0E0_dp) thetac=-1.0E0_dp
    theta = acos(thetac)

    if (extended) then
       ! calculate cross product of cross products ***
       xcc = y12*z23 - z12*y23
       ycc = z12*x23 - x12*z23
       zcc = x12*y23 - y12*x23
       ! calculate scalar triple product ***
       tcc = xcc*xvec2 + ycc*yvec2 + zcc*zvec2
       ! determine angle between -180 and 180, not 0 to 180
       if (tcc .lt. 0.0E0_dp) theta = -theta
    end if

    return
  end subroutine dihedral_angle

!> \brief Calculates vibrational potential using linear interpolation between two points.
  function lininter_vib(r,typ) result(tabulated_vib)
    use util_math,only:polint
    use util_search,only:LOCATE
    real::tabulated_vib
    real,intent(in)::r
    integer,intent(in)::typ
    integer::low,high

    low=locate(vib(:,typ),vibsplits(typ),r,2)
    high=low+1
    if (vib(low,typ).gt.r.or.vib(high,typ).lt.r) then
       write(io_output,*) 'problem in lininter_vib!'
       write(io_output,*) 'len', r, ' vibtyp', typ
       write(io_output,*) 'low ', low, vib(low, typ)
       write(io_output,*) 'high ',high, vib(high, typ)
       write(io_output,*)
    end if
    call polint(vib(low:high,typ),tabvib(low:high,typ),2,r,tabulated_vib)
    return
  end function lininter_vib

!> \brief Calculates 1-3 nonbonded 'bending' potential using linear interpolation between two points.
  function lininter_bend(r,typ) result(tabulated_bend)
    use util_math,only:polint
    use util_search,only:LOCATE
    real::tabulated_bend
    real,intent(in)::r
    integer,intent(in)::typ
    integer::low,high

    low=locate(bend(:,typ),bendsplits(typ),r,2)
    high=low+1
    if (bend(low,typ).gt.r.or.bend(high,typ).lt.r) then
       write(io_output,*) 'problem in lininter_bend!'
       write(io_output,*) 'r', r, ' bendtyp', typ
       write(io_output,*) 'low ', low, bend(low, typ)
       write(io_output,*) 'high ',high, bend(high, typ)
       write(io_output,*)
    end if
    call polint(bend(low:high,typ),tabbend(low:high,typ),2,r,tabulated_bend)
    return
  end function lininter_bend

!> \brief Calculates torsional potential using interpolation.
  function inter_tor(thetarad,typ) result(tabulated_tor)
    use util_math,only:polint,splint
    use util_search,only:LOCATE
    real::tabulated_tor
    real,intent(in)::thetarad
    integer,intent(in)::typ
    real::theta
    integer::low,high

    theta=raddeg*thetarad
    if (L_spline) then
       ! spline interpolation
       call splint(deg(:,typ),tabtorso(:,typ),torderiv2(:,typ),splpnts(typ),theta,tabulated_tor)
    else if (L_linear) then
       ! linear interpolation between two points
       low=locate(deg(:,typ),splpnts(typ),theta,2)
       high=low+1
       if (deg(low,typ).gt.theta.or.deg(high,typ).lt.theta) then
          write(io_output,*) 'problem in inter_tor_linear!'
          write(io_output,*) 'theta ',thetarad, ' [rad], ',theta, ' [deg]. tortyp', typ
          write(io_output,*) 'low ', low, deg(low, typ)
          write(io_output,*) 'high ',high, deg(high, typ)
          write(io_output,*)
       end if
       call polint(deg(low:high,typ),tabtorso(low:high,typ),2,theta,tabulated_tor)
    end if
    return
  end function inter_tor

!> branched and linear molecules with connectivity table -
!> go through entire chain -
!> calculate all bonds vectors and lengths
!DEC$ ATTRIBUTES FORCEINLINE :: calc_connectivity
  subroutine calc_connectivity(i,imolty)
    use sim_system,only:nunit,nugrow,rxu,ryu,rzu,ijvib,invib
    integer,intent(in)::i,imolty

    real::rxui,ryui,rzui
    integer::ii,iivib,jj

    do ii = 1, nunit(imolty)
       rxui=rxu(i,ii)
       ryui=ryu(i,ii)
       rzui=rzu(i,ii)
       do iivib = 1, invib(imolty,ii)
          jj = ijvib(imolty,ii,iivib)
          rxvec(ii,jj) = rxu(i,jj) - rxui
          ryvec(ii,jj) = ryu(i,jj) - ryui
          rzvec(ii,jj) = rzu(i,jj) - rzui
          distij2(ii,jj) = ( rxvec(ii,jj)**2 + ryvec(ii,jj)**2 + rzvec(ii,jj)**2 )
          distanceij(ii,jj) = sqrt(distij2(ii,jj))

          if ( nunit(imolty) .ne. nugrow(imolty) )then
! account for explct atoms in opposite direction
             rxvec(jj,ii)   = -rxvec(ii,jj)
             ryvec(jj,ii)   = -ryvec(ii,jj)
             rzvec(jj,ii)   = -rzvec(ii,jj)
             distanceij(jj,ii) = distanceij(ii,jj)
          end if
       end do
    end do
  end subroutine calc_connectivity

!DEC$ ATTRIBUTES FORCEINLINE :: U_torsion
  function U_torsion(i,imolty,ist,lupdate_connectivity) result(vtg)
    use sim_system,only:nunit,intor,ittor,ijtor2,ijtor3,ijtor4
    real::vtg
    integer,intent(in)::i,imolty,ist
    logical,intent(in)::lupdate_connectivity

    integer::j,jjtor,ip1,ip2,ip3

    if (lupdate_connectivity) call calc_connectivity(i,imolty)

    vtg=0.0E0_dp
    do j = ist, nunit(imolty)
       do jjtor = 1, intor(imolty,j)
          ip3 = ijtor4(imolty,j,jjtor)
          if ( ip3 .lt. j ) then
             ip1 = ijtor2(imolty,j,jjtor)
             ip2 = ijtor3(imolty,j,jjtor)
             vtg = vtg + vtorso(rxvec(j,ip1),ryvec(j,ip1),rzvec(j,ip1),rxvec(ip1,ip2),ryvec(ip1,ip2),rzvec(ip1,ip2)&
              ,rxvec(ip2,ip3),ryvec(ip2,ip3),rzvec(ip2,ip3),ittor(imolty,j,jjtor))
          end if
       end do
    end do
  end function U_torsion

!> \brief calculate all stretching, bending, and torsional potentials
!> that have an end-bead with an index smaller than the current bead
!DEC$ ATTRIBUTES FORCEINLINE :: U_bonded
  subroutine U_bonded(i,imolty,vvib,vbend,vtg)
    use sim_system,only:nunit,invib,itvib,ijvib,inben,itben,ijben2,ijben3,L_vib_table,L_bend_table
    real,intent(out)::vvib,vbend,vtg
    integer,intent(in)::i,imolty

    real::theta,thetac,rbend,rbendsq
    integer::j,jjvib,ip1,ip2,it,jjben

    call calc_connectivity(i,imolty)

! stretching -
    vvib=0.0E0_dp
    do j = 2, nunit(imolty)
       do jjvib = 1, invib(imolty,j)
          ip1 = ijvib(imolty,j,jjvib)
          it  = itvib(imolty,j,jjvib)
          if ( ip1.lt. j .and. L_vib_table) then
             vvib = vvib + lininter_vib(distanceij(ip1,j),it)
! write(io_output,*) 'TABULATED VVIB: ', tabulated_vib,
!   &         distanceij(ip1,j), ip1, j
          end if
          if ( ip1 .lt. j .and..not.L_vib_table) vvib = vvib + brvibk(it) * (distanceij(ip1,j) - brvib(it))**2
       end do
    end do

! bending -
! molecule with bond bending
    vbend=0.0E0_dp
    do j = 2, nunit(imolty)
       do jjben = 1, inben(imolty,j)
          ip2 = ijben3(imolty,j,jjben)
          if ( ip2 .lt. j ) then
             ip1 = ijben2(imolty,j,jjben)
             it  = itben(imolty,j,jjben)
             if (brbenk(it).lt.-0.1E0_dp) cycle
             thetac = ( rxvec(ip1,j)*rxvec(ip1,ip2) + ryvec(ip1,j)*ryvec(ip1,ip2)&
              + rzvec(ip1,j)*rzvec(ip1,ip2) ) / ( distanceij(ip1,j)*distanceij(ip1,ip2) )
             if(thetac.gt.1.0E0_dp) then
                theta = 0.0E0_dp
             else if (thetac.lt.-1.0E0_dp) then
                theta = onepi
             else
                theta = acos(thetac)
             end if
             ! if bend table exists compute bending energy from that. Otherwise
             ! if we have a freely jointed chain then the energy is still zero.
             ! Otherwise compute the energy. Max value of freely jointed chain
             ! chosen to match the geometry subroutine in CBMC.
             if (L_bend_table) then
                rbendsq=distij2(ip1,j)+distij2(ip1,ip2)-2.0E0_dp*distanceij(ip1,j)*distanceij(ip1,ip2)*thetac
                rbend = sqrt(rbendsq)
                vbend = vbend + lininter_bend(rbend,it)
             else if (.not.(brbenk(it).lt.-0.1E0_dp)) then
                vbend = vbend +  brbenk(it) * (theta-brben(it))**2
             end if

          end if
       end do
    end do

! torsions -
! molecule with dihedral potenials ###
    vtg=U_torsion(i,imolty,2,.false.)

  end subroutine U_bonded

!> \brief Read in tabulated potential for bonded interactions (stretching, bending, and torsion) and set up linear interpolation
  subroutine read_table(file_tab,ntab,r,tab,splits,lists)
    use util_files,only:get_iounit
    use util_search,only:indexOf
    character(LEN=*),intent(in)::file_tab
    integer,intent(out)::ntab
    integer,allocatable,intent(inout)::splits(:)
    real,allocatable,intent(inout)::r(:,:),tab(:,:)
    type(LookupTable),intent(inout)::lists

    integer::io_tab,mmm,t,i,jerr

    io_tab=get_iounit()
    open(unit=io_tab,access='sequential',action='read',file=file_tab,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open tabulated potential file: '//file_tab,myid+1)

    read(io_tab,*) ntab
    do mmm=1,ntab
       read(io_tab,*) t
       t=indexOf(lists,t)
       if (t.eq.0) call err_exit(__FILE__,__LINE__,'read_table: undefined bead',jerr)
       i=1
       do
          if (i.gt.size(r,1)) then
             call reallocate(r,1,2*size(r,1),1,size(r,2))
             call reallocate(tab,1,2*size(tab,1),1,size(tab,2))
          end if
          read(io_tab,*,end=100) r(i,t), tab(i,t)
          if (r(i,t).eq.1000) exit
          ! write(io_tab+10,*) i,v(i,t),tab(i,t)
          i=i+1
       end do
100    splits(t)=i-1
       call reallocate(r,1,splits(t),1,size(r,2))
       call reallocate(tab,1,splits(t),1,size(tab,2))
    end do
    close(io_tab)
  end subroutine read_table

!> \brief Read in vibrational, 1-3 nonbonded 'bending', and torsional potential.
!> \since KM 12/02/08 vib
!> \since KM 12/03/08 1-3 'bending'
  subroutine read_tabulated_ff_bonded()
    use util_math,only:spline
    use sim_system,only:L_vib_table,L_bend_table
    integer,parameter::grid_size=1500
    integer::jerr,ttor

    if (dihedrals%size.gt.0) then
       if (ANY(torsion_type(1:dihedrals%size).eq.-1)) then
          if (L_spline.and.L_linear) call err_exit(__FILE__,__LINE__&
           ,'L_spline and L_linear should have one and only one to be true if using tabulated potential for torsions',myid+1)

          if (allocated(splpnts)) deallocate(splpnts,deg,tabtorso,stat=jerr)
          allocate(splpnts(1:dihedrals%size),deg(1:grid_size,1:dihedrals%size),tabtorso(1:grid_size,1:dihedrals%size),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__&
           ,'init_tabulated_potential_bonded: allocation failed for vib_table 1',myid+1)
          call read_table('fort.40',nttor,deg,tabtorso,splpnts,dihedrals)
          if (L_spline) then
             if (allocated(torderiv2)) deallocate(torderiv2,stat=jerr)
             allocate(torderiv2(1:grid_size,1:dihedrals%size),stat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__&
              ,'init_tabulated_potential_bonded: allocation failed for tor_table 2',myid+1)
             do ttor=1,nttor
                call spline(deg(:,ttor),tabtorso(:,ttor),splpnts(ttor),1.0E31_dp,1.0E31_dp,torderiv2(:,ttor))
             end do
          else if (.not.L_linear) then
             call err_exit(__FILE__,__LINE__&
              ,'Must set one of L_spline and L_linear to be true if using tabulated potential for torsions',myid+1)
          end if
       end if
    end if

    if (bonds%size.gt.0.and.L_vib_table) then
       if (allocated(vibsplits)) deallocate(vibsplits,vib,tabvib,stat=jerr)
       allocate(vibsplits(1:bonds%size),vib(1:grid_size,1:bonds%size),tabvib(1:grid_size,1:bonds%size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_bonded: allocation failed for vib_table',myid+1)
       call read_table('fort.41',ntabvib,vib,tabvib,vibsplits,bonds)
    end if

    if (angles%size.gt.0.and.L_bend_table) then
       if (allocated(bendsplits)) deallocate(bendsplits,bend,tabbend,stat=jerr)
       allocate(bendsplits(1:angles%size),bend(1:grid_size,1:angles%size),tabbend(1:grid_size,1:angles%size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_bonded: allocation failed for bend_table',myid+1)
       call read_table('fort.42',ntabbend,bend,tabbend,bendsplits,angles)
    end if
  end subroutine read_tabulated_ff_bonded

  subroutine allocate_energy_bonded()
    use sim_system,only:numax
    integer::jerr
    if (allocated(rxvec)) deallocate(rxvec,ryvec,rzvec,distij2,distanceij,stat=jerr)
    allocate(rxvec(numax,numax),ryvec(numax,numax),rzvec(numax,numax),distij2(numax,numax),distanceij(numax,numax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_energy_bonded: allocation failed',jerr)
    end if
  end subroutine allocate_energy_bonded
end MODULE energy_intramolecular
