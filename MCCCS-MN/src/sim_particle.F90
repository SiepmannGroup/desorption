MODULE sim_particle
  use var_type,only:dp,default_string_length
  use sim_system,only:lneigh,lneighbor,lgaro
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  private
  public::BeadType,AtomType,MoleculeType,check_neighbor_list,init_neighbor_list,add_neighbor_list,add_neighbor_list_molecule&
   ,update_neighbor_list_molecule,save_neighbor_list,restore_neighbor_list,allocate_neighbor_list,ctrmas,lnn,lnn_t,neighbor&
   ,neigh_cnt,ndij,nxij,nyij,nzij,neighi,neigh_icnt,ndiji,nxiji,nyiji,nziji,extend_molecule_type,mp_bcast_molecule&
   ,update_coord_in_tree

  type BeadType
     integer::type
     real::coord(3)
  end type BeadType

  type AtomType
     logical::lqchg,llj
     integer::atomID
     real::mass,charge,eps,sig
     character(LEN=default_string_length)::name& !<atomic symobl
      ,desc !<full description of the atom
  end type AtomType

  type MoleculeType
     integer::nbead
     type(BeadType),allocatable::bead(:)
  end type MoleculeType

  !> \brief Neighbor list implemented in adjacency matrix
  !> \details The purpose of this implementation is to speed up energy calculations by
  !> skipping over atoms that are not neighbors. Maximum displacements of translation
  !> and rotation are capped at upnn(sq) and upnndg, respectively.
  logical,allocatable::lnn(:,:),lnn_t(:,:)
  real,allocatable::upnn(:),upnnsq(:),upnndg(:,:) !< \bug It will likely not work correctly for complex molecules.


  !> \brief neighbor list implemented in adjacency list
  !> \details This neighbor list is used to find neighbors of a selected molecule,
  !> as in AVBMC or in the calculation of many-body interactions.
  integer,allocatable::neighbor(:,:),neigh_cnt(:),neighboro(:,:),neigh_o(:)
  real,allocatable::ndij(:,:),nxij(:,:),nyij(:,:),nzij(:,:),ndijo(:,:),nxijo(:,:),nyijo(:,:),nzijo(:,:)
  integer,parameter::maxneigh=20
  integer::neighi(maxneigh),neigh_icnt
  real::ndiji(maxneigh),nxiji(maxneigh),nyiji(maxneigh),nziji(maxneigh)

contains
  SUBROUTINE extend_molecule_type(p,i)
    type(MoleculeType),intent(inout)::p
    integer,intent(in)::i
    type(MoleculeType)::p_temp
    integer::array_size
    if (.not.allocated(p%bead)) then
       allocate(p%bead(i:i+999))
    else
       array_size=size(p%bead)
       if (i.le.array_size) return
       allocate(p_temp%bead(i+999))
       p_temp%bead(1:array_size) = p%bead
       deallocate(p%bead)
       call move_alloc(p_temp%bead,p%bead)
    endif
  END SUBROUTINE extend_molecule_type

  SUBROUTINE mp_bcast_molecule(msg,source,comm)
    use util_runtime,only:err_exit
    use util_mp,only:mp_bcast
    TYPE(MoleculeType) :: msg
    INTEGER, INTENT(IN) :: source, comm

#ifdef __MPI__
    INTEGER,PARAMETER::COUNT=2
    INTEGER,DIMENSION(COUNT)::lena,typa
    INTEGER::ierr,mpiMoleculeType
    INTEGER(KIND=MPI_ADDRESS_KIND)::baseaddr,loca(COUNT)

    ! call MPI_SIZEOF(msg%bead(1)%coord(1),sizeReal,ierr)
    ! if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_sizeof failed for msg%bead(1)%coord(1)',ierr)

    ! call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,sizeReal,mpiReal,ierr)
    ! if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_type_match_size failed for mpiReal',ierr)

    lena=(/1,3/)
#ifdef __DOUBLE_PRECISION__
    typa=(/MPI_INTEGER,MPI_DOUBLE_PRECISION/)
#else
    typa=(/MPI_INTEGER,MPI_REAL/)
#endif

    call mpi_get_address(msg%bead(1),baseaddr,ierr)
    if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_get_address failed for baseaddr',ierr)
    call mpi_get_address(msg%bead(1)%type,loca(1),ierr)
    if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_get_address failed for loca(1)',ierr)
    call mpi_get_address(msg%bead(1)%coord,loca(2),ierr)
    if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_get_address failed for loca(2)',ierr)
    loca=loca-baseaddr

    call MPI_TYPE_CREATE_STRUCT(COUNT,lena,loca,typa,mpiMoleculeType,ierr)
    if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_type_create_struct failed for mpiMoleculeType',ierr)

    call MPI_TYPE_COMMIT(mpiMoleculeType,ierr)
    if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_type_commit failed for mpiMoleculeType',ierr)

    call MPI_BCAST(msg%bead,msg%nbead,mpiMoleculeType,source,comm,ierr)
    if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'mp_bcast_molecule: mpi_bcast failed for msg%bead',ierr)

    call MPI_TYPE_FREE(mpiMoleculeType,ierr)
#endif
  END SUBROUTINE mp_bcast_molecule

  subroutine check_neighbor_list(ibox,imol)
    use sim_system,only:rmtrax,rmtray,rmtraz,rmrotx,rmroty,rmrotz,Armtrax,Armtray,Armtraz,io_output
    integer,intent(in)::ibox,imol

    if (rmtrax(imol,ibox).gt.upnn(ibox)) then
       write(io_output,*) ' rmtrax greater than upnn',ibox,imol
       rmtrax(imol,ibox) = upnn(ibox)
    end if
    if (rmtray(imol,ibox).gt.upnn(ibox)) then
       write(io_output,*) ' rmtray greater than upnn',ibox,imol
       rmtray(imol,ibox) = upnn(ibox)
    end if
    if (rmtraz(imol,ibox).gt.upnn(ibox)) then
       write(io_output,*) ' rmtraz greater than upnn',ibox,imol
       rmtraz(imol,ibox) = upnn(ibox)
    end if

    if (rmrotx(imol,ibox).gt.upnndg(imol,ibox)) then
       write(io_output,*) ' rmrotx greater than upnndg',ibox,imol
       rmrotx(imol,ibox) = upnndg(imol,ibox)
    end if
    if (rmroty(imol,ibox).gt.upnndg(imol,ibox)) then
       write(io_output,*) ' rmroty greater than upnndg',ibox,imol
       rmroty(imol,ibox) = upnndg(imol,ibox)
    end if
    if (rmrotz(imol,ibox).gt.upnndg(imol,ibox)) then
       write(io_output,*) ' rmrotz greater than upnndg',ibox,imol
       rmrotz(imol,ibox) = upnndg(imol,ibox)
    end if

    if (Armtrax.gt.upnn(ibox)) then
       Armtrax = upnn(ibox)
       write(io_output,*) '### problem : for target accept ', 'ratio Armtrax should be smaller than upnn'
    end if
    if (Armtray.gt.upnn(ibox)) then
       Armtray = upnn(ibox)
       write(io_output,*) '### problem : for target accept ', 'ratio Armtray should be smaller than upnn'
    end if
    if (Armtraz.gt.upnn(ibox)) then
       Armtraz = upnn(ibox)
       write(io_output,*) '### problem : for target accept ', 'ratio Armtraz should be smaller than upnn'
    end if
  end subroutine check_neighbor_list

  subroutine init_neighbor_list(ibox)
    use sim_system,only:nchain,nboxi,rcut,rcutnn,nmolty,nunit,brvib
    integer,intent(in)::ibox
    integer::imol,j
    real::umatch

    do j=1,nchain
       if (nboxi(j).eq.ibox) then
          ! set j-part of logical map to .false.
          lnn(j,:)=.false.
          lnn(:,j)=.false.

          ! reset adjacency list of chain j
          neigh_cnt(j) = 0
       end if
    end do

    if (lneigh) then
       ! calculate max. angular displacement that doesn't violate upnn
       ! calculate max. all-trans chain length ( umatch )
       upnn(ibox) = ( rcutnn(ibox) - rcut(ibox) ) / 3.0E0_dp
       upnnsq(ibox)=upnn(ibox)*upnn(ibox)
       do imol=1,nmolty
          umatch = 0.0E0_dp
          do j = 1, nunit(imol) - 1
             umatch = umatch + brvib(1)
          end do
          upnndg(imol,ibox) = asin( upnn(ibox) / umatch )

          call check_neighbor_list(ibox,imol)
       end do
    end if
  end subroutine init_neighbor_list

  subroutine add_neighbor_list(i,j,rij,xij,yij,zij)
    integer,intent(in)::i,j
    real,intent(in)::rij,xij,yij,zij

    neigh_cnt(i)=neigh_cnt(i)+1
    neighbor(neigh_cnt(i),i)=j
    neigh_cnt(j)=neigh_cnt(j)+1
    neighbor(neigh_cnt(j),j)=i
    ndij(neigh_cnt(i),i) = rij
    ndij(neigh_cnt(j),j) = rij
    nxij(neigh_cnt(i),i) = xij
    nyij(neigh_cnt(i),i) = yij
    nzij(neigh_cnt(i),i) = zij
    nxij(neigh_cnt(j),j) = -xij
    nyij(neigh_cnt(j),j) = -yij
    nzij(neigh_cnt(j),j) = -zij
  end subroutine add_neighbor_list

  subroutine add_neighbor_list_molecule(j,rij,xij,yij,zij)
    integer,intent(in)::j
    real,intent(in)::rij,xij,yij,zij

    neigh_icnt=neigh_icnt+1
    neighi(neigh_icnt)=j
    ndiji(neigh_icnt) = rij
    nxiji(neigh_icnt) = xij
    nyiji(neigh_icnt) = yij
    nziji(neigh_icnt) = zij
  end subroutine add_neighbor_list_molecule

  subroutine update_neighbor_list_molecule(i)
    integer,intent(in)::i
    integer::ic,j,ip
    logical::lneighij

    if (lneigh) then
       ! set i-part of logical map to .false.
       lnn(i,:)=lnn_t(:,i)
       lnn(:,i)=lnn_t(:,i)
    end if

    if (lneighbor.or.lgaro) then
       do ic = 1, neigh_cnt(i)
          j = neighbor(ic,i)
          do ip = 1,neigh_cnt(j)
             if (neighbor(ip,j).eq.i) then
                neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                ndij(ip,j) = ndij(neigh_cnt(j),j)
                nxij(ip,j) = nxij(neigh_cnt(j),j)
                nyij(ip,j) = nyij(neigh_cnt(j),j)
                nzij(ip,j) = nzij(neigh_cnt(j),j)
                neigh_cnt(j) = neigh_cnt(j)-1
                exit
             end if
          end do
       end do
       neigh_cnt(i) = neigh_icnt
       do ic = 1,neigh_icnt
          j = neighi(ic)
          neighbor(ic,i)=j
          ndij(ic,i) = ndiji(ic)
          nxij(ic,i) = nxiji(ic)
          nyij(ic,i) = nyiji(ic)
          nzij(ic,i) = nziji(ic)
          lneighij = .false.
          do ip = 1,neigh_cnt(j)
             if (neighbor(ip,j).eq.i) then
                lneighij = .true.
             end if
          end do
          if (.not.lneighij) then
             neigh_cnt(j) = neigh_cnt(j)+1
             neighbor(neigh_cnt(j),j) = i
             ndij(neigh_cnt(j),j) = ndiji(ic)
             nxij(neigh_cnt(j),j) = -nxiji(ic)
             nyij(neigh_cnt(j),j) = -nyiji(ic)
             nzij(neigh_cnt(j),j) = -nziji(ic)
          end if
       end do
    end if
  end subroutine update_neighbor_list_molecule

  subroutine save_neighbor_list(i)
    integer,intent(in)::i
    integer::j

    if (lneigh) lnn_t(:,i)=lnn(:,i)

    if (lneighbor.or.lgaro) then
       neigh_o(i) = neigh_cnt(i)
       do j=1,neigh_o(i)
          neighboro(j,i) = neighbor(j,i)
          ndijo(j,i) = ndij(j,i)
          nxijo(j,i) = nxij(j,i)
          nyijo(j,i) = nyij(j,i)
          nzijo(j,i) = nzij(j,i)
       end do
    end if
  end subroutine save_neighbor_list

  subroutine restore_neighbor_list(i)
    integer,intent(in)::i
    integer::j

    if (lneigh) then
       lnn(:,i)=lnn_t(:,i)
       lnn(i,:)=lnn_t(:,i)
    end if

    if (lneighbor.or.lgaro) then
       neigh_cnt(i) = neigh_o(i)
       do j=1,neigh_cnt(i)
          neighbor(j,i) = neighboro(j,i)
          ndij(j,i) = ndijo(j,i)
          nxij(j,i) = nxijo(j,i)
          nyij(j,i) = nyijo(j,i)
          nzij(j,i) = nzijo(j,i)
       end do
    end if
  end subroutine restore_neighbor_list

  subroutine allocate_neighbor_list()
    use util_runtime,only:err_exit
    use util_string,only:integer_to_string
    use sim_system,only:nbox,rcut,rcutnn,nbxmax,ntmax,nmax
    integer::ibox,jerr

    if (lneigh) then
       do ibox = 1,nbox
          ! If using neighbour list make sure the rcut & rcutnn is the same
          ! for all the boxes
          !> \bug Needs to be checked. Does not seem necessary.
          if ((abs(rcut(1)-rcut(ibox)).gt.1.0E-10_dp).and.(abs(rcutnn(1)-rcutnn(ibox)).gt.1.0E-10_dp)) then
             call err_exit(__FILE__,__LINE__,'Keep rcut and rcutnn for all the boxes same',-1)
          end if

          if (rcut(ibox).ge.rcutnn(ibox)) then
             call err_exit(__FILE__,__LINE__,' rcut greater equal rcutnn for box'//integer_to_string(ibox),-1)
          end if
       end do
    end if

    if (allocated(lnn)) deallocate(lnn,lnn_t,upnn,upnnsq,upnndg,neighbor,neigh_cnt,neighboro,neigh_o,ndij,nxij,nyij,nzij,ndijo,nxijo,nyijo,nzijo,stat=jerr)
    allocate(lnn(nmax,nmax),lnn_t(nmax,nmax),upnn(nbxmax),upnnsq(nbxmax),upnndg(ntmax,nbxmax),neighbor(maxneigh,nmax)&
     ,neigh_cnt(nmax),neighboro(maxneigh,nmax),neigh_o(nmax),ndij(maxneigh,nmax),nxij(maxneigh,nmax),nyij(maxneigh,nmax)&
     ,nzij(maxneigh,nmax),ndijo(maxneigh,nmax),nxijo(maxneigh,nmax),nyijo(maxneigh,nmax),nzijo(maxneigh,nmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'allocate_neighbor_list: allocation failed',jerr)
  end subroutine allocate_neighbor_list

!> \brief Find the center of mass of a chain and returns it to
!> the periodic box if it has left the box.
!> \param lall controls whether this is done for just chain \a j or for the entire box
!> \param ibox is the box of the particle
!> \param ldx is whether or not the COM > x-direction boxlength
!> \param ldy is whether or not the COM > y-direction boxlength
!> \param ldz is whether or not the COM > z-direction boxlength
!> \param j is the particle number
!> \param mtype 1 = calling from translation, should need only 1 operation to fold back,
!> unless the simulation box is non-orthorhombic in which case 2 operations may be needed
!> 2 = calling from rotation, 2 operations \n
!> 3 = calling from swap, 0 operation needed \n
!> 4 = N/A \n
!> 5 = calling from volume moves, 1 operation \n
!> 6 = calling from readdat, any number of operations \n
!> 7 = calling from config, 2 operations \n
!> 8 = calling from swatch, 2 opeartions \n
!> 9 = calling from energy for the dummy atom, 2 operations \n
!> 10= calling from atom_translation, 1 operation
!> j_box, optional, used for COM-kdtree, when updating molecule coordinates for swap/swatch moves,
!>          the old box information is still needed
  subroutine ctrmas(lall,ibox,j,mtype, jbox)
    use util_runtime,only:err_exit
    use sim_system
    use sim_cell

    logical,intent(in)::lall
    integer,intent(in)::ibox,j,mtype
    integer,intent(in),optional::jbox

    integer::stt,edd,i,imolty,iunit,iadjust,ii,itype,inboxx,inboxy,inboxz,iwarn,old_box
    real::rbx,rby,rbz,nxcm,nycm,nzcm,dx,dy,dz,sx,sy,sz,nxcm2,nycm2,nzcm2,dmaxsq,rxuij,ryuij,rzuij,rijsq
    logical::lintbx,ldx,ldy,ldz

    if (present(jbox)) then
        old_box = jbox
    else
        old_box = ibox
    end if

    rbx = boxlx(ibox)
    rby = boxly(ibox)
    rbz = boxlz(ibox)

    if (lall) then
       stt = 1
       edd = nchain
    else
       stt = j
       edd = j
    end if

    if ((mtype.eq.1).or.(mtype.eq.2).or.(mtype.eq.5)) then
       if (mtype.eq.1.and.(lsolid(ibox).and..not.lrect(ibox))) then
          iwarn = 2
       else if (mtype.eq.2) then
          ! kea 6/3/09 --- necessary for non-COM rotations
          iwarn = 2
       else
          iwarn = 1
       end if
    else if (mtype.eq.6) then
       iwarn = 0
    else
       iwarn = 2
    end if

    do i = stt, edd
       imolty = moltyp(i)
       iunit = nunit(imolty)
       iadjust = 1

       ! Check if the chain i is in the correct box
       if (nboxi(i) .eq. ibox) then
          ! Determine new center of mass for chain i
          lintbx = .false.
25        nxcm = 0.0E0_dp
          nycm = 0.0E0_dp
          nzcm = 0.0E0_dp
          do ii = 1, iunit
             itype = ntype(imolty,ii)
             nxcm = nxcm + rxu(i,ii) * mass(itype)
             nycm = nycm + ryu(i,ii) * mass(itype)
             nzcm = nzcm + rzu(i,ii) * mass(itype)
          end do
          nxcm = nxcm / masst(imolty)
          nycm = nycm / masst(imolty)
          nzcm = nzcm / masst(imolty)

          ldx = .false.
          ldy = .false.
          ldz = .false.
          dx = 0.0E0_dp
          dy = 0.0E0_dp
          dz = 0.0E0_dp

          if (lsolid(ibox) .and. .not. lrect(ibox)) then
             sx = nxcm*hmati(ibox,1)+nycm*hmati(ibox,4) +nzcm*hmati(ibox,7)
             sy = nxcm*hmati(ibox,2)+nycm*hmati(ibox,5) +nzcm*hmati(ibox,8)
             sz = nxcm*hmati(ibox,3)+nycm*hmati(ibox,6) +nzcm*hmati(ibox,9)

             if ( sx .lt. -1.0E-10_dp ) then
                sx = sx + 1.0E0_dp
                ldx = .true.
             else if ( sx .gt. 1E0_dp ) then
                sx = sx - 1.0E0_dp
                ldx = .true.
             end if
             if ( sy .lt. -1.0E-10_dp ) then
                sy = sy + 1.0E0_dp
                ldy = .true.
             else if ( sy .gt. 1E0_dp ) then
                sy = sy - 1.0E0_dp
                ldy = .true.
             end if
             if ( sz .lt. -1.0E-10_dp ) then
                sz = sz + 1.0E0_dp
                ldz = .true.
             else if ( sz .gt. 1E0_dp ) then
                sz = sz - 1.0E0_dp
                ldz = .true.
             end if
             sxcm(i) = sx
             sycm(i) = sy
             szcm(i) = sz

             if ( ldx .or. ldy .or. ldz ) then
                if ( mtype .eq. 5 ) then
                   write(io_output,*) 'sx, sy, sz:',sx,sy,sz
                end if
                nxcm2 = sx*hmat(ibox,1)+sy*hmat(ibox,4) +sz*hmat(ibox,7)
                nycm2 = sx*hmat(ibox,2)+sy*hmat(ibox,5) +sz*hmat(ibox,8)
                nzcm2 = sx*hmat(ibox,3)+sy*hmat(ibox,6) +sz*hmat(ibox,9)
                dx = nxcm2-nxcm
                dy = nycm2-nycm
                dz = nzcm2-nzcm

                do ii = 1, iunit
                   rxu(i,ii) = rxu(i,ii) + dx
                   ryu(i,ii) = ryu(i,ii) + dy
                   rzu(i,ii) = rzu(i,ii) + dz
                end do
                nxcm = nxcm2
                nycm = nycm2
                nzcm = nzcm2
             end if
          else
             if ( lintbx ) then
                if ( nxcm .gt. rbx ) then
                   inboxx = int(nxcm/rbx)
                   ldx = .true.
                else if ( nxcm .lt. -1.0E-10_dp ) then
                   inboxx = int(nxcm/rbx) - 1
                   ldx = .true.
                end if
                if ( nycm .gt. rby ) then
                   inboxy = int(nycm/rby)
                   ldy = .true.
                else if ( nycm .lt. -1.0E-10_dp )then
                   inboxy = int(nycm/rby) - 1
                   ldy = .true.
                end if
                if (lpbcz) then
                   if ( nzcm .gt. rbz ) then
                      inboxz = int(nzcm/rbz)
                      ldz = .true.
                   else if ( nzcm .lt. -1.0E-10_dp ) then
                      inboxz = int(nzcm/rbz) - 1
                      ldz = .true.
                   end if
                end if
                if ( ldx ) then
                   dx = -real(inboxx,dp)*rbx
                end if
                if ( ldy ) then
                   dy = -real(inboxy,dp)*rby
                end if
                if ( ldz ) then
                   dz = -real(inboxz,dp)*rbz
                end if
             else
                if (nxcm .lt. -1.0E-10_dp) then
                   dx = rbx
                   ldx = .true.
                else if (nxcm .gt. rbx) then
                   dx = -rbx
                   ldx = .true.
                end if
                if (nycm .lt. -1.0E-10_dp) then
                   dy = rby
                   ldy = .true.
                else if (nycm .gt. rby) then
                   dy = -rby
                   ldy = .true.
                end if
                if (lpbcz) then
                   if (nzcm .lt. -1.0E-10_dp) then
                      dz = rbz
                      ldz = .true.
                   else if (nzcm .gt. rbz) then
                      dz = -rbz
                      ldz = .true.
                   end if
                end if
             end if

             ! Update coordinates in kdtree
             if ((ldx .or. ldy .or. ldz) .and. (.not. lcutcm) .and. lkdtree .and. lkdtree_box(ibox)) then
                 do ii = 1, iunit
                     rxu_update(ii) = rxu(i, ii) + dx
                     ryu_update(ii) = ryu(i, ii) + dy
                     rzu_update(ii) = rzu(i, ii) + dz
                 end do

                 if ((mtype .eq. 5) .or. (mtype .eq. 9)) then
                     ! if energy call or volume move, do not update the kdtree
                     ! for energy call, it's a fictious bead
                     ! for volume move, we'll reconstruct the tree later
                     call update_coord_in_tree(i, iunit, ibox, ibox, .false., .true.)
                 else
                     call update_coord_in_tree(i, iunit, ibox, ibox, .true., .true.)
                 end if
             end if

             ! Update coordinates in r*xu array
             if (ldx) then
                do ii = 1, iunit
                   rxu(i,ii) = rxu(i,ii) + dx
                end do
             end if
             if (ldy) then
                do ii = 1, iunit
                   ryu(i,ii) = ryu(i,ii) + dy
                end do
             end if
             if (ldz) then
                do ii = 1, iunit
                   rzu(i,ii) = rzu(i,ii) + dz
                end do
             end if

          end if

          if (ldx .or. ldy .or. ldz ) then
             if ( (iadjust .ge. iwarn) ) then
                if (mtype .eq. 1) then
                   write(io_output,*) 'translational move'
                   if ((boxlx(ibox) .lt. rmtrax(imolty,ibox)) .or.&
                       (boxly(ibox) .lt. rmtray(imolty,ibox)) .or.&
                       (boxlz(ibox) .lt. rmtraz(imolty,ibox))) then
                       write(io_output,FMT='(A)') 'Future error possibly due to max displacement &
                                                    & for translation being > boxlength. &
                                                    & Consider decreasing iratio.'
                   end if
                end if
                if (mtype .eq. 2) write(io_output,*) 'rotational move'
                if (mtype .eq. 3) write(io_output,*) 'swap move'
                if (mtype .eq. 4) write(io_output,*) 'switch move'
                if (mtype .eq. 5) write(io_output,*) 'volume move'
                if (mtype .eq. 6) write(io_output,*) 'readdat move'
                if (mtype .eq. 7) write(io_output,*) 'config move'
                if (mtype .eq. 8) write(io_output,*) 'swatch move'
                if (mtype .eq. 9) write(io_output,*) 'energy call'
                write(io_output,*) 'ibox,i,iunit,boxlen',ibox,i,iunit,rbx,rby,rbz
                lintbx = .true.
                write(io_output,*) 'nxcm,nycm,nzcm',nxcm,nycm,nzcm
                write(io_output,*) 'dx,dy,dz',dx,dy,dz
                if (iwarn .ne. 0) call err_exit(__FILE__,__LINE__,'ctrmas failure',myid+1)
             end if
             iadjust = iadjust + 1
             goto 25
          end if

          ! Update COM-based kd-tree
          if ((.not. lall) .and. lcutcm .and. lkdtree .and. lkdtree_box(ibox) .and. (mtype .ne. 5) .and. (mtype .ne. 9)) then
              ! if energy call or volume move, do not update the kdtree
              ! for energy call, it's a fictious bead
              ! for volume move, we'll reconstruct the tree later
              rxu_update(1) = xcm(i)
              ryu_update(1) = ycm(i)
              rzu_update(1) = zcm(i)
          end if

          ! assign the new center of mass
          xcm(i) = nxcm
          ycm(i) = nycm
          zcm(i) = nzcm

          ! Update COM-based kd-tree if there is a coordinate change
          if ((.not. lall) .and. lcutcm .and. lkdtree .and. lkdtree_box(ibox) .and. (mtype .ne. 5) .and. (mtype .ne. 9)) then
              if ((abs(xcm(i)-rxu_update(1)) .gt. 1e-7) .or. &
                    (abs(ycm(i)-ryu_update(1)) .gt. 1e-7) .or. &
                    (abs(zcm(i)-rzu_update(1)) .gt. 1e-7)) then
                  call update_coord_in_tree(i, 1, old_box, ibox, .true., .true.)
              end if
          end if

          if ( (lcutcm .or. ldual) .and. iwarn .ne. 1 ) then
             dmaxsq = 0.0E0_dp
             do ii=1,iunit
                rxuij = rxu(i,ii)-nxcm
                ryuij = ryu(i,ii)-nycm
                rzuij = rzu(i,ii)-nzcm
                ! minimum image the ctrmas pair separations ---
                if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                if ( rijsq .gt. dmaxsq ) dmaxsq = rijsq
             end do
             rcmu(i) = sqrt(dmaxsq)+ 1.0E-10_dp
             ! write(io_output,*) 'rcmu(i)',rcmu(i)
          end if
       else if (.not. lall) then
          write(io_output,*) 'prob with box in ctrmas'
       end if
    end do
  end subroutine ctrmas

  ! Update the particle coordinates
  !< i: the index of the molecule
  !< iunit: the total number of beads
  !< box_to_delete: the box where the bead needs to be deleted from
  !< box_to_insert: the box where the bead needs to be inserted, should be the same as box_to_delete if it's NOT swap/swatch move
  !< l_update_tree: see explanation below
  !< lctrmas: see explanation below
  subroutine update_coord_in_tree(i, iunit, box_to_delete, box_to_insert, l_update_tree, lctrmas)
      use util_runtime,only:err_exit
      use sim_system
      use util_kdtree,only: construct_kdtree,search_node_in_tree,delete_node,insert_node, allocate_nxyz, check_tree_coord
      logical :: l_update_tree !< whether to update kdtree IF THE KDTREE IS USED
                             ! (not used when some ghost particle coordinates are updated; or the tree has not yet been initialized)
      logical, intent(in) :: lctrmas !< if called from ctrmas, delete all the beads first before insertion
                                     !< otherwise do it at the same time
      !< the molecule number; the number of units; and the box number to delete and insert
      integer, intent(in) :: i, iunit, box_to_delete, box_to_insert
      integer :: j, ibox
      real :: dx, dy, dz, x_coord, y_coord, z_coord

      !kdtree variables
      type(tree), pointer :: kd_tree
      integer, allocatable :: nx(:), ny(:), nz(:)
      integer :: ix, iy, iz, iDim
      real, dimension(3) :: old_coord, new_coord
      logical :: lAdd
      real :: xmin, xmax, ymin, ymax, zmin, zmax, rbcut_plus_buffer, xmin_all, xmax_all, ymin_all, ymax_all, zmin_all, zmax_all
      type(tree_node), pointer :: node_found

      ! if kdtree is used, ready to update kdtree first
      ! if kdtree is NOT used,
      ! or if it's the ctrmas call before kdtree is initialized, update only the coordinate
      if (.not. allocated(mol_tree)) return

      ! if both boxes don't use kdtree, don't need to update kdtree
      if ((.not. lkdtree_box(box_to_delete)) .and. (.not. lkdtree_box(box_to_insert))) return

      ! otherwise, kdtree needs to be updated if l_update_tree = T
      if (lkdtree .and. l_update_tree) then
         ! allocate nx, ny and nz used for later PCB use
         call allocate_nxyz(nx, ny, nz)

         ! Update beads
         ! Define tree related stuff regarding the box_to_delete
         ! Tree information for box_to_delete is updated here to avoid multiple occurrence of this part
         if (lkdtree_box(box_to_delete)) then
            ibox = box_to_delete
            kd_tree => mol_tree(ibox)%tree
            rbcut_plus_buffer = rcut(ibox) + kdtree_buffer_len(ibox)
            xmin = kd_tree%bound(1)%lower - rbcut_plus_buffer !< these are min and max's for old coordinates
            xmax = kd_tree%bound(1)%upper + rbcut_plus_buffer
            ymin = kd_tree%bound(2)%lower - rbcut_plus_buffer
            ymax = kd_tree%bound(2)%upper + rbcut_plus_buffer
            zmin = kd_tree%bound(3)%lower - rbcut_plus_buffer
            zmax = kd_tree%bound(3)%upper + rbcut_plus_buffer
            xmin_all = kd_tree%bound_all(1)%lower
            xmax_all = kd_tree%bound_all(1)%upper
            ymin_all = kd_tree%bound_all(2)%lower
            ymax_all = kd_tree%bound_all(2)%upper
            zmin_all = kd_tree%bound_all(3)%lower
            zmax_all = kd_tree%bound_all(3)%upper
         end if

         if ((.not. lctrmas) .and. (box_to_delete .eq. box_to_insert)) then
             ! if not called from ctrmas, or not called from swap/swatch, do the insertion and the deletion at the same time
             if (lkdtree_box(box_to_delete)) then
                 ! if kd-tree is used for box_to_delete and box_to_insert, delete and insert at the same time
                 ! note: at this point, lkdtree_box(box_to_delete) == lkdtree_box(box_to_insert) == T
                 ! start deletion and insertion
                 do ix = 1, size(nx)
                     dx = nx(ix) * boxlx(ibox)
                     do iy = 1, size(ny)
                         dy = ny(iy) * boxly(ibox)
                         do iz = 1, size(nz)
                             dz = nz(iz) * boxlz(ibox)
                             do j = 1, iunit

                                 if (lcutcm) then
                                     old_coord(1) = rxu_update(1) + dx
                                     old_coord(2) = ryu_update(1) + dy
                                     old_coord(3) = rzu_update(1) + dz
                                     new_coord(1) = xcm(i) + dx
                                     new_coord(2) = ycm(i) + dy
                                     new_coord(3) = zcm(i) + dz
                                 else
                                     old_coord(1) = rxu(i, j) + dx
                                     old_coord(2) = ryu(i, j) + dy
                                     old_coord(3) = rzu(i, j) + dz
                                     new_coord(1) = rxu_update(j) + dx
                                     new_coord(2) = ryu_update(j) + dy
                                     new_coord(3) = rzu_update(j) + dz
                                 end if

                                 ! Don't need to do anything if old and new coord are the same (e.g. in CBMC)
                                 if ((abs(old_coord(1)-new_coord(1)) .gt. 1e-10) .or. (abs(old_coord(2)-new_coord(2)) .gt. 1e-10) &
                                    .or. (abs(old_coord(3)-new_coord(3)) .gt. 1e-10)) then
                                     if ((ix .eq. 1) .and. (iy .eq. 1) .and. (iz .eq. 1)) then
                                         kd_tree => delete_node(kd_tree, old_coord)

                                         do iDim = 1, 3
                                             if (new_coord(iDim) .gt. kd_tree%bound(iDim)%upper) &
                                                kd_tree%bound(iDim)%upper = new_coord(iDim)
                                             if (new_coord(iDim) .lt. kd_tree%bound(iDim)%lower) &
                                                kd_tree%bound(iDim)%lower = new_coord(iDim)
                                         end do
                                         lAdd = .true.
                                     else
                                         if ((old_coord(1) .gt. xmin_all) .and. (old_coord(1) .lt. xmax_all) &
                                            .and. (old_coord(2) .gt. ymin_all) .and. (old_coord(2) .lt. ymax_all) &
                                            .and. (old_coord(3) .gt. zmin_all) .and. (old_coord(3) .lt. zmax_all)) then
                                             node_found => search_node_in_tree(kd_tree, kd_tree%tree_root, old_coord, 1)
                                             if (associated(node_found)) then
                                                 kd_tree => delete_node(kd_tree, old_coord)
                                             end if
                                         end if

                                         if ((new_coord(1) .gt. (kd_tree%bound(1)%lower-rbcut_plus_buffer)) &
                                             .and. (new_coord(1) .lt. (kd_tree%bound(1)%upper+rbcut_plus_buffer)) &
                                             .and. (new_coord(2) .gt. (kd_tree%bound(2)%lower-rbcut_plus_buffer)) &
                                             .and. (new_coord(2) .lt. (kd_tree%bound(2)%upper+rbcut_plus_buffer)) &
                                             .and. (new_coord(3) .gt. (kd_tree%bound(3)%lower-rbcut_plus_buffer)) &
                                             .and. (new_coord(3) .lt. (kd_tree%bound(3)%upper+rbcut_plus_buffer))) then

                                             lAdd = .true.
                                         else
                                             lAdd = .false.
                                         end if
                                     end if

                                     if (lAdd) then
                                         kd_tree => insert_node(kd_tree, new_coord, i, j, nx(ix), ny(iy), nz(iz))
                                     end if
                                 end if
                             end do !< j = 1, iunit
                         end do !< iz
                     end do !< iy
                 end do !< ix
             end if !< if lkdtree_box(box_to_delete)

         else
             ! if ctrmas, swap of swatch move
             ! do the deletion and insertion separately
             if (lkdtree_box(box_to_delete)) then
                  do j = 1, iunit

                     if (lcutcm) then
                         x_coord = rxu_update(1)
                         y_coord = ryu_update(1)
                         z_coord = rzu_update(1)
                     else
                         x_coord = rxu(i, j)
                         y_coord = ryu(i, j)
                         z_coord = rzu(i, j)
                     end if

                     do ix = 1, size(nx)
                         old_coord(1) = x_coord + nx(ix) * boxlx(ibox)
                         do iy = 1, size(ny)
                             old_coord(2) = y_coord + ny(iy) * boxly(ibox)
                             do iz = 1, size(nz)
                                 old_coord(3) = z_coord + nz(iz) * boxlz(ibox)

                                 if ((ix .eq. 1) .and. (iy .eq. 1) .and. (iz .eq. 1)) then
                                     kd_tree => delete_node(kd_tree, old_coord)
                                 else
                                     if ((old_coord(1) .gt. xmin_all) .and. (old_coord(1) .lt. xmax_all) &
                                         .and. (old_coord(2) .gt. ymin_all) .and. (old_coord(2) .lt. ymax_all) &
                                         .and. (old_coord(3) .gt. zmin_all) .and. (old_coord(3) .lt. zmax_all)) then
                                         node_found => search_node_in_tree(kd_tree, kd_tree%tree_root, old_coord, 1)
                                         if (associated(node_found)) kd_tree => delete_node(kd_tree, old_coord)

                                     end if
                                 end if
                             end do
                         end do
                     end do
                 end do
             end if

             ! insertion, only if the box_to_insert uses kdtree
             if (lkdtree_box(box_to_insert)) then
                 ! if swap/swatch, box_to_delete and box_to_insert can be different, update tree information
                 if (box_to_insert .ne. box_to_delete) then
                     ibox = box_to_insert
                     kd_tree => mol_tree(ibox)%tree
                     rbcut_plus_buffer = rcut(ibox) + kdtree_buffer_len(ibox)
                 end if

                 if (lcutcm) then
                     x_coord = xcm(i)
                     y_coord = ycm(i)
                     z_coord = zcm(i)
                 end if

                 do j = 1, iunit

                     if (.not. lcutcm) then
                         x_coord = rxu_update(j)
                         y_coord = ryu_update(j)
                         z_coord = rzu_update(j)
                     end if

                     do ix = 1, size(nx)
                         new_coord(1) = x_coord + nx(ix) * boxlx(ibox)
                         do iy = 1, size(ny)
                             new_coord(2) = y_coord + ny(iy) * boxly(ibox)
                             do iz = 1, size(nz)
                                 new_coord(3) = z_coord + nz(iz) * boxlz(ibox)

                                 if ((ix .eq. 1) .and. (iy .eq. 1) .and. (iz .eq. 1)) then
                                     if (.not. lcutcm) then
                                         do iDim = 1, 3
                                             if (new_coord(iDim) .gt. kd_tree%bound(iDim)%upper) &
                                                 kd_tree%bound(iDim)%upper = new_coord(iDim)
                                             if (new_coord(iDim) .lt. kd_tree%bound(iDim)%lower) &
                                                 kd_tree%bound(iDim)%lower = new_coord(iDim)
                                         end do
                                     end if

                                     kd_tree => insert_node(kd_tree, new_coord, i, j, nx(ix), ny(iy), nz(iz))
                                 else
                                     if ((new_coord(1) .gt. (kd_tree%bound(1)%lower-rbcut_plus_buffer)) &
                                        .and. (new_coord(1) .lt. (kd_tree%bound(1)%upper+rbcut_plus_buffer)) &
                                        .and. (new_coord(2) .gt. (kd_tree%bound(2)%lower-rbcut_plus_buffer)) &
                                        .and. (new_coord(2) .lt. (kd_tree%bound(2)%upper+rbcut_plus_buffer)) &
                                        .and. (new_coord(3) .gt. (kd_tree%bound(3)%lower-rbcut_plus_buffer)) &
                                        .and. (new_coord(3) .lt. (kd_tree%bound(3)%upper+rbcut_plus_buffer))) then

                                        kd_tree => insert_node(kd_tree, new_coord, i, j, nx(ix), ny(iy), nz(iz))

                                     end if
                                 end if
                             end do
                         end do
                     end do
                 end do
             end if !< if lkdtree_box(box_to_insert) == T

         end if

       end if !<  if lkdtree

  end subroutine update_coord_in_tree

end MODULE sim_particle
