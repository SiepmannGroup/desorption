module zeolite
  use var_type,only:dp,default_string_length
  use const_math,only:eps,onepi,twopi,fourpi
  use const_phys,only:qqfact,N_Avogadro
  use util_math,only:erfunc
  use util_runtime,only:err_exit
  use util_timings, only: time_now
  use util_files,only:get_iounit
  use sim_system
  use sim_cell
  use sim_particle
  use sim_zeolite
  use parser_pdb
  use parser_cssr
  use parser_cif
  use energy_kspace,only:calp
#ifdef __OPENMP__
  use omp_lib
#endif
  implicit none
  private
  save
  public::zeocoord,suzeo,exzeo,upper_limit_zeo

  real,parameter::overlapValue=1.0E+20_dp,LJScaling=2E4_dp
  integer,parameter::boxZeo=1

  logical,allocatable::lunitcell(:)
  logical::ltailcZeo=.true.,ltestztb=.false.,lpore_volume=.false.,lsurface_area=.false.,printztb=.false.

  integer::nlayermax,n_pieces_ztb=1,num_points_interpolation=4,volume_probe=124,volume_nsample=20,area_probe=124,area_nsample=100
  real,allocatable::my_zgrid(:,:,:),zgrid(:,:,:),egrid(:,:),yjtmp(:),yktmp(:),yltmp(:),xt(:),yt(:),zt(:)
  real::requiredPrecision=1.0E-2_dp,upperLimit=1.0E+5_dp,upper_limit_zeo
  character(LEN=default_path_length)::file_zeocoord='zeolite.cssr',file_ztb='zeolite.ztb',file_supercell=''

  type(MoleculeType)::zeo
  type(ZeoliteBeadType)::ztype
  type(CellMaskType)::zcell !< the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
  type(ZeoliteUnitCellGridType)::zunit
  type(ZeolitePotentialType)::zpot

contains
!> \brief Interface to old code using mask data types
  subroutine initZeo()
    integer::i,j
    zcell%solid=>lsolid(boxZeo)
    zcell%ortho=>lrect(boxZeo)
    zcell%pbc(1)=lpbcx
    zcell%pbc(2)=lpbcy
    zcell%pbc(3)=lpbcz
    zcell%cut=>rcut(boxZeo)
    zcell%vol=>cell_vol(boxZeo)
    zcell%calp=>calp(boxZeo)
    zcell%boxl(1)%val=>boxlx(boxZeo)
    zcell%boxl(2)%val=>boxly(boxZeo)
    zcell%boxl(3)%val=>boxlz(boxZeo)
    do i=1,3
       zcell%ang(i)%val=>cell_ang(boxZeo,i)
       zcell%height(i)%val=>min_width(boxZeo,i)
       do j=1,3
          zcell%hmat(j,i)%val=>hmat(boxZeo,3*(i-1)+j)
          zcell%hmati(j,i)%val=>hmati(boxZeo,3*(i-1)+j)
       end do
    end do
  end subroutine initZeo

  subroutine zeocoord(file_in,lprint)
    use util_search,only:indexOf
    use util_mp,only:mp_bcast
    character(LEN=*),intent(in)::file_in
    LOGICAL,INTENT(IN)::lprint

    integer::io_input,jerr

    namelist /zeolite_in/ file_zeocoord,dgr,file_supercell,file_ztb,n_pieces_ztb,requiredPrecision,num_points_interpolation,upperLimit,ltailcZeo&
     ,ltestztb,lpore_volume,volume_probe,volume_nsample,lsurface_area,area_probe,area_nsample,printztb

    if (myid.eq.rootid) then
       io_input=get_iounit()
       open(unit=io_input,access='sequential',action='read',file=file_in,form='formatted',iostat=jerr,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open zeolite input file '//trim(file_in),myid+1)

       read(UNIT=io_input,NML=zeolite_in,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: zeolite_in',jerr)
       close(io_input)
    end if
    upper_limit_zeo = upperLimit

    call mp_bcast(file_zeocoord,rootid,groupid)
    call mp_bcast(dgr,1,rootid,groupid)
    call mp_bcast(file_supercell,rootid,groupid)
    call mp_bcast(file_ztb,rootid,groupid)
    call mp_bcast(n_pieces_ztb,1,rootid,groupid)
    call mp_bcast(requiredPrecision,1,rootid,groupid)
    call mp_bcast(num_points_interpolation,1,rootid,groupid)
    call mp_bcast(upperLimit,1,rootid,groupid)
    call mp_bcast(upper_limit_zeo,1,rootid,groupid)
    call mp_bcast(ltailcZeo,1,rootid,groupid)
    call mp_bcast(ltestztb,1,rootid,groupid)
    call mp_bcast(lpore_volume,1,rootid,groupid)
    call mp_bcast(volume_probe,1,rootid,groupid)
    call mp_bcast(volume_nsample,1,rootid,groupid)
    call mp_bcast(lsurface_area,1,rootid,groupid)
    call mp_bcast(area_probe,1,rootid,groupid)
    call mp_bcast(area_nsample,1,rootid,groupid)

    if (lprint.and.(ltestztb.or.lpore_volume.or.lsurface_area)) then
       write(io_output,FMT='(A)',advance='no') 'zeocoord: will'
       if (ltestztb) write(io_output,FMT='(A)') ' test accuracy of tabulated potential;'
       if (lpore_volume) write(io_output,FMT='(A,I0,A,I0,A)') ' calculate pore volume using atom type ',volume_probe&
        ,' as probe, with ',volume_nsample,' sample points;'
       if (lsurface_area) write(io_output,FMT='(A,I0,A,I0,A)') ' calculate surface area using atom type ',area_probe&
        ,' as probe, with ',area_nsample,' sample points around each framework atom;'
    end if

    volume_probe=indexOf(atoms,volume_probe)
    area_probe=indexOf(atoms,area_probe)
    if (((ltestztb.or.lpore_volume).and.volume_probe.eq.0).or.(lsurface_area.and.area_probe.eq.0)) call err_exit(__FILE__,__LINE__&
     ,'zeocoord: atom parameters for volume_probe or area_probe undefined',myid+1)

    call initZeo()

    if (myid.eq.rootid) then
       if (index(file_zeocoord,'.cif').gt.0) then
          call readCIF(file_zeocoord,zeo,lunitcell,ztype,zcell,zunit,lprint)
       else if (index(file_zeocoord,'.pdb').gt.0) then
          call readPDB(file_zeocoord,zeo,lunitcell,ztype,zcell,zunit,lprint)
       else if (index(file_zeocoord,'.cssr').gt.0) then
          call readCSSR(file_zeocoord,zeo,lunitcell,ztype,zcell,zunit,lprint)
       end if
    end if

    call mp_bcast(zunit%dup,3,rootid,groupid)
    call mp_bcast(zunit%boxl,3,rootid,groupid)
    call mp_bcast(zunit%ngrid,3,rootid,groupid)
    call mp_bcast(zunit%hmat,9,rootid,groupid)
    call mp_bcast(zunit%hmati,9,rootid,groupid)
    call mp_bcast(ztype%ntype,1,rootid,groupid)
    call mp_bcast(zeo%nbead,1,rootid,groupid)
    if (myid.ne.rootid) then
       allocate(ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype),ztype%num(ztype%ntype)&
        ,zeo%bead(1:zeo%nbead),lunitcell(1:zeo%nbead),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'zeocoord: allocation failed',-1)
    end if
    call mp_bcast(ztype%name,rootid,groupid)
    call mp_bcast(ztype%type,ztype%ntype,rootid,groupid)
    call mp_bcast(ztype%num,ztype%ntype,rootid,groupid)
    call mp_bcast(ztype%radiisq,ztype%ntype,rootid,groupid)
    call mp_bcast(cell_vol(boxZeo),1,rootid,groupid)
    call mp_bcast(cell_ang(boxZeo,:),3,rootid,groupid)
    call mp_bcast(min_width(boxZeo,:),3,rootid,groupid)
    call mp_bcast(hmat(boxZeo,:),9,rootid,groupid)
    call mp_bcast(hmati(boxZeo,:),9,rootid,groupid)
    call mp_bcast_molecule(zeo,rootid,groupid)
    call mp_bcast(lunitcell,zeo%nbead,rootid,groupid)

    if (lprint.and.len_trim(file_supercell).gt.0) call writePDB(file_supercell,zeo,ztype,zcell)
  end subroutine zeocoord

  subroutine addAllBeadTypes()
    integer::imol,iunit,igtype,i,idi,idj,ntij,jerr,list(nntype),mstart,mend
    real::sig6

    ! find all bead types and store them in an array
    if (ltestztb.or.lpore_volume) then
       nmolty=nmolty+1
       nunit(nmolty)=1
       ntype(nmolty,1)=volume_probe
    end if

    if (lsurface_area) then
       nmolty=nmolty+1
       nunit(nmolty)=1
       ntype(nmolty,1)=area_probe
    end if



    zpot%ntype=0 ! number of bead types present in all molec in fort.4
    do imol=1,nmolty
       do iunit=1,nunit(imol)
          ! check if this bead type is already accounted for
          do igtype=1,zpot%ntype
             if (list(igtype).eq.ntype(imol,iunit)) exit
          end do
          ! add to list if not accounted for
          if (igtype.gt.zpot%ntype) then
             zpot%ntype=igtype
             list(igtype)=ntype(imol,iunit)
          end if
       end do
    end do

    if (ltestztb.or.lpore_volume) nmolty=nmolty-1
    if (lsurface_area) nmolty=nmolty-1

    mend=num_points_interpolation/2
    mstart=mend+1-num_points_interpolation
    if (allocated(egrid)) deallocate(egrid,zpot%param,zpot%table,yjtmp,yktmp,yltmp,xt,yt,zt,stat=jerr)
    allocate(egrid(0:product(zunit%ngrid(1:3))-1,zpot%ntype),zpot%param(3,ztype%ntype,zpot%ntype)&
     ,zpot%table(zpot%ntype),yjtmp(mstart:mend),yktmp(mstart:mend),yltmp(mstart:mend),xt(mstart:mend),yt(mstart:mend)&
     ,zt(mstart:mend),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'addAllBeadTypes: allocation failed',myid+1)

    ! loop over all bead types in fort.4
    do igtype=1,zpot%ntype
       idi=list(igtype) ! get the bead index
       zpot%table(igtype)=idi
       if (lij(idi).or.lqchg(idi)) then
          do i=1,ztype%ntype
             idj=ztype%type(i)
             if (lij(idi).and.lij(idj)) then
                ntij = (idi-1)*nntype + idj
                ! LJ parameters are scaled to reduce round-off errors, see exzeof()
                sig6=vvdW(3,ntij)**3/LJScaling
                zpot%param(2,i,igtype)=vvdW(1,ntij)*sig6
                zpot%param(1,i,igtype)=zpot%param(2,i,igtype)*sig6
             else
                zpot%param(1:2,i,igtype)=0._dp
             end if
             if (lqchg(idi).and.lqchg(idj)) then
                zpot%param(3,i,igtype)=qelect(idi)*qelect(idj)
             else
                zpot%param(3,i,igtype)=0._dp
             end if
          end do
       end if
    end do
  end subroutine addAllBeadTypes

!DEC$ ATTRIBUTES FORCEINLINE :: idx
  integer function idx(i,j,k)
    integer,intent(in)::i,j,k

    idx=(k*zunit%ngrid(2)+j)*zunit%ngrid(1)+i
  end function idx

  subroutine combine_energies(zgrid,egrid,n_points_per_piece,ltailcZeotmp,rcuttmp)
    integer,intent(in)::n_points_per_piece
    real,intent(in)::zgrid(1:3,1:ztype%ntype,0:n_points_per_piece-1)
    real,intent(out)::egrid(0:n_points_per_piece-1,1:zpot%ntype)
    real,intent(in)::rcuttmp
    logical,intent(in)::ltailcZeotmp
    integer::i,j,k,idi,idj
    real::rci3,rci9,rho

    if (ltailczeo.and..not.ltailcZeotmp) then
       rci3=1._dp/rcuttmp**3
       rci9=LJScaling*LJScaling*rci3**3
       rci3=LJScaling*rci3
    end if

    egrid=0._dp
    forall(k=0:n_points_per_piece-1,zgrid(1,1,k).ge.overlapValue) egrid(k,:)=overlapValue

    do j=1,ztype%ntype
       idj=ztype%type(j)
       if (lij(idj).or.lqchg(idj)) then
          rho=ztype%num(j)/zcell%vol
          do i=1,zpot%ntype
             idi=zpot%table(i)
             if (lij(idj).and.lij(idi)) then
                where(egrid(:,i).lt.overlapValue) egrid(:,i)=egrid(:,i)+zgrid(1,j,:)*zpot%param(1,j,i)&
                 -zgrid(2,j,:)*zpot%param(2,j,i)
                if (ltailczeo.and..not.ltailcZeotmp) then
                   where(egrid(:,i).lt.overlapValue) egrid(:,i)=egrid(:,i)+twopi*rho*(rci9*zpot%param(1,j,i)/9_dp&
                    -rci3*zpot%param(2,j,i)/3_dp)
                end if
             end if
             if (lqchg(idj).and.lqchg(idi)) then
                where(egrid(:,i).lt.overlapValue) egrid(:,i)=egrid(:,i)&
                 +qqfact*qelect(idi)*qelect(idj)*zgrid(3,j,:)
             end if
          end do
       end if
    end do
  end subroutine combine_energies

  subroutine write_energies(egrid)
    real,intent(in) :: egrid(0:,1:)
    real :: ri(3), scoord(3)
    character(LEN=128)::filename_egrid
    integer io_egrid, nynx, nznynx
    integer i, j, k, r, ii


    nynx = product(zunit%ngrid(1:2))
    nznynx = nynx*zunit%ngrid(3)
    write(io_output,"(A)") 'Writing energy grids in zeolite to files for each bead'

    do ii=1,zpot%ntype
        io_egrid = get_iounit()
        write(filename_egrid,"(I2,A,A)") atoms%list(zpot%table(ii)),'_', chemid(zpot%table(ii))
        filename_egrid= adjustl(filename_egrid)
        filename_egrid='energy_grid_'//trim(filename_egrid)//'.out'
        open(unit=io_egrid,file=filename_egrid,status='replace')
        write(io_egrid,*) zunit%ngrid(:)
        do r = 0, nznynx - 1
            k=r/(nynx)
            j=(r-k*nynx)/zunit%ngrid(1)
            i=r-k*nynx-j*zunit%ngrid(1)
            scoord(1) = (real(i,dp)/zunit%ngrid(1))/zunit%dup(1)
            scoord(2) = (real(j,dp)/zunit%ngrid(2))/zunit%dup(2)
            scoord(3) = (real(k,dp)/zunit%ngrid(3))/zunit%dup(3)
            ri(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
            ri(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
            ri(3)=scoord(3)*zcell%hmat(3,3)%val
            write(io_egrid,"(3I4,3F8.3,E11.3)")  i, j, k, ri, egrid(r, ii)
        enddo
    close(io_egrid)
    enddo

  end subroutine write_energies

  subroutine suzeo(lprint)
    use util_search,only:indexOf
    use util_string,only:format_n
    use util_mp,only:mp_bcast,mp_set_displs,mp_allgather
    logical,intent(in)::lprint
    character(LEN=128)::atom
    integer::rcounts(numprocs),displs(numprocs),io_ztb,jerr,i,j,k,st,st0,p,r,n_points_per_piece,my_start,my_end,blocksize&
     ,nznynx,nynx,ngridtmp(3),zntypetmp
    real::wzeo,zunittmp(3),zangtmp(3),rcuttmp
    logical::lewaldtmp,ltailczeotmp,lshifttmp,is_ztb_complete

    do i=1,ztype%ntype
       ztype%type(i)=indexOf(atoms,ztype%type(i))
       if (ztype%type(i).eq.0) then
          call err_exit(__FILE__,__LINE__,'zeocoord: type '//integer_to_string(i)//' undefined: '//ztype%name(i),myid+1)
       end if
    end do

    ! Calculate zeolite density
    wzeo=dot_product(ztype%num(1:ztype%ntype),mass(ztype%type(1:ztype%ntype)))
    if (lprint) then
       write(io_output,"(/,' Tabulated framework potential: ',/&
                       &,' --------------------------------------------------',/,"&
                       //format_n(ztype%ntype,"(' number of ',A,':',I0,3X)")//&
                       ",/,' number of framework atom types   = ',i10,/&
                       &,' number of framework atoms        = ',i10,/&
                       &,' number of atoms in the unit cell = ',i10,/&
                       &,' framework mass                   = ',e16.9,' grams',/&
                       &,' framework volume                 = ', f16.5,' Angst**3',/&
                       &,' one adsorbed molecule in sim box = ',e16.9 ,' mol/kg',/&
                       &,' unit-cell size: ',f7.4,' x ',f7.4,' x ',f7.4,/&
                       &,'   x-dir       : ',i5,'  size: ',f7.4,/&
                       &,'   y-dir       : ',i5,'  size: ',f7.4,/&
                       &,'   z-dir       : ',i5,'  size: ',f7.4,/)")&
                       (trim(ztype%name(i)),ztype%num(i),i=1,ztype%ntype)&
                       ,ztype%ntype,zeo%nbead,count(lunitcell)&
                       ,wzeo/N_Avogadro,zcell%vol,1000.0_dp/wzeo&
                       ,zunit%boxl(1),zunit%boxl(2),zunit%boxl(3)&
                       ,zunit%ngrid(1),zunit%boxl(1)/zunit%ngrid(1)&
                       ,zunit%ngrid(2),zunit%boxl(2)/zunit%ngrid(2)&
                       ,zunit%ngrid(3),zunit%boxl(3)/zunit%ngrid(3)
    end if

    if (lsurface_area) call zsurface()

    ! tabulation of the zeolite potential
    call addAllBeadTypes()
    call setpbc(boxZeo)

    nlayermax=0
    nynx=product(zunit%ngrid(1:2))
    nznynx=nynx*zunit%ngrid(3)
    n_points_per_piece=ceiling(real(nznynx,dp)/n_pieces_ztb,dp)
    allocate(my_zgrid(3,ztype%ntype,0:ceiling(real(n_points_per_piece,dp)/numprocs,dp)-1),zgrid(3,ztype%ntype,0:n_points_per_piece-1),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'suzeo: allocation failed',myid+1)

    if (myid.eq.rootid) then
       st=0
       lewaldtmp=lewald
       ltailczeotmp=ltailcZeo
       lshifttmp=lshift
       rcuttmp=zcell%cut

       io_ztb=get_iounit()
       open(unit=io_ztb,access='stream',action='read',file=file_ztb,form='unformatted',iostat=jerr,status='old')
       if (jerr.eq.0) then
          ! read zeolite table from disk
          is_ztb_complete=.false.
          if (lprint) write(io_output,*) 'read in tabulated potential:'
          read(io_ztb,end=100) zunittmp,zangtmp,ngridtmp,zntypetmp,lewaldtmp,ltailczeotmp,lshifttmp,rcuttmp
          if (ANY(abs(zunittmp-zunit%boxl).gt.eps).or.ANY(ngridtmp.ne.zunit%ngrid).or.(zntypetmp.ne.ztype%ntype))&
           call err_exit(__FILE__,__LINE__,'problem 1 in zeolite potential table',myid+1)
          do i=1,ztype%ntype
             read(io_ztb,end=100) atom
             if (trim(ztype%name(i))/=trim(atom)) then
                write(io_output,*) i,' atom should be ',trim(atom)
                call err_exit(__FILE__,__LINE__,'problem 2 in zeolite potential table',myid+1)
             end if
          end do
          do p=0,n_pieces_ztb-1
             n_points_per_piece=ceiling(real(nznynx-st,dp)/(n_pieces_ztb-p),dp)
             read(io_ztb,end=100) zgrid(1:3,1:ztype%ntype,0:n_points_per_piece-1)
             call combine_energies(zgrid(1:3,1:ztype%ntype,0:n_points_per_piece-1),egrid(st:st+n_points_per_piece-1,:),n_points_per_piece,ltailcZeotmp,rcuttmp)
             st=st+n_points_per_piece
          end do
          is_ztb_complete=.true.

100       if (.not.is_ztb_complete) then
             jerr=1
             close(io_ztb)
             write(io_output,*) st,' points read in.'
          end if
       end if
    end if

    call mp_bcast(jerr,1,rootid,groupid)
    call mp_bcast(ltailcZeotmp,1,rootid,groupid)
    call mp_bcast(rcuttmp,1,rootid,groupid)
    call mp_bcast(egrid(:,:),product(zunit%ngrid(1:3))*zpot%ntype,rootid,groupid)

    if (jerr.ne.0) then
       call mp_bcast(st,1,rootid,groupid)
       call mp_bcast(lewaldtmp,1,rootid,groupid)
       call mp_bcast(lshifttmp,1,rootid,groupid)
       ! make a tabulated potential of the zeolite
       if (myid.eq.rootid) then
          write(io_output,*) 'make tabulated potential'
          open(unit=io_ztb,access='stream',action='readwrite',file=file_ztb,form='unformatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) then
             write(*,*) 'error in trying to create tabulated potential: iostat =', jerr
             call err_exit(__FILE__,__LINE__,'cannot create file for tabulated potential',myid+1)
          end if
          if (st.eq.0) then
             write(io_ztb) zunit%boxl,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val,zunit%ngrid,ztype%ntype,lewald,ltailcZeo&
              ,lshift,zcell%cut
             do i=1,ztype%ntype
                write(io_ztb) ztype%name(i)
             end do
          else
             read(io_ztb) zunittmp,zangtmp,ngridtmp,zntypetmp,lewaldtmp,ltailczeotmp,lshifttmp,rcuttmp
             do i=1,ztype%ntype
                read(io_ztb) atom
             end do
          end if
          write(io_output,*) 'time 1:',time_now()
       end if

       st0=st
       st=0
       do p=0,n_pieces_ztb-1
          n_points_per_piece=ceiling(real(nznynx-st,dp)/(n_pieces_ztb-p),dp)
          if (st.lt.st0) then
             call mp_bcast(egrid(st:st+n_points_per_piece-1,:),zpot%ntype*n_points_per_piece,rootid,groupid)
             if (myid.eq.rootid) read(io_ztb) zgrid(:,:,0:n_points_per_piece-1)
             st=st+n_points_per_piece
             cycle
          end if
          blocksize = n_points_per_piece / numprocs
          rcounts = blocksize
          blocksize = n_points_per_piece - blocksize*numprocs
          if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
          call mp_set_displs(rcounts,displs,blocksize,numprocs)
          my_start = displs(myid+1) + st
          my_end = my_start + rcounts(myid+1) - 1

!$omp parallel &
!$omp private(r,i,j,k) shared(my_start,my_end,my_zgrid) default(shared)
!$omp do
          do r=my_start,my_end
             k=r/(nynx)
             j=(r-k*nynx)/zunit%ngrid(1)
             i=r-k*nynx-j*zunit%ngrid(1)
             ! pass to exzeof arguments in fractional coordinates with respect to the unit cell
             call exzeof(my_zgrid(:,:,r-my_start),real(i,dp)/zunit%ngrid(1),real(j,dp)/zunit%ngrid(2),real(k,dp)/zunit%ngrid(3))
          end do
!$omp end parallel

          call mp_allgather(my_zgrid,zgrid(:,:,0:n_points_per_piece-1),rcounts,displs,groupid)
          if (myid.eq.rootid) write(io_ztb) zgrid(:,:,0:n_points_per_piece-1)
          call combine_energies(zgrid(1:3,1:ztype%ntype,0:n_points_per_piece-1),egrid(st:st+n_points_per_piece-1,:),n_points_per_piece,ltailcZeotmp,rcuttmp)
          st=st+n_points_per_piece
       end do

       if (myid.eq.rootid) then
          write(io_output,*) 'time 2:',time_now()
          if (ltailcZeo) write(io_output,*) 'maxlayer = ',nlayermax
       end if
    end if

    if (ltestztb.or.lpore_volume) call ztest()
    if (printztb) call write_energies(egrid(:,:))

    deallocate(lunitcell,my_zgrid,zgrid,zeo%bead,ztype%type,ztype%num,ztype%radiisq,ztype%name,zpot%param)
    if (lprint) then
       close(io_ztb)
       write(io_output,'(4(A,L1),A,G10.3,A)') 'lewald[',lewaldtmp,'] ltailcZeo[',ltailcZeotmp,'] lshift['&
        ,lshifttmp,'] AddTailc[',ltailcZeo,'] rcut[',rcuttmp,']'
    end if

  end subroutine suzeo

!> \param i,j,k are in fractional coordinates with respect to the unit cell
!>
!> \remarks U_LJ = A/r^12 + B/r^6, scale A by 4*10^8, B by 2*10^4 to reduce round-off error
  subroutine exzeof(tab,i,j,k)
    real,intent(out)::tab(3,ztype%ntype)
    real,intent(in)::i,j,k

    integer::izeo,layer,ii,jj,kk,iztype
    real::rcutsq,vac,vbc,vb,r,scoord(3),ri(3),dr(3),r2,vnew(2,ztype%ntype)

    tab=0._dp
    rcutsq = zcell%cut*zcell%cut
    vbc=LJScaling/rcutsq**3
    vac=vbc*vbc

    if (lewald.or.(.not.ltailcZeo)) then
       ! further scale i,j,k with respect to the whole cell
       !!!
       scoord(1)=real(i,dp)/zunit%dup(1)
       scoord(2)=real(j,dp)/zunit%dup(2)
       scoord(3)=real(k,dp)/zunit%dup(3)
       ri(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
       ri(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
       ri(3)=scoord(3)*zcell%hmat(3,3)%val
       do izeo=1,zeo%nbead
          iztype=zeo%bead(izeo)%type
          dr=ri-zeo%bead(izeo)%coord
          call mimage(dr(1),dr(2),dr(3),boxZeo)
          r2=dot_product(dr,dr)
          if (r2.le.ztype%radiisq(iztype)) then
             tab=overlapValue
             return
          else if (r2 .lt. rcutsq) then
             if (.not.ltailcZeo) then
                vb=LJScaling/r2**3
                if (lshift) then
                   tab(2,iztype)=tab(2,iztype)+vb-vbc
                   tab(1,iztype)=tab(1,iztype)+vb*vb-vac
                else
                   tab(2,iztype)=tab(2,iztype)+vb
                   tab(1,iztype)=tab(1,iztype)+vb*vb
                end if
             end if
             r=sqrt(r2)
             if (lewald) then
                tab(3,iztype)=tab(3,iztype)+erfunc(calp(boxZeo)*r)/r
             else
                tab(3,iztype)=tab(3,iztype)+1._dp/r
             end if
          end if
       end do

       if (lewald) call recipzeo(tab(3,:),ri)

    end if

    ! Calculate the Lennard-Jones interactions, include as many layers
    ! of neighboring unit cells as needed for the specified precision
    if (ltailcZeo) then
       layer=0
       do
          vnew=0._dp
          do izeo=1,zeo%nbead
             if (lunitcell(izeo)) then
                iztype=zeo%bead(izeo)%type
                do ii=-layer,layer
                   do jj=-layer,layer
                      do kk=-layer,layer
                         if (abs(ii).eq.layer .or. abs(jj).eq.layer .or.abs(kk).eq.layer) then
                            !!!
                            scoord(1) = zeo%bead(izeo)%coord(1)*zcell%hmati(1,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(1,2)%val&
                             +zeo%bead(izeo)%coord(3)*zcell%hmati(1,3)%val+real(ii-i,dp)/zunit%dup(1)
                            scoord(2) = zeo%bead(izeo)%coord(1)*zcell%hmati(2,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(2,2)%val&
                             +zeo%bead(izeo)%coord(3)*zcell%hmati(2,3)%val+real(jj-j,dp)/zunit%dup(2)
                            scoord(3) = zeo%bead(izeo)%coord(1)*zcell%hmati(3,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(3,2)%val&
                             +zeo%bead(izeo)%coord(3)*zcell%hmati(3,3)%val+real(kk-k,dp)/zunit%dup(3)
                            dr(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
                            dr(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
                            dr(3)=scoord(3)*zcell%hmat(3,3)%val
                            r2=dot_product(dr,dr)
                            if (r2.le.ztype%radiisq(iztype)) then
                               tab=overlapValue
                               return
                            end if
                            vb=LJScaling/r2**3
                            if (lshift) then
                               vnew(2,iztype)=vnew(2,iztype)+vb-vbc
                               vnew(1,iztype)=vnew(1,iztype)+vb*vb-vac
                            else
                               vnew(2,iztype)=vnew(2,iztype)+vb
                               vnew(1,iztype)=vnew(1,iztype)+vb*vb
                            end if
                         end if
                      end do
                   end do
                end do
             end if
          end do
          tab(1:2,:)=tab(1:2,:)+vnew
          layer=layer+1
          if (layer.gt.nlayermax) nlayermax=layer
          if (abs(sum(vnew(1,:)-vnew(2,:))).lt.requiredPrecision) exit
       end do
    end if

  end subroutine exzeof

  subroutine recipzeo(tab,ri)
    real,intent(out)::tab(ztype%ntype)
    real,intent(in)::ri(3)

    integer::kmax(3),i,j,l,m,n,kmin(2:3)
    real::hmatik(3,3),ki(3),alpsqr4,hmaxsq,ksqr,arg,sums(ztype%ntype),vrecipz(ztype%ntype)

    ! *** Set up the reciprocal space vectors ***
    vrecipz = 0.0E+0_dp

    forall(i=1:3,j=1:3) hmatik(j,i) = twopi*zcell%hmati(j,i)%val
    kmax(1) = int(zcell%hmat(1,1)%val*zcell%calp)+1
    kmax(2) = int(zcell%hmat(2,2)%val*zcell%calp)+1
    kmax(3) = int(zcell%hmat(3,3)%val*zcell%calp)+1
    if (zcell%solid.and..not.zcell%ortho) kmax = kmax+1

    alpsqr4 = 4.0E0_dp*zcell%calp*zcell%calp
    hmaxsq = alpsqr4*onepi*onepi

    ! Generate the reciprocal-space
    ! here -kmax(1),-kmax(1)+1,...,-1 are skipped, so no need to divide by 2 for the prefactor
    do l = 0,kmax(1)
       if ( l .eq. 0 ) then
          kmin(2) = 0
       else
          kmin(2) = -kmax(2)
       end if
       do m = kmin(2), kmax(2)
          if (l .eq. 0 .and. m .eq. 0) then
             kmin(3) = 1
          else
             kmin(3) = -kmax(3)
          end if
          do n = kmin(3), kmax(3)
             ki(1) = real(l,dp)*hmatik(1,1)+real(m,dp)*hmatik(2,1)+real(n,dp)*hmatik(3,1)
             ki(2) = real(l,dp)*hmatik(1,2)+real(m,dp)*hmatik(2,2)+real(n,dp)*hmatik(3,2)
             ki(3) = real(l,dp)*hmatik(1,3)+real(m,dp)*hmatik(2,3)+real(n,dp)*hmatik(3,3)
             ksqr = dot_product(ki,ki)
             ! if ( ksqr .lt. hmaxsq ) then
             ! sometimes these are about equal, which can cause different
             ! behavior on 32 and 64 bit machines without this .and. statement
             if ( hmaxsq-ksqr.gt.eps ) then
                sums = 0.0E0_dp
                do i = 1,zeo%nbead
                   arg=dot_product(ki,ri-zeo%bead(i)%coord)
                   sums(zeo%bead(i)%type)=sums(zeo%bead(i)%type)+cos(arg)
                end do
                vrecipz=vrecipz+sums*exp(-ksqr/alpsqr4)/ksqr
             end if
          end do
       end do
    end do

    tab=tab+vrecipz*8._dp*onepi/zcell%vol

  end subroutine recipzeo

  function exzeo(xi,yi,zi,idi,ignoreTable)
    use util_math,only:polint
    real::exzeo
    real,intent(in)::xi,yi,zi
    integer,intent(in)::idi
    logical,intent(in)::ignoreTable

    logical::lignore
    integer::mstart,mend,j,j0,jp,k,k0,kp,l,l0,lp,igtype,idj
    real::scoord(3),r(3),tab(3,ztype%ntype),rci3,rci9,rho

    ! find the correct bead type
    do igtype=1,zpot%ntype
       if (zpot%table(igtype).eq.idi) exit
    end do
    if (igtype.gt.zpot%ntype) then
       call err_exit(__FILE__,__LINE__,'exzeo: no such bead type',myid+1)
    end if

    lignore=ignoreTable

    ! fold coordinates into the unit cell, result in fractional coordinates
    !!!
    scoord(1)=(xi*zunit%hmati(1,1)+yi*zunit%hmati(1,2)+zi*zunit%hmati(1,3))
    scoord(2)=(xi*zunit%hmati(2,1)+yi*zunit%hmati(2,2)+zi*zunit%hmati(2,3))
    scoord(3)=(xi*zunit%hmati(3,1)+yi*zunit%hmati(3,2)+zi*zunit%hmati(3,3))
    scoord(1)=scoord(1)-floor(scoord(1))
    scoord(2)=scoord(2)-floor(scoord(2))
    scoord(3)=scoord(3)-floor(scoord(3))
    ! get the Cartesian coordinates of the point in the unit cell
    r(1)=scoord(1)*zunit%hmat(1,1)+scoord(2)*zunit%hmat(1,2)+scoord(3)*zunit%hmat(1,3)
    r(2)=scoord(2)*zunit%hmat(2,2)+scoord(3)*zunit%hmat(2,3)
    r(3)=scoord(3)*zunit%hmat(3,3)
    ! get the index of the grid that the point resides in
    j = scoord(1)*zunit%ngrid(1)
    k = scoord(2)*zunit%ngrid(2)
    l = scoord(3)*zunit%ngrid(3)

    exzeo=upperLimit
    if (.not.lignore) then
       ! calculation using a grid
       mend=num_points_interpolation/2
       mstart=mend+1-num_points_interpolation

       ! test if in the reasonable regime
       if (egrid(idx(j,k,l),igtype).ge.upperLimit) return
       ! block centered around: j,k,l
       ! set up hulp array: (allow for going beyond unit cell
       ! for polynom fitting)
       do l0=mstart,mend
          lp=l+l0
          scoord(3)=real(lp,dp)/zunit%ngrid(3)
          ! ---    store x,y,z values around xi,yi,zi in arrays
          zt(l0)=scoord(3)*zunit%hmat(3,3)
          if (lp.lt.0) lp=lp+zunit%ngrid(3)
          if (lp.ge.zunit%ngrid(3)) lp=lp-zunit%ngrid(3)
          do k0=mstart,mend
             kp=k+k0
             scoord(2)=real(kp,dp)/zunit%ngrid(2)
             yt(k0)=scoord(2)*zunit%hmat(2,2)+scoord(3)*zunit%hmat(2,3)
             if (kp.lt.0) kp=kp+zunit%ngrid(2)
             if (kp.ge.zunit%ngrid(2)) kp=kp-zunit%ngrid(2)
             do j0=mstart,mend
                jp=j+j0
                scoord(1)=real(jp,dp)/zunit%ngrid(1)
                xt(j0)=scoord(1)*zunit%hmat(1,1)+scoord(2)*zunit%hmat(1,2)+scoord(3)*zunit%hmat(1,3)
                if (jp.lt.0) jp=jp+zunit%ngrid(1)
                if (jp.ge.zunit%ngrid(1)) jp=jp-zunit%ngrid(1)
                yjtmp(j0)=egrid(idx(jp,kp,lp),igtype)
                if (yjtmp(j0).ge.upperLimit) return
             end do
             call polint(xt,yjtmp,num_points_interpolation,r(1),yktmp(k0))
          end do
          call polint(yt,yktmp,num_points_interpolation,r(2),yltmp(l0))
       end do
       call polint(zt,yltmp,num_points_interpolation,r(3),exzeo)
    else
       ! calculating interaction energy with the zeolite framework explicitly
       call exzeof(tab,scoord(1),scoord(2),scoord(3))
       if (tab(1,1).ge.upperLimit) return

       if (ltailc) then
          rci3=1._dp/zcell%cut**3
          rci9=LJScaling*LJScaling*rci3**3
          rci3=LJScaling*rci3
       end if

       exzeo=0.0E0_dp
       do j=1,ztype%ntype
          idj=ztype%type(j)
          if (lij(idj).or.lqchg(idj)) then
             rho=ztype%num(j)/zcell%vol
             if (lij(idj).and.lij(idi)) then
                exzeo=exzeo+tab(1,j)*zpot%param(1,j,igtype)-tab(2,j)*zpot%param(2,j,igtype)
                if (ltailc) then
                   exzeo=exzeo+twopi*rho*(rci9*zpot%param(1,j,igtype)/9_dp-rci3*zpot%param(2,j,igtype)/3_dp)
                end if
             end if
             if (lqchg(idj).and.lqchg(idi)) then
                exzeo=exzeo+qqfact*qelect(idi)*qelect(idj)*tab(3,j)
            end if
          end if
       end do
    end if

  end function exzeo

!> \brief Test accuracy of tabulated zeolite potential and calculate void volume
!>
!> \par Literature
!> O. Talu and A.L. Myers, "Molecular simulation of adsorption: Gibbs dividing surface and comparison with experiment", AICHE J., 47(5), 1160-1168 (2001).
  subroutine ztest()
    use util_random,only:random
    integer::i,tel
    real::errTot,errRel,errAbs,err,BoltTabulated,eBoltTabulated,BoltExplicit,eBoltExplicit,xi,yi,zi,Utabulated,Uexplicit,weight

    ! test accuracy
    if (myid.eq.rootid) write(io_output,'(A,/,A)') ' Test accuracy of tabulated zeolite potential & Calculate void volume'&
     ,' -------------------------------------------------'
    tel=0
    errTot=0
    errRel=0
    errAbs=0
    BoltTabulated=0
    eBoltTabulated=0
    BoltExplicit=0
    eBoltExplicit=0
    do i=1,volume_nsample
       xi=random(-1)*zunit%boxl(1)
       yi=random(-1)*zunit%boxl(2)
       zi=random(-1)*zunit%boxl(3)
       Utabulated=exzeo(xi,yi,zi,volume_probe,ignoreTable=.false.)
       if (ltestztb) then
          Uexplicit=exzeo(xi,yi,zi,volume_probe,ignoreTable=.true.)
       else
          Uexplicit=Utabulated
       end if
       weight=exp(-Uexplicit*beta)
       if (weight.gt.1.0E-5_dp) then
          tel=tel+1
          BoltExplicit=BoltExplicit+weight
          eBoltExplicit=eBoltExplicit+Uexplicit*weight
          weight=exp(-Utabulated*beta)
          BoltTabulated=BoltTabulated+weight
          eBoltTabulated=eBoltTabulated+Utabulated*weight
          err=abs(Uexplicit-Utabulated)
          if (errAbs.lt.err) errAbs=err
          err=abs(err/Uexplicit)
          if (errRel.lt.err) errRel=err
          errTot=errTot+err
          if (myid.eq.rootid .and. err.gt.3.0E-2_dp) then
             write(io_output,'("WARNING: interpolation error at (",3(F8.5,1X),")=",G20.7," > 3%")') xi,yi,zi,err
             write(io_output,*) Utabulated,Uexplicit
          end if
       end if
    end do

    if (myid.eq.rootid) then
       write(io_output,'(A,I0,A,I0,A)') ' test over ',tel,' out of ',volume_nsample,' random positions '
       write(io_output,'(A,G16.9)') ' average error: ',errTot/tel
       write(io_output,'(A,G16.9)') ' maximum relative error: ',errRel
       write(io_output,'(A,G16.9)') ' maximum absolute error [K]: ',errAbs
       write(io_output,'(A,G16.9,A,G16.9,A)') ' void fraction: ',BoltTabulated/volume_nsample, '(tabulated), '&
        ,BoltExplicit/volume_nsample,'(explicit)'
       write(io_output,'(A,G16.9,A,G16.9,A)') ' void volume in [Angstrom^3]: ',BoltTabulated*zcell%vol/volume_nsample&
        , '(tabulated), ',BoltExplicit*zcell%vol/volume_nsample,'(explicit)'
       write(io_output,'(A,G16.9,A,G16.9,A)') ' void volume in [cm^3/g]: '&
        ,BoltTabulated*zcell%vol/volume_nsample*N_Avogadro*1E-24_dp/dot_product(ztype%num(1:ztype%ntype)&
        ,mass(ztype%type(1:ztype%ntype))), '(tabulated), '&
        ,BoltExplicit*zcell%vol/volume_nsample*N_Avogadro*1E-24_dp/dot_product(ztype%num(1:ztype%ntype)&
        ,mass(ztype%type(1:ztype%ntype))),'(explicit)'
       write(io_output,'(A,G16.9,A,G16.9,A)') ' Boltzmann averaged energy in [K]: ',eBoltTabulated/BoltTabulated,'(tabulated), '&
        ,eBoltExplicit/BoltExplicit,'(explicit)'
       write(io_output,*)
    end if

    return
  end subroutine ztest

!> \brief Calculate geometric surface area
!>
!> \par Literature
!> 1. O.K. Farha, A.O. Yazaydin, I. Eryazici, C.D. Malliakas, B.G. Hauser, M.G. Kanatzidis, S.T. Nguyen, R.Q. Snurr, and J.T. Hupp,
!>   "De novo synthesis of a metal-organic framework material featuring ultra-high surface area and gas storage capacities", Nature Chem., xx(x), xxx-xxx (2010). \n
!> 2. T. Duren, L. Sarkisov, O.M. Yaghi, and R.Q. Snurr, "Design of New Materials for Methane Storage", Langmuir, 20(7), 2683-2689 (2004). \n
!> 3. K.S. Walton and R.Q. Snurr, "Applicability of the BET method for determining surface areas of microporous metal-organic frameworks", J. Am. Chem. Soc., 129(27), 8552-8556 (2007). \n
!> 4. T. Duren, F. Millange, G. Ferey, K.S. Walton, and R.Q. Snurr, "Calculating geometric surface areas as a characterization tool for metal-organic frameworks", J. Phys. Chem., 111(42), 15350-15356 (2007). \n
!> 5. T. Duren, Y.S. Bae, and R.Q. Snurr, "Using molecular simulation to characterise metal-organic frameworks for adsorption applications", Chem. Soc. Rev., 38(5), 1237-1247 (2009). \n
!> 6. Y.S. Bae, A.O. Yazaydin, and R.Q. Snurr, "Evaluation of the BET method for determining surface areas of MOFs and zeolites that contain ultra-micropores", Langmuir, 26(8), 5475-5483 (2010). \n
  subroutine zsurface()
    use util_random,only:random,sphere
    integer::i,j,k,ntij,ncount
    real::stotal,sigij,rsq,scoord(3),coord(3),sdr(3),dr(3)
    real,allocatable::position(:,:)

    Write(io_output,'(A,/,A)') 'Calculate geometric accessible surface area',' -------------------------------------------------'

    allocate(position(3,zeo%nbead))

    DO i=1, zeo%nbead
       if (.not.lunitcell(i)) cycle
       call absoluteToFractional(scoord,zeo%bead(i)%coord,zcell)
       position(:,i)=scoord*zunit%dup
       if (ANY(position(:,i).ge.1)) then
          write(io_output,*) i,position(:,i),scoord,zunit%dup
          call err_exit(__FILE__,__LINE__,'',myid+1)
       end if
    END DO

    stotal=0.0E0_dp
    Do i=1,zeo%nbead ! Loop over all framework atoms
       if (.not.lunitcell(i)) cycle
       ntij = (area_probe-1)*nntype + ztype%type(zeo%bead(i)%type)
       sigij=vvdW(2,ntij)

       ncount=0
       Do j=1,area_nsample ! Number of trial positions around each framework atom
          call sphere(coord(1),coord(2),coord(3),-1)

          ! Make this vector of length (sigma_atom+sigma_probe)/2.0 and centered at particle i
          coord=coord*sigij+zeo%bead(i)%coord
          CALL foldToUnitCell(coord,zunit,scoord)

          ! Check for overlap
          Do k=1,zeo%nbead
             if(.not.lunitcell(k).or.k.eq.i) cycle
             sdr = scoord - position(:,k)
             sdr = sdr - int(2.0*sdr) ! apply PBC

             call fractionalToAbsolute(dr,sdr/zunit%dup,zcell)
             rsq=dot_product(dr,dr)
             ntij = (area_probe-1)*nntype + ztype%type(zeo%bead(k)%type)
             If(rsq.lt.0.998001E0_dp*vvdW(3,ntij)) exit
          End Do

          If (k.le.zeo%nbead) cycle
          ncount=ncount+1
       End Do

       ! Surface area for sphere i in real units [Angstrom^2]
       stotal=stotal+sigij**2*real(ncount,dp)/real(area_nsample,dp)
    End Do

    ! Report results
    stotal=stotal*fourpi*product(zunit%dup)
    Write(io_output,'(A,F12.2)') ' Total surface area in [Angstroms^2]: ', stotal
    Write(io_output,'(A,F12.2)') ' Total surface area per volume in [m^2/cm^3]: ',stotal/zcell%vol*1.0E4_dp
    Write(io_output,'(A,F12.2,/)') ' Total surface area per volume in [m^2/g]: '&
     ,stotal*N_Avogadro*1E-20_dp/dot_product(ztype%num(1:ztype%ntype),mass(ztype%type(1:ztype%ntype)))

    deallocate(position)

  end subroutine zsurface
end module zeolite
