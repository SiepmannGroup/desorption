! *****************************************************************************
!> \brief Handles CSSR (SERC Daresbury Laboratory's Cambridge Structure Search and Retrieval) files
!>
!> \note Bonding connectivities and charges are ignored.
!> \verbatim
!> LINE      FORMAT                           DEFINITION
!> -------------------------------------------------------------------------------------
!>  1    38X,3F8.3                            Cell dimensions of the a, b, c vectors in Angstrom.
!>  2    21X,3F8.3,4X,'SPGR =',I3,1X,A11      Cell angles of alpha, beta, gamma in degree, space group number, space group name.
!>  3    2I4,1X,A60                           Number of atoms, coordinate system flag (0=fractional coordinates, 1=orthogonal coordinates in Angstrom), brief title/description.
!>  4    A53                                  Long title/description.
!>  5    I4,1X,A4,2X,3(F9.5,1X),8I4,1X,F7.3   Atom ID, atom name(element symbol optionally followed by its ID), x, y, z coordinates, bonding connectivities (max 8), charge.
!> \endverbatim
! *****************************************************************************
MODULE parser_cssr
  use var_type,only:default_string_length
  use util_runtime,only:err_exit
  use util_files,only:get_iounit,readLine
  use sim_cell
  use sim_particle
  use sim_zeolite,only:ZeoliteUnitCellGridType,ZeoliteBeadType,setUpAtom,setUpCellStruct,foldToCenterCell,fractionalToAbsolute
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: readCSSR

CONTAINS
  SUBROUTINE readCSSR(fileCSSR,zeo,lunitcell,ztype,zcell,zunit,lprint)
    character(LEN=*),intent(in)::fileCSSR
    type(MoleculeType),intent(out)::zeo
    logical,allocatable,intent(out)::lunitcell(:)
    type(ZeoliteBeadType),intent(out)::ztype
    type(CellMaskType),intent(out)::zcell !< the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(out)::zunit
    LOGICAL,INTENT(IN)::lprint

    integer::IOCSSR,jerr,i,cflag,pos
    character(LEN=default_string_length)::line,atom
    real::coord(3)

    IOCSSR=get_iounit()
    open(unit=IOCSSR,access='sequential',action='read',file=fileCSSR,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open zeolite CSSR file',-1)
    end if

    CALL readLine(IOCSSR,line,.false.,jerr)
    IF(jerr.ne.0) call err_exit(__FILE__,__LINE__,'wrong CSSR file format',-1)

    read(line(2:),*) zunit%dup(1),zunit%dup(2),zunit%dup(3),ztype%ntype
    allocate(ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype),ztype%num(ztype%ntype),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readCSSR: ztype allocation failed',-1)

    do i=1,ztype%ntype
       CALL readLine(IOCSSR,line,.false.,jerr)
       IF(jerr.ne.0) call err_exit(__FILE__,__LINE__,'wrong CSSR file format',-1)
       read(line(2:),*) ztype%name(i),ztype%type(i),ztype%radiisq(i)
       ztype%radiisq(i)=ztype%radiisq(i)*ztype%radiisq(i)
       ztype%num(i)=0
    end do

    read(IOCSSR,*) zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val
    read(IOCSSR,*) zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val
    call setUpCellStruct(zcell,zunit,lprint,fileCSSR)

    read(IOCSSR,*) zeo%nbead,cflag
    read(IOCSSR,*)
    allocate(zeo%bead(zeo%nbead),lunitcell(zeo%nbead),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readCSSR: zeo%nbead allocation failed',-1)

    do i = 1,zeo%nbead
       ! In order to accommodate number of framework atoms up to five digits, I5 is used instead of I4,1X, which should also work if the latter form is actually used.
       read(IOCSSR,'(i5,a4,2x,3(f9.5,1x))') pos,atom,coord
       if (cflag.eq.0) then
          coord=coord-floor(coord)
          call fractionalToAbsolute(zeo%bead(i)%coord,coord,zcell)
       else
          zeo%bead(i)%coord=coord
          CALL foldToCenterCell(zeo%bead(i)%coord,zcell,coord)
       end if
       call setUpAtom(atom,i,zeo,lunitcell,ztype,zcell,zunit)
    end do

!     do i=1,ztype%ntype
!        if (ztype%num(i).eq.0) then
!           ztype%type(i)=ztype%type(ztype%ntype)
!           ztype%num(i)=ztype%num(ztype%ntype)
!           ztype%radiisq(i)=ztype%radiisq(ztype%ntype)
!           ztype%name(i)=ztype%name(ztype%ntype)
!           ztype%ntype=ztype%ntype-1
!           if (i.ge.ztype%ntype) exit
!        end if
!     end do

    close(IOCSSR)

  END SUBROUTINE readCSSR
END MODULE parser_cssr
