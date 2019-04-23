! *****************************************************************************
!> \brief Handles PDB (Protein Data Bank) files
!>
!> \par Protein Data Bank Contents Guide:
!> Atomic Coordinate Entry Format Version 3.3 (July, 2011)\n
!>   http://www.wwpdb.org/documentation/format33/v3.3.html
!> \verbatim
!> COLUMNS        DATA  TYPE    FIELD        DEFINITION
!> -------------------------------------------------------------------------------------
!>  1 -  6        Record name   "ATOM  "
!>  7 - 11        Integer       serial       Atom serial number.
!> 13 - 16        Atom          name         Atom name.
!> 17             Character     altLoc       Alternate location indicator.
!> 18 - 20        Residue name  resName      Residue name.
!> 22             Character     chainID      Chain identifier.
!> 23 - 26        Integer       resSeq       Residue sequence number.
!> 27             AChar         iCode        Code for insertion of residues.
!> 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
!> 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
!> 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
!> 55 - 60        Real(6.2)     occupancy    Occupancy.
!> 61 - 66        Real(6.2)     tempFactor   Temperature factor.
!> 77 - 78        LString(2)    element      Element symbol, right-justified.
!> 79 - 80        LString(2)    charge       Charge on the atom.
!> -------------------------------------------------------------------------------------
!>  1 -  6        Record name    "CONECT"
!>  7 - 11        Integer        serial       Atom serial number
!> 12 - 16        Integer        serial       Serial number of bonded atom
!> 17 - 21        Integer        serial       Serial number of bonded atom
!> 22 - 26        Integer        serial       Serial number of bonded atom
!> 27 - 31        Integer        serial       Serial number of bonded atom
!> \endverbatim
! *****************************************************************************
MODULE parser_pdb
  use var_type,only:default_string_length
  use util_runtime,only:err_exit
  use util_files,only:get_iounit,readLine
  use sim_cell
  use sim_particle
  use sim_zeolite,only:ZeoliteUnitCellGridType,ZeoliteBeadType,setUpAtom,setUpCellStruct,foldToCenterCell
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: readPDB,writePDB,writePDBmovie

CONTAINS
  SUBROUTINE readPDB(filePDB,zeo,lunitcell,ztype,zcell,zunit,lprint)
    character(LEN=*),intent(in)::filePDB
    type(MoleculeType),intent(out)::zeo
    logical,allocatable,intent(out)::lunitcell(:)
    type(ZeoliteBeadType),intent(out)::ztype
    type(CellMaskType),intent(out)::zcell !< the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(out)::zunit
    LOGICAL,INTENT(IN)::lprint

    INTEGER::IOPDB,jerr,i,resID
    CHARACTER(LEN=default_string_length)::line,atomName,resName,molName,elem
    logical::uninitialized
    real::scoord(3),occup,beta

    IOPDB=get_iounit()
    open(unit=IOPDB,access='sequential',action='read',file=filePDB,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open zeolite PDB file',-1)
    end if

    CALL readLine(IOPDB,line,.false.,jerr)
    IF(jerr.ne.0) call err_exit(__FILE__,__LINE__,'wrong PDB file format',-1)

    read(line(2:),*) zeo%nbead,zunit%dup(1),zunit%dup(2),zunit%dup(3),ztype%ntype
    allocate(zeo%bead(zeo%nbead),lunitcell(zeo%nbead),ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype)&
     ,ztype%num(ztype%ntype),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readPDB: allocation failed',-1)

    do i=1,ztype%ntype
       CALL readLine(IOPDB,line,.false.,jerr)
       IF(jerr.ne.0) call err_exit(__FILE__,__LINE__,'wrong PDB file format',-1)
       read(line(2:),*) ztype%name(i),ztype%type(i),ztype%radiisq(i)
       ztype%radiisq(i)=ztype%radiisq(i)*ztype%radiisq(i)
       ztype%num(i)=0
    end do

    uninitialized=.true.
    i=0
    DO
       CALL readLine(IOPDB,line,.true.,jerr)
       IF(jerr.ne.0) EXIT

       SELECT CASE (line(1:6))
       CASE ("CRYST1")
          read(line(7:),*) zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val,zcell%ang(1)%val,zcell%ang(2)%val&
           ,zcell%ang(3)%val
          call setUpCellStruct(zcell,zunit,lprint,filePDB)
          uninitialized=.false.
       CASE ("ATOM","HETATM")
          if (uninitialized) call err_exit(__FILE__,__LINE__,'PDB: CRYST1 needs to be defined before ATOM',-1)
          i = i + 1
          READ(line(13:16),*) atomname
          READ(line(18:20),*,IOSTAT=jerr) resName
          READ(line(23:26),*,IOSTAT=jerr) resID
          READ(line(31:38),*,IOSTAT=jerr) zeo%bead(i)%coord(1)
          READ(line(39:46),*,IOSTAT=jerr) zeo%bead(i)%coord(2)
          READ(line(47:54),*,IOSTAT=jerr) zeo%bead(i)%coord(3)
          READ(line(55:60),*,IOSTAT=jerr) occup
          READ(line(61:66),*,IOSTAT=jerr) beta
          READ(line(73:76),*,IOSTAT=jerr) molName
          READ(line(77:78),*,IOSTAT=jerr) elem

          IF (LEN_TRIM(elem).eq.0) THEN
             ! Element is assigned on the basis of the atomName
             elem = atomname
          END IF
          IF (LEN_TRIM(molName).eq.0) THEN
             ! If molname is missing (as in the PDB generated by VMD) let's
             ! use the resname for the molname
             molName =  resname
          END IF

          CALL foldToCenterCell(zeo%bead(i)%coord,zcell,scoord)
          call setUpAtom(elem,i,zeo,lunitcell,ztype,zcell,zunit)
       CASE ("END","TER")
          EXIT
       CASE ("REMARK")
       CASE DEFAULT
       END SELECT
    END DO

    if (i.ne.zeo%nbead) call err_exit(__FILE__,__LINE__,'PDB: Number of atoms incorrect',-1)

  END SUBROUTINE readPDB

  SUBROUTINE writePDB(filePDB,zeo,ztype,zcell)
    use const_math,only:raddeg
    character(LEN=*),intent(in)::filePDB
    type(MoleculeType),intent(in)::zeo
    type(ZeoliteBeadType),intent(in)::ztype
    type(CellMaskType),intent(in)::zcell !< the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype

    INTEGER::IOPDB,jerr,i

    IOPDB=get_iounit()
    open(unit=IOPDB,access='sequential',action='write',file=filePDB,form='formatted',iostat=jerr,status='unknown')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open file for writing (PDB)',-1)
    end if

    WRITE(IOPDB,'(A6,3(F9.3),3(F7.2),1X,A10)') "CRYST1",zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val&
     ,zcell%ang(1)%val*raddeg,zcell%ang(2)%val*raddeg,zcell%ang(3)%val*raddeg,"P1        "

    DO i=1,zeo%nbead
       WRITE(IOPDB,'(A6,I5,1X,A4,1X,A3,2X,I4,4X,3(F8.3),2(F6.2),10X,A2)') "ATOM  ",i,ADJUSTL(TRIM(ztype%name(zeo%bead(i)%type)))&
        ,"ZEO",1,zeo%bead(i)%coord,1.0,0.0,ADJUSTL(TRIM(ztype%name(zeo%bead(i)%type)))
    END DO

    WRITE(IOPDB,'(A3,/,A3)') "TER","END"

  END SUBROUTINE writePDB

  subroutine writePDBmovie(unitPDB,boxnum)

      use sim_system
      integer,intent(in)::unitPDB,boxnum
      integer::i,j,ntj,imolty,atomcount,chaincount
      character(LEN=4)::resname
      write(unitPDB,'(A6,3(F9.3),3(F7.2),1X,A10)') "CRYST1",boxlx(boxnum),boxly(boxnum),boxlz(boxnum)&
       ,90.0,90.0,90.0,"P1        "
      
      atomcount  = 0
      chaincount = 0
      do i = 1,nchain
          if (nboxi(i).eq.boxnum) then
              chaincount = chaincount + 1
              imolty = moltyp(i)
              do j = 1,nunit(imolty)
                  atomcount = atomcount + 1
                  ntj = ntype(imolty,j)
                  resname=to_upper(molecname(imolty)(1:4))
                  write(unitPDB,'(A6,I5,1X,A4,1X,A3,2X,I4,4X,3(F8.3),2(F6.2),10X,A2)')"ATOM  ",atomcount&
                   ,ADJUSTL(TRIM(chemid(ntj))),ADJUSTL(TRIM(resname)),chaincount,rxu(i,j),ryu(i,j),rzu(i,j),1.0,0.0&
                   ,ADJUSTL(TRIM(chemid(ntj)))
              end do
          end if
      end do
      
      write(unitPDB,'(A3)') "END"

  end subroutine writePDBmovie

  function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

  end function to_upper       
        
END MODULE parser_pdb
