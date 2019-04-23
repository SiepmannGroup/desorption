MODULE energy_3body
  use var_type,only:default_string_length
  use const_math,only:degrad
  use util_runtime,only:err_exit
  use util_string,only:uppercase
  use util_files,only:readLine
  use util_memory,only:reallocate
  use util_search
  use sim_system
  use sim_cell
  implicit none
  private
  save
  public::hasThreeBody,U3System,U3MolSys,readThreeBody

  integer,parameter::nparam(1)=(/4/)
  character(LEN=default_path_length)::file_triplets='triplets.lst'
  logical::hasThreeBody=.false.,hasTable=.true.,lbuild_triplet_table=.true.
  namelist /threebody/ file_triplets,lbuild_triplet_table
  real,allocatable::coeffs(:,:) !< coeffs(i,j): j = triplet type, i = 1(LJ cutoff**2)/2(Coulomb cutoff**2)/3(1st coefficient)/4/etc.; j = atom2*nAtomType**2+atom1*nAtomType+atom3 (atom1<=atom3)
!< Or, atom3=mod(itype,nAtomType), atom1=mod((itype-atom3)/nAtomType,nAtomType), atom2=(itype-atom3-atom1*nAtomType)/nAtomType**2
  integer,allocatable::nTriplets(:),ttype(:)& !<functional type
   ,triplets(:,:,:) !< triplet(i,j,k): k = box #, j = j-th triplet, i = 0(type)/1(1st molecule)/2/3/4(unit in 1st molecule)/5/6
  type(LookupTable)::threebodies

CONTAINS
  integer function idx(beadtype)
    integer,intent(in)::beadtype(3)
    integer::atomIdx(3),nAtomType,i

    do i=1,3
       atomIdx(i)=indexOf(atoms,beadtype(i))
    end do

    if (atomIdx(1).gt.atomIdx(3)) then
       i=atomIdx(3)
       atomIdx(3)=atomIdx(1)
       atomIdx(1)=i
    end if

    nAtomType=atoms%size
    !write(*,*) atomIdx,nAtomType
    idx=atomIdx(2)*nAtomType**2+atomIdx(1)*nAtomType+atomIdx(3)
  end function idx

  subroutine readThreeBody(file_ff)
    character(LEN=*),INTENT(IN)::file_ff
    integer,parameter::initial_size=10
    character(LEN=default_string_length)::line
    integer::io_ff,jerr,nEntries,i,beadtype(3)

    !write(*,*) 'three_body::readInput'
    io_ff=get_iounit()
    open(unit=io_ff,access='sequential',action='read',file=file_ff,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open three_body input file',jerr)

    read(UNIT=io_ff,NML=threebody,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: threebody',jerr)

    ! Looking for section THREEBODY
    REWIND(io_ff)
    DO
       call readLine(io_ff,line,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) exit

       if (UPPERCASE(line(1:9)).eq."THREEBODY") then
          call checkAtom()
          allocate(triplets(0:6,1:initial_size,1:nbox),coeffs(1:4,1:initial_size),nTriplets(1:nbox),ttype(1:initial_size),STAT=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'initiateThreeBody: allocation error',jerr)
          call initiateTable(threebodies,initial_size)
          nEntries=0
          do
             call readLine(io_ff,line,skipComment=.true.,iostat=jerr)
             if (UPPERCASE(line(1:13)).eq.'END THREEBODY') exit
             nEntries=nEntries+1
             read(line,*) beadtype
             i=addToTable(threebodies,idx(beadtype),expand=.true.)
             if (i.gt.ubound(ttype,1)) then
                call reallocate(coeffs,1,4,1,2*ubound(coeffs,2))
                call reallocate(ttype,1,2*ubound(ttype,1))
             end if
             read(line,*) beadtype,ttype(i),coeffs(1:nparam(ttype(i)),i)
             write(*,*) threebodies%list(i),ttype(i),coeffs(:,i)
             coeffs(1,i)=coeffs(1,i)**2 !cutoffLJ
             coeffs(2,i)=coeffs(2,i)**2 !cutoffCoul
             coeffs(3,i)=coeffs(3,i)*degrad !degree to radians
          end do
          if (nEntries.gt.0) hasThreeBody=.true.
       end if
    END DO
    close(io_ff)

    if (hasThreeBody) call buildTripletTable()
  end subroutine readThreeBody

  subroutine buildTripletTable()
    integer::ibox,lchain,cchain,rchain,lmolty,cmolty,rmolty,lunit,cunit,runit,nAtomType,itype,io_triplets,nboxtmp,jerr
    real::ar(3),br(3)

    io_triplets=get_iounit()
    if (.not.lbuild_triplet_table) then
       open(unit=io_triplets,access='sequential',action='read',file=file_triplets,form='formatted',iostat=jerr,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot read triplets file',jerr)
       read(io_triplets,*) nboxtmp,nTriplets
       if (nbox.ne.nboxtmp) call err_exit(__FILE__,__LINE__,'triplets file not generated with current input',myid+1)
       do ibox=1,nbox
          read(io_triplets,*) triplets(:,1:nTriplets(ibox),ibox)
       end do
       write(io_output,*) 'Get ',nTriplets,' triplets'
    else
       nAtomType=atoms%size
       do ibox=1,nbox
          if (lpbc) call setpbc(ibox)
          nTriplets(ibox)=0
          do cchain=1,nchain
             if (nboxi(cchain).ne.ibox) cycle
             do lchain=1,nchain
                if (nboxi(lchain).ne.ibox) cycle
                do rchain=lchain,nchain
                   if (nboxi(rchain).ne.ibox) cycle
                   lmolty=moltyp(lchain)
                   cmolty=moltyp(cchain)
                   rmolty=moltyp(rchain)
                   do cunit=1,nunit(cmolty)
                      do lunit=1,nunit(lmolty)
                         do runit=1,nunit(rmolty)
                            if ((lchain.eq.cchain.and.lunit.eq.cunit).or.(lchain.eq.rchain.and.lunit.eq.runit)&
                             .or.(cchain.eq.rchain.and.cunit.eq.runit)) cycle
                            itype=indexOf(threebodies,idx((/ntype(lmolty,lunit),ntype(cmolty,cunit),ntype(rmolty,runit)/)))
                            if (itype.eq.0) cycle !no 3-body defined for this triplet
                            ! write(*,'("Box ",I1,"(",I4,",",I4,",",I4,")->[",I2,",",I2,",",I2,"]:",I2,",",I2,",",I2,". itype=",I5,". Triplet ",I2," of size ",I2)')&
                            !  ibox,lchain,cchain,rchain,lunit,cunit,runit,indexOf(atoms,ntype(lmolty,lunit))&
                            !  ,indexOf(atoms,ntype(cmolty,cunit)),indexOf(atoms,ntype(rmolty,runit))&
                            !  ,itype,nTriplets(ibox),size(triplets,2)
                            ar(1)=rxu(lchain,lunit)-rxu(cchain,cunit)
                            ar(2)=ryu(lchain,lunit)-ryu(cchain,cunit)
                            ar(3)=rzu(lchain,lunit)-rzu(cchain,cunit)
                            br(1)=rxu(rchain,runit)-rxu(cchain,cunit)
                            br(2)=ryu(rchain,runit)-ryu(cchain,cunit)
                            br(3)=rzu(rchain,runit)-rzu(cchain,cunit)
                            if (lpbc) then
                               call mimage(ar(1),ar(2),ar(3),ibox)
                               call mimage(br(1),br(2),br(3),ibox)
                            end if
                            !write(*,'("ar=",D16.9,". br=",D16.9)') dot_product(ar,ar),dot_product(br,br)
                            if (dot_product(ar,ar).le.coeffs(1,itype).and.dot_product(br,br).le.coeffs(1,itype)) then
                               nTriplets(ibox)=nTriplets(ibox)+1
                               write(*,'("Box ",I1,": (",I4,",",I4,",",I4,")->[",I2,",",I2,",",I2,"]: tritype=",I2,". Triplet ",&
                                &I8," of size ",I8,". ar=",D16.9,". br=",D16.9)')&
                                ibox,lchain,cchain,rchain,lunit,cunit,runit,itype,nTriplets(ibox)&
                                ,size(triplets,2),dot_product(ar,ar),dot_product(br,br)
                               if (nTriplets(ibox).gt.size(triplets,2)) call reallocate(triplets,0,6,1,2*size(triplets,2),1,nbox)
                               triplets(0,nTriplets(ibox),ibox)=itype
                               triplets(1,nTriplets(ibox),ibox)=lchain
                               triplets(2,nTriplets(ibox),ibox)=cchain
                               triplets(3,nTriplets(ibox),ibox)=rchain
                               triplets(4,nTriplets(ibox),ibox)=lunit
                               triplets(5,nTriplets(ibox),ibox)=cunit
                               triplets(6,nTriplets(ibox),ibox)=runit
                            end if
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do

       open(unit=io_triplets,access='sequential',action='write',file=file_triplets,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open triplets file',jerr)
       write(io_triplets,*) nbox,nTriplets
       do ibox=1,nbox
          write(io_triplets,*) triplets(:,1:nTriplets(ibox),ibox)
       end do
       close(io_triplets)
    end if
  end subroutine buildTripletTable

!> \brief Return the angle between vector a, b in radians
  function bondAngle(a,b,positive)
    real::bondAngle
    real,intent(in)::a(3),b(3)
    logical,intent(in),optional::positive !< whether to return positive angles only, default: .false. [in]

    bondAngle=acos(dot_product(a,b)/sqrt(dot_product(a,a)*dot_product(b,b)))

    if (present(positive)) then
       if (positive.and.bondAngle.lt.0) bondAngle=-bondAngle
    end if

  end function bondAngle

  function U3(ar,br,type)
    integer,intent(in)::type
    real::U3
    real,intent(in)::ar(3),br(3)

    if (ttype(type).eq.1) then
       U3=coeffs(4,type)*(bondAngle(ar,br)-coeffs(3,type))**2
    end if
  end function U3

  function U3System(ibox)
    real::U3System
    integer,intent(in)::ibox
    integer::type,i
    real::ar(3),br(3)

    if (lpbc) call setpbc(ibox)

    U3System=0.
    if (hasTable) then
      do i=1,nTriplets(ibox)
         ar(1)=rxu(triplets(1,i,ibox),triplets(4,i,ibox))-rxu(triplets(2,i,ibox),triplets(5,i,ibox))
         ar(2)=ryu(triplets(1,i,ibox),triplets(4,i,ibox))-ryu(triplets(2,i,ibox),triplets(5,i,ibox))
         ar(3)=rzu(triplets(1,i,ibox),triplets(4,i,ibox))-rzu(triplets(2,i,ibox),triplets(5,i,ibox))
         br(1)=rxu(triplets(3,i,ibox),triplets(6,i,ibox))-rxu(triplets(2,i,ibox),triplets(5,i,ibox))
         br(2)=ryu(triplets(3,i,ibox),triplets(6,i,ibox))-ryu(triplets(2,i,ibox),triplets(5,i,ibox))
         br(3)=rzu(triplets(3,i,ibox),triplets(6,i,ibox))-rzu(triplets(2,i,ibox),triplets(5,i,ibox))
         if (lpbc) then
            call mimage(ar(1),ar(2),ar(3),ibox)
            call mimage(br(1),br(2),br(3),ibox)
         end if
         type=triplets(0,i,ibox)
         U3System=U3System+U3(ar,br,type)
       end do
    end if
  end function U3System

  function U3MolSys(ichain,istart,iuend,flagon)
    real::U3MolSys
    integer,intent(in)::ichain,istart,iuend,flagon
    integer::type,i,ibox,iunit
    real::ar(3),br(3)

    ibox=nboxi(ichain)
    U3MolSys=0.
    if (hasTable) then
      do i=1,nTriplets(ibox)
         if (triplets(1,i,ibox).eq.ichain) then
            iunit=triplets(4,i,ibox)
            if (iunit.gt.iuend.or.iunit.lt.istart) cycle
            ar(1)=rxuion(iunit,flagon)-rxu(triplets(2,i,ibox),triplets(5,i,ibox))
            ar(2)=ryuion(iunit,flagon)-ryu(triplets(2,i,ibox),triplets(5,i,ibox))
            ar(3)=rzuion(iunit,flagon)-rzu(triplets(2,i,ibox),triplets(5,i,ibox))
            br(1)=rxu(triplets(3,i,ibox),triplets(6,i,ibox))-rxu(triplets(2,i,ibox),triplets(5,i,ibox))
            br(2)=ryu(triplets(3,i,ibox),triplets(6,i,ibox))-ryu(triplets(2,i,ibox),triplets(5,i,ibox))
            br(3)=rzu(triplets(3,i,ibox),triplets(6,i,ibox))-rzu(triplets(2,i,ibox),triplets(5,i,ibox))
         else if (triplets(3,i,ibox).eq.ichain) then
            iunit=triplets(6,i,ibox)
            if (iunit.gt.iuend.or.iunit.lt.istart) cycle
            ar(1)=rxuion(iunit,flagon)-rxu(triplets(2,i,ibox),triplets(5,i,ibox))
            ar(2)=ryuion(iunit,flagon)-ryu(triplets(2,i,ibox),triplets(5,i,ibox))
            ar(3)=rzuion(iunit,flagon)-rzu(triplets(2,i,ibox),triplets(5,i,ibox))
            br(1)=rxu(triplets(1,i,ibox),triplets(4,i,ibox))-rxu(triplets(2,i,ibox),triplets(5,i,ibox))
            br(2)=ryu(triplets(1,i,ibox),triplets(4,i,ibox))-ryu(triplets(2,i,ibox),triplets(5,i,ibox))
            br(3)=rzu(triplets(1,i,ibox),triplets(4,i,ibox))-rzu(triplets(2,i,ibox),triplets(5,i,ibox))
         else if (triplets(2,i,ibox).eq.ichain) then
            iunit=triplets(5,i,ibox)
            if (iunit.gt.iuend.or.iunit.lt.istart) cycle
            ar(1)=rxu(triplets(1,i,ibox),triplets(4,i,ibox))-rxuion(iunit,flagon)
            ar(2)=ryu(triplets(1,i,ibox),triplets(4,i,ibox))-ryuion(iunit,flagon)
            ar(3)=rzu(triplets(1,i,ibox),triplets(4,i,ibox))-rzuion(iunit,flagon)
            br(1)=rxu(triplets(3,i,ibox),triplets(6,i,ibox))-rxuion(iunit,flagon)
            br(2)=ryu(triplets(3,i,ibox),triplets(6,i,ibox))-ryuion(iunit,flagon)
            br(3)=rzu(triplets(3,i,ibox),triplets(6,i,ibox))-rzuion(iunit,flagon)
         else
            cycle
         end if
         if (lpbc) then
            !setpbc(ibox) should have been invoked by the calling subroutine
            call mimage(ar(1),ar(2),ar(3),ibox)
            call mimage(br(1),br(2),br(3),ibox)
         end if
         type=triplets(0,i,ibox)
         U3MolSys=U3MolSys+U3(ar,br,type)
       end do
    end if
  end function U3MolSys
end MODULE energy_3body
