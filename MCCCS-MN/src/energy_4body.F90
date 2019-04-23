MODULE energy_4body
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
  public::hasFourBody,U4System,U4MolSys,readFourBody

  integer,parameter::nparam(1)=(/5/)
  character(LEN=default_path_length)::file_quadruplets='quadruplets.lst'
  logical::hasFourBody=.false.,hasTable=.true.,lbuild_quadruplet_table=.true.
  namelist /fourbody/ file_quadruplets,lbuild_quadruplet_table
  real,allocatable::coeffs(:,:) !< coeffs(i,j): j = quadruplet type, i = 1(LJ cutoff**2)/2(Coulomb cutoff**2)/3(1st coefficient)/4/etc.; j = atom2*nAtomType**3+atom3*nAtomType**2+atom1*nAtomType+atom4 (atom2<atom3 or atom1<=atom4 if atom2==atom3)
!< Or, atom4=mod(itype,nAtomType), atom1=mod((itype-atom4)/nAtomType,nAtomType), atom3=(itype-atom4-atom1*nAtomType)/nAtomType**2, atom2=(itype-atom4-atom1*nAtomType-atom3*nAtomType**2)/nAtomType**3
  integer,allocatable::nQuadruplets(:),qtype(:)& !<functional type
   ,quadruplets(:,:,:) !< quadruplet(i,j,k): k = box #, j = j-th quadruplet, i = 0(type)/1(1st molecule)/2/3/4/5(unit in 1st molecule)/6/7/8

  type(LookupTable)::fourbodies

CONTAINS
  integer function idx(beadtype)
    integer,intent(in)::beadtype(4)
    integer::atomIdx(4),nAtomType,i

    do i=1,4
       atomIdx(i)=indexOf(atoms,beadtype(i))
    end do

    if (atomIdx(2).gt.atomIdx(3) .or. (atomIdx(2).eq.atomIdx(3).and.atomIdx(1).gt.atomIdx(4))) then
       i=atomIdx(3)
       atomIdx(3)=atomIdx(2)
       atomIdx(2)=i
       i=atomIdx(4)
       atomIdx(4)=atomIdx(1)
       atomIdx(1)=i
    end if

    nAtomType=atoms%size
    idx=atomIdx(2)*nAtomType**3+atomIdx(3)*nAtomType**2+atomIdx(1)*nAtomType+atomIdx(4)
  end function idx

  subroutine readFourBody(file_ff)
    character(LEN=*),INTENT(IN)::file_ff
    integer,parameter::initial_size=10
    character(LEN=default_string_length)::line
    integer::io_ff,jerr,nEntries,i,beadtype(4)

    !write(*,*) 'four_body::readInput'
    io_ff=get_iounit()
    open(unit=io_ff,access='sequential',action='read',file=file_ff,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open four_body input file',jerr)

    read(UNIT=io_ff,NML=fourbody,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: fourbody',jerr)

    ! Looking for section FOURBODY
    REWIND(io_ff)
    DO
       call readLine(io_ff,line,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) exit

       if (UPPERCASE(line(1:8)).eq."FOURBODY") then
          call checkAtom()
          allocate(quadruplets(0:8,initial_size,1:nbox),coeffs(1:5,1:initial_size),nQuadruplets(1:nbox),qtype(1:initial_size)&
           ,STAT=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'initiateFourBody: allocation error',jerr)
          call initiateTable(fourbodies,initial_size)
          nEntries=0
          do
             call readLine(io_ff,line,skipComment=.true.,iostat=jerr)
             if (UPPERCASE(line(1:12)).eq.'END FOURBODY') exit
             nEntries=nEntries+1
             read(line,*) beadtype
             i=addToTable(fourbodies,idx(beadtype),expand=.true.)
             if (i.gt.ubound(qtype,1)) then
                call reallocate(coeffs,1,5,1,2*ubound(coeffs,2))
                call reallocate(qtype,1,2*ubound(qtype,1))
             end if
             read(line,*) beadtype,qtype(i),coeffs(1:nparam(qtype(i)),i)
             write(*,*) fourbodies%list(i),qtype(i),coeffs(:,i)
             coeffs(1,i)=coeffs(1,i)**2 !cutoffLJ
             coeffs(2,i)=coeffs(2,i)**2 !cutoffCoul
             coeffs(3,i)=coeffs(3,i)*degrad !degree to radians
          end do
          if (nEntries.gt.0) hasFourBody=.true.
       end if
    END DO
    close(io_ff)

    if (hasFourBody) call buildQuadrupletTable()
  end subroutine readFourBody

  subroutine buildQuadrupletTable()
    integer::ibox,lchain,clchain,crchain,rchain,lmolty,clmolty,crmolty,rmolty,lunit,clunit,crunit,runit,nAtomType,itype&
     ,io_quadruplets,nboxtmp,jerr
    real::ar(3),br(3)

    io_quadruplets=get_iounit()
    if (.not.lbuild_quadruplet_table) then
       open(unit=io_quadruplets,access='sequential',action='read',file=file_quadruplets,form='formatted',iostat=jerr,status='old')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot read quadruplets file',jerr)
       read(io_quadruplets,*) nboxtmp,nQuadruplets
       if (nbox.ne.nboxtmp) call err_exit(__FILE__,__LINE__,'quadruplets file not generated with current input',myid+1)
       write(*,*) nbox,nQuadruplets
       do ibox=1,nbox
          read(io_quadruplets,*) quadruplets(:,1:nQuadruplets(ibox),ibox)
          write(*,*) quadruplets(:,1:nQuadruplets(ibox),ibox)
       end do
       write(io_output,*) 'Get ',nQuadruplets,' quadruplets'
    else
       nAtomType=atoms%size
       do ibox=1,nbox
          if (lpbc) call setpbc(ibox)
          nQuadruplets(ibox)=0
          do clchain=1,nchain
             if (nboxi(clchain).ne.ibox) cycle
             do crchain=clchain,nchain
                if (nboxi(crchain).ne.ibox) cycle
                do lchain=1,nchain
                   if (nboxi(lchain).ne.ibox) cycle
                   do rchain=lchain,nchain
                      if (nboxi(rchain).ne.ibox) cycle
                      lmolty=moltyp(lchain)
                      clmolty=moltyp(clchain)
                      crmolty=moltyp(crchain)
                      rmolty=moltyp(rchain)
                      do clunit=1,nunit(clmolty)
                         do crunit=1,nunit(crmolty)
                            do lunit=1,nunit(lmolty)
                               do runit=1,nunit(rmolty)
                                  if ((lchain.eq.clchain.and.lunit.eq.clunit)&
                                   .or.(lchain.eq.crchain.and.lunit.eq.crunit)&
                                   .or.(lchain.eq.rchain.and.lunit.eq.runit)&
                                   .or.(clchain.eq.crchain.and.clunit.eq.crunit)&
                                   .or.(clchain.eq.rchain.and.clunit.eq.runit)&
                                   .or.(crchain.eq.rchain.and.crunit.eq.runit)) cycle
                                  itype=indexOf(fourbodies,idx((/ntype(lmolty,lunit),ntype(clmolty,clunit),ntype(crmolty,crunit)&
                                   ,ntype(rmolty,runit)/)))
                                  if (itype.eq.0) cycle !no 4-body defined for this quadruplet
                                  ! write(*,'("Box ",I1,"(",I4,",",I4,",",I4,")->[",I2,",",I2,",",I2,"]:",I2,",",I2,",",I2,". itype=",I5,". Quadruplet ",I2," of size ",I2)')&
                                  !  ibox,lchain,cchain,rchain,lunit,cunit,runit,indexOf(atoms,ntype(lmolty,lunit))&
                                  !  ,indexOf(atoms,ntype(cmolty,cunit)),indexOf(atoms,ntype(rmolty,runit))&
                                  !  ,itype,nQuadruplets(ibox),size(quadruplets,2)
                                  ar(1)=rxu(lchain,lunit)-rxu(clchain,clunit)
                                  ar(2)=ryu(lchain,lunit)-ryu(clchain,clunit)
                                  ar(3)=rzu(lchain,lunit)-rzu(clchain,clunit)
                                  br(1)=rxu(rchain,runit)-rxu(crchain,crunit)
                                  br(2)=ryu(rchain,runit)-ryu(crchain,crunit)
                                  br(3)=rzu(rchain,runit)-rzu(crchain,crunit)
                                  if (lpbc) then
                                     call mimage(ar(1),ar(2),ar(3),ibox)
                                     call mimage(br(1),br(2),br(3),ibox)
                                  end if
                                  !write(*,'("ar=",D16.9,". br=",D16.9)') dot_product(ar,ar),dot_product(br,br)
                                  if (dot_product(ar,ar).le.coeffs(1,itype).and.dot_product(br,br).le.coeffs(1,itype)) then
                                     nQuadruplets(ibox)=nQuadruplets(ibox)+1
                                     write(*,'("Box ",I1,": (",I4,",",I4,",",I4,",",I4,")->[",I2,",",I2,",",I2,",",I2,"]&
                                      &: quatype=",I2,". Quadruplet ",I8," of size ",I8,". ar=",D16.9,". br=",D16.9)')&
                                      ibox,lchain,clchain,crchain,rchain,lunit,clunit,crunit,runit&
                                      ,indexOf(atoms,ntype(lmolty,lunit)),itype,nQuadruplets(ibox)&
                                      ,size(quadruplets,2),dot_product(ar,ar),dot_product(br,br)
                                     if (nQuadruplets(ibox).gt.size(quadruplets,2))&
                                      call reallocate(quadruplets,0,8,1,2*size(quadruplets,2),1,nbox)
                                     quadruplets(0,nQuadruplets(ibox),ibox)=itype
                                     quadruplets(1,nQuadruplets(ibox),ibox)=lchain
                                     quadruplets(2,nQuadruplets(ibox),ibox)=clchain
                                     quadruplets(3,nQuadruplets(ibox),ibox)=crchain
                                     quadruplets(4,nQuadruplets(ibox),ibox)=rchain
                                     quadruplets(5,nQuadruplets(ibox),ibox)=lunit
                                     quadruplets(6,nQuadruplets(ibox),ibox)=clunit
                                     quadruplets(7,nQuadruplets(ibox),ibox)=crunit
                                     quadruplets(8,nQuadruplets(ibox),ibox)=runit
                                  end if
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do

       open(unit=io_quadruplets,access='sequential',action='write',file=file_quadruplets,form='formatted',iostat=jerr&
        ,status='unknown')
       if (jerr.ne.0) then
          call err_exit(__FILE__,__LINE__,'cannot open quadruplets file',jerr)
       end if
       write(io_quadruplets,*) nbox,nQuadruplets
       do ibox=1,nbox
          write(io_quadruplets,*) quadruplets(:,1:nQuadruplets(ibox),ibox)
       end do
       close(io_quadruplets)
    end if
  end subroutine buildQuadrupletTable

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

  function U4(ar,br,type)
    integer,intent(in)::type
    real::U4
    real,intent(in)::ar(3),br(3)

    if (qtype(type).eq.1) then
       U4=coeffs(4,type)*(bondAngle(ar,br)-coeffs(3,type))**2
    end if
  end function U4

  function U4System(ibox)
    real::U4System
    integer,intent(in)::ibox
    integer::type,i
    real::ar(3),br(3)

    if (lpbc) call setpbc(ibox)

    U4System=0.
    if (hasTable) then
      do i=1,nQuadruplets(ibox)
         ar(1)=rxu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-rxu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         ar(2)=ryu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-ryu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         ar(3)=rzu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-rzu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         br(1)=rxu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-rxu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         br(2)=ryu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-ryu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         br(3)=rzu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-rzu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         if (lpbc) then
            call mimage(ar(1),ar(2),ar(3),ibox)
            call mimage(br(1),br(2),br(3),ibox)
         end if
         type=quadruplets(0,i,ibox)
         U4System=U4System+U4(ar,br,type)
       end do
    end if
  end function U4System

  function U4MolSys(ichain,istart,iuend,flagon)
    real::U4MolSys
    integer,intent(in)::ichain,istart,iuend,flagon
    integer::type,i,ibox,iunit
    real::ar(3),br(3)

    ibox=nboxi(ichain)
    U4MolSys=0.
    if (hasTable) then
      do i=1,nQuadruplets(ibox)
         if (quadruplets(1,i,ibox).eq.ichain) then
            iunit=quadruplets(4,i,ibox)
            if (iunit.gt.iuend.or.iunit.lt.istart) cycle
            ar(1)=rxuion(iunit,flagon)-rxu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            ar(2)=ryuion(iunit,flagon)-ryu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            ar(3)=rzuion(iunit,flagon)-rzu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            br(1)=rxu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-rxu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            br(2)=ryu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-ryu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            br(3)=rzu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-rzu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         else if (quadruplets(3,i,ibox).eq.ichain) then
            iunit=quadruplets(6,i,ibox)
            if (iunit.gt.iuend.or.iunit.lt.istart) cycle
            ar(1)=rxuion(iunit,flagon)-rxu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            ar(2)=ryuion(iunit,flagon)-ryu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            ar(3)=rzuion(iunit,flagon)-rzu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            br(1)=rxu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-rxu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            br(2)=ryu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-ryu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
            br(3)=rzu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-rzu(quadruplets(2,i,ibox),quadruplets(5,i,ibox))
         else if (quadruplets(2,i,ibox).eq.ichain) then
            iunit=quadruplets(5,i,ibox)
            if (iunit.gt.iuend.or.iunit.lt.istart) cycle
            ar(1)=rxu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-rxuion(iunit,flagon)
            ar(2)=ryu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-ryuion(iunit,flagon)
            ar(3)=rzu(quadruplets(1,i,ibox),quadruplets(4,i,ibox))-rzuion(iunit,flagon)
            br(1)=rxu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-rxuion(iunit,flagon)
            br(2)=ryu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-ryuion(iunit,flagon)
            br(3)=rzu(quadruplets(3,i,ibox),quadruplets(6,i,ibox))-rzuion(iunit,flagon)
         else
            cycle
         end if
         if (lpbc) then
            !setpbc(ibox) should have been invoked by the calling subroutine
            call mimage(ar(1),ar(2),ar(3),ibox)
            call mimage(br(1),br(2),br(3),ibox)
         end if
         type=quadruplets(0,i,ibox)
         U4MolSys=U4MolSys+U4(ar,br,type)
       end do
    end if
  end function U4MolSys
end MODULE energy_4body
