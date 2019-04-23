MODULE sim_initia
  use util_runtime,only:err_exit
  use sim_system
  implicit none
  private
  public::get_molecule_config,setup_molecule_config,setup_system_config, setup_cbmc_bend
  ! Q.Paul C.--added setup_cbmc_bend in the above line for tabulated CBMC bending growth

  real,allocatable,dimension(:,:)::samx,samy,samz

CONTAINS
  subroutine get_molecule_config(file_struct,linit)
    use util_mp,only:mp_bcast
    character(LEN=*),intent(in)::file_struct
    logical,intent(in)::linit

    integer::jerr,io_struct,imolty,m,ichain
    logical::lusefile

    if (allocated(samx)) deallocate(samx,samy,samz,stat=jerr)
    allocate(samx(ntmax,numax),samy(ntmax,numax),samz(ntmax,numax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'store_molecule_config: allocating sam{x,y,z} failed',jerr)

    lusefile=.false.
    do imolty=1,nmolty
       if ((lbranch(imolty).or.lrigid(imolty)).and.(linit.or.(lgrand.and.temtyp(imolty).eq.0))) then
          lusefile=.true.
          exit
       end if
    end do

    if (lusefile.and.myid.eq.rootid) then
       ! read sample structure from unit 78 -
       io_struct=get_iounit()
       open(unit=io_struct,access='sequential',action='read',file=file_struct,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open structure file '//trim(file_struct),myid+1)
    end if

    do imolty=1,nmolty
       if (lbranch(imolty).or.lrigid(imolty)) then
          if (lusefile) then
             if (myid.eq.rootid) then
                read(io_struct,*) !first read comment line for imolty
                do m = 1, nunit(imolty)
                   read(io_struct,*) samx(imolty,m),samy(imolty,m),samz(imolty,m)
                end do
             end if
             call mp_bcast(samx(imolty,:),nunit(imolty),rootid,groupid)
             call mp_bcast(samy(imolty,:),nunit(imolty),rootid,groupid)
             call mp_bcast(samz(imolty,:),nunit(imolty),rootid,groupid)
          else if (temtyp(imolty).gt.0) then
             ichain=parall(imolty,1)
             do m = 1, nunit(imolty)
                samx(imolty,m) = rxu(ichain,m)
                samy(imolty,m) = ryu(ichain,m)
                samz(imolty,m) = rzu(ichain,m)
             end do
          end if
          do m = 2, nunit(imolty)
             samx(imolty,m) = samx(imolty,m) - samx(imolty,1)
             samy(imolty,m) = samy(imolty,m) - samy(imolty,1)
             samz(imolty,m) = samz(imolty,m) - samz(imolty,1)
          end do
          samx(imolty,1) = 0._dp
          samy(imolty,1) = 0._dp
          samz(imolty,1) = 0._dp
       end if
    end do

    if (lusefile.and.myid.eq.rootid) close(io_struct)
  end subroutine get_molecule_config

  subroutine setup_molecule_config(imolty,ichain)
    use moves_cbmc,only:rosenbluth,schedule
    integer,intent(in)::imolty,ichain

    real::ddum
    integer::m,ifrom,nbxp1
    logical::lterm

    moltyp(ichain) = imolty
    do m = 1,nunit(imolty)
       qqu(ichain,m) = qelect(ntype(imolty,m))
    end do

    if (lrigid(imolty)) then
       do m = 1,nunit(imolty)
          rxu(ichain,m)=samx(imolty,m)
          ryu(ichain,m)=samy(imolty,m)
          rzu(ichain,m)=samz(imolty,m)
       end do
       return
    end if

    if (nunit(imolty) .ne. nugrow(imolty)) then
       call err_exit(__FILE__,__LINE__,'Cannot grow molecule. Please provide a structure',myid+1)
    end if

    ! put the first bead at the origin
    rxnew(1) = 0.0E0_dp
    rynew(1) = 0.0E0_dp
    rznew(1) = 0.0E0_dp

    ! determine the growth schedule
    call schedule(nugrow(imolty),imolty,ifrom,1,0,2,0)

    ! actually grow the structure
    nbxp1 = nbox   + 1
    lideal(nbxp1) = .true.
    nboxi(ichain) = nbxp1
    call rosenbluth(.true.,lterm,ichain,ichain,imolty,ifrom,nbxp1,nugrow(imolty),ddum,.false.,ddum,2)
    if (lterm) call err_exit(__FILE__,__LINE__&
     ,'Error in gen_molecule_config growing molecule: maybe increasing nchoi would help?',myid+1)

    ! assign the coordinates
    do m = 1,nunit(imolty)
       rxu(ichain,m)=rxnew(m)
       ryu(ichain,m)=rynew(m)
       rzu(ichain,m)=rznew(m)
    end do
  end subroutine setup_molecule_config

subroutine setup_system_config(file_struct)
  use const_math,only:degrad
  use util_random,only:random
  use energy_intramolecular,only:vtorso
  use moves_cbmc,only:explct

  character(LEN=*),intent(in)::file_struct

  integer::temtypma(ntmax),mcmt(ntmax,nbxmax),bmap(numax),imap(numax),ibox,unitc,imolty,m,nn,n,kc,ic,jc,intemp,ibuild,zzz,prev&
   ,m1,m2,ibtype
  real::ux(nbxmax),uy(nbxmax),uz(nbxmax),xtemp(numax),ytemp(numax),ztemp(numax),xshift,dic,rot,angold,angnew,xynext&
   ,xnext,ynext,znext,vdummy
  logical::lgrown(ntmax),lgrow,lacc(numax)

  if (myid.eq.rootid) then
     write(io_output,'(/,A,/,A)') 'Generating Initial Structures','------------------------------------------'
  end if

  ncmt = 0
  nchbox = 0

  ! calculation of unit cell dimensions
  do ibox = 1,nbox
     nchbox(ibox)=sum(ininch(1:nmolty,ibox))
     unitc = inix(ibox)*iniy(ibox)*iniz(ibox)
     if ( nchbox(ibox) .gt. unitc ) then
        call err_exit(__FILE__,__LINE__,'unit cell too small in box '//integer_to_string(ibox),myid+1)
     end if

     ux(ibox) = boxlx(ibox) / inix(ibox)
     uy(ibox) = boxly(ibox) / iniy(ibox)
     uz(ibox) = boxlz(ibox) / iniz(ibox)
     if (myid.eq.rootid) then
        write(io_output,'(A,I0,A,2(I0,"X"),I0)') 'Box ',ibox,': ',inix(ibox),iniy(ibox),iniz(ibox)
        write(io_output,'(A,2(G16.9,"X"),G16.9)') 'Dimension: ',boxlx(ibox),boxly(ibox),boxlz(ibox)
        write(io_output,'(A,2(G16.9,"X"),G16.9)') 'Spacing: ',ux(ibox),uy(ibox),uz(ibox)
     end if
  end do

  mcmt(1:nmolty,1:nbox)=ininch(1:nmolty,1:nbox)
  temtypma(1) = temtyp(1)
  do imolty = 2, nmolty
     temtypma(imolty) = temtypma(imolty-1) + temtyp(imolty)
  end do

  call get_molecule_config(file_struct,linit=.true.)

  ! calculate coordinates
  do imolty = 1, nmolty
     lgrown(imolty) = .true.
     if (.not.lbranch(imolty)) then
        ! if lbranch is false but the molecule is not linear attempt
        ! to grow it with cbmc
        lgrow = .false.
        do m = 1,nunit(imolty)
           if (invib(imolty,m) .gt. 2) then
              lgrow = .true.
              exit
           end if
        end do
        if (lgrow) then
           if (myid.eq.rootid) then
              write(io_output,'(A,I0)') 'growing a sample structure with CBMC for molecule type ',imolty
           end if

           call setup_molecule_config(imolty,nchain+1)

           do m=1,nunit(imolty)
              samx(imolty,m) = rxu(nchain+1,m)
              samy(imolty,m) = ryu(nchain+1,m)
              samz(imolty,m) = rzu(nchain+1,m)
           end do
        else
           lgrown(imolty) = .false.
        end if
     end if
  end do

  do_ibox:do ibox = 1,nbox
     n = 0
     intemp = 1
     do kc = 0, iniz(ibox)-1
        if ( mod(kc,2) .eq. 0) then
           xshift = 0.0E0_dp
        else
           xshift = dshift(ibox)
        end if

        do ic = 0, inix(ibox)-1
           do jc = 0, iniy(ibox)-1
              if ( mod(jc,2) .eq. 0) then
                 dic = 0.0E0_dp
              else
                 dic = 0.5E0_dp
              end if

              n=n+1
              if (n.gt.nchbox(ibox)) cycle do_ibox

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! inimix = 0 : take molecules at random
! inimix > 0 : take molecules in order (first type I etc.)
! inimix < 0 : take molecules in alternating order
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if (inimix(ibox).gt.0) then
                 do while (mcmt(intemp,ibox).le.0)
                    intemp = intemp + 1
                 end do
              else if (inimix(ibox).lt.0) then
                 intemp = intemp + 1
                 if (intemp.gt.nmolty) intemp = 1
                 do while (mcmt(intemp,ibox).le.0)
                    intemp = intemp + 1
                    if (intemp.gt.nmolty) intemp = 1
                 end do
              else
                 intemp=int(random(-1)*nmolty)+1
                 do while (mcmt(intemp,ibox).le.0)
                    intemp=int(random(-1)*nmolty)+1
                 end do
              end if
              mcmt(intemp,ibox)=mcmt(intemp,ibox)-1

              nn=temtypma(intemp)
              temtypma(intemp)=temtypma(intemp)-1

              rxu(nn,1) = (real(ic,dp)+dic)*ux(ibox)+xshift
              ryu(nn,1) = real(jc,dp)*uy(ibox)
              rzu(nn,1) = real(kc,dp)*uz(ibox)+zshift(ibox)
              nboxi(nn) = ibox
              ncmt(ibox,intemp) = ncmt(ibox,intemp) + 1

              if (lgrown(intemp)) then
                 ibuild = nunit(intemp)
              else
                 ibuild = nugrow(intemp)
              end if

              if (.not.lgrown(intemp)) then
! *************************************************
! start determination of linear chain positions
! allowing for numbering out of order
! note: doesn't exactly create equilibrium structure with
! respect to bond angles or torsions, but that will shake out
! with CBMC anyways.  Should at least take away the overlaps of
! the previous method
! *************************************************
! first need to determine re-mapped bead order- search through connectivity
!
! call the results map(i) where i=1 is one chain end, and its
! value is equal to the bead number of that end
!
! for example, methanol oxygen first, then hydrogen, then CH3
!
! H---O--CH3
!
! bead numbers:   2 - 1 - 3
!
! bmap(1) = 2
! bmap(2) = 1
! bmap(3) = 3
!
! the inverse map is just the opposite:
!
! imap(1) = 2
! imap(2) = 1
! imap(3) = 3

                 ! initialize accounted for variable
                 lacc = .false.

                 ! first find the end with the lowest number
                 zzz = ibuild
                 do m = 1,ibuild
                    if (invib(intemp,m) .le. 1) then
                       if (m .le. zzz) then
                          zzz = m
                       end if
                    else if (invib(intemp,m) .gt. 2) then
                       call err_exit(__FILE__,__LINE__&
                        ,'initia only works for linear molecules! Maybe you should make a fort.78 file and use lbranch?',myid+1)
                    end if
                 end do

                 bmap(1) = zzz
                 imap(zzz) = 1
                 lacc(zzz) = .true.

                 ! now determine the rest
                 do m = 2,ibuild
                    prev = bmap(m-1)
                    do zzz = 1,ibuild
                       if (ijvib(intemp,zzz,1).eq.prev .or. ijvib(intemp,zzz,2).eq.prev) then
                          if (.not.lacc(zzz)) then
                             bmap(m) = zzz
                             imap(zzz) = m
                             lacc(zzz) = .true.
                          end if
                       end if
                    end do
                 end do

                 ! now use old method with re-mapped numbers:
                 ! put first end at origin:
                 xtemp(1) = 0.0E0_dp
                 ytemp(1) = 0.0E0_dp
                 ztemp(1) = 0.0E0_dp

                 ! now we need to loop over all the other beads:
                 do m = 2, ibuild
                    m1 = m - 1
                    m2 = m - 2

                    if ( inirot(ibox) .eq. 0 ) then
                       rot = random(-1) * 360.0E0_dp * degrad
                    else if ( inirot(ibox) .gt. 0 ) then
                       rot = real(inirot(ibox),dp) * degrad
                    else
                       if ( mod(jc,2) .eq. 0 ) then
                          rot = real(inirot(ibox),dp) * degrad
                       else
                          rot = -(real(inirot(ibox),dp) * degrad)
                       end if
                    end if

                    if ( inben(intemp,bmap(m1)) .gt. 0 ) then
                       ibtype = itben(intemp,bmap(m1),1)
                       angold = brben(ibtype) / 2.0E0_dp
                       if ( m .gt. 2 ) ibtype = itben(intemp,bmap(m2),1)
                       angnew = brben(ibtype) - angold
                       ! write(io_output,*) 'angold',angold*raddeg,'   angnew',angnew*raddeg
                       angold = angnew

                       ! need to search for proper bond length
                       do zzz = 1,invib(intemp,bmap(m))
                          if (ijvib(intemp,bmap(m),zzz).eq.bmap(m1)) then
                             ibtype = itvib(intemp,bmap(m),zzz)
                          end if
                       end do

                       ztemp(m) = sin(angnew) * brvib(ibtype)
                       xynext = cos(angnew) * brvib(ibtype)
                       ! write(io_output,*) 'znext',znext,'   xynext',xynext
                    else
                       ! need to search for proper bond length
                       do zzz = 1,invib(intemp,bmap(m))
                          if (ijvib(intemp,bmap(m),zzz).eq.bmap(m1)) then
                             ibtype = itvib(intemp,bmap(m),zzz)
                          end if
                       end do
                       ztemp(m) = brvib(ibtype)
                       xynext = 0.0E0_dp
                    end if

                    if ( mod(m,2) .eq. 0 ) then
                       xtemp(m) = cos(rot) * xynext
                       ytemp(m) = sin(rot) * xynext
                    else
                       xtemp(m) = -(cos(rot) * xynext)
                       ytemp(m) = -(sin(rot) * xynext)
                    end if

                    xtemp(m) = xtemp(m1) + xtemp(m)
                    ytemp(m) = ytemp(m1) + ytemp(m)
                    ztemp(m) = ztemp(m1) + ztemp(m)
                 end do
              end if
! *************************************************
! end linear determination
! *************************************************

              do m = 2, ibuild
                 m1 = m - 1
                 ! write(io_output,*) 'intemp',intemp
                 if ( lgrown(intemp) ) then
                    ! branched molecule with sample structure -
                    xnext = samx(intemp,m) -samx(intemp,m1)
                    ynext = samy(intemp,m) -samy(intemp,m1)
                    znext = samz(intemp,m) -samz(intemp,m1)
                 else
                    ! linear molecule determined above- replacing old code that follows.
                    xnext = xtemp(imap(m)) - xtemp(imap(m1))
                    ynext = ytemp(imap(m)) - ytemp(imap(m1))
                    znext = ztemp(imap(m)) - ztemp(imap(m1))
                 end if

                 rxu(nn,m) = rxu(nn,m1) + xnext
                 ryu(nn,m) = ryu(nn,m1) + ynext
                 rzu(nn,m) = rzu(nn,m1) + znext
              end do
           end do
        end do
     end do
  end do do_ibox

  do n=1,nchain
     imolty = moltyp(n)
     if ( nugrow(imolty) .ne. nunit(imolty) ) then
        call explct(n,vdummy,.true.,.false.)
     end if
     ! set up intial charges on the atoms
     do m = 1,nunit(imolty)
        qqu(n,m) = qelect(ntype(imolty,m))
     end do
  end do

  return
end subroutine setup_system_config

! --- Q.Paul C. --- read input file for tabulated CBMC bending growth
subroutine setup_cbmc_bend(file_cbmc_bend)
  character(LEN=*),intent(in)::file_cbmc_bend
  integer::jerr,io_cbmc_bend
  integer::i,j,k,count
  integer::trash
  integer::lin_num,br_num,itwobranch
  integer::ilin,ibr
  integer::max_lin_bend_dim,max_br_bend_dim1,max_br_bend_dim2,max_br_bend_dim3
        
  io_cbmc_bend=get_iounit()
  open(unit=io_cbmc_bend,access='sequential',action='read',file=file_cbmc_bend &
     ,form='formatted',iostat=jerr,status='unknown')
  if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open tabulated alkane bending table' &
     //trim(file_cbmc_bend),myid+1)
    
! ---First read general configuration
     read(io_cbmc_bend,*)
     read(io_cbmc_bend,*) lin_num,br_num,itwobranch
     allocate(lin_bend_type(lin_num))
     allocate(lin_bend_dim(lin_num))
     allocate(br_bend_type(br_num))
     allocate(br_bend_dim1(br_num))
     allocate(br_bend_dim2(br_num))
     allocate(br_bend_dim3(br_num))
     
! ---Then read linear part
     if ( lin_num .ge. 1) then
        read(io_cbmc_bend,*) 
        do ilin=1,lin_num
           read(io_cbmc_bend,*) 
           read(io_cbmc_bend,*) lin_bend_type(ilin),lin_bend_dim(ilin)
           if ( ilin .eq. 1) then
              max_lin_bend_dim=lin_bend_dim(ilin)
           else
              max_lin_bend_dim=max(lin_bend_dim(ilin),max_lin_bend_dim)
           end if
        end do
     
     
        allocate(lin_bend_table(ilin,max_lin_bend_dim))
        allocate(lin_bend_prob(ilin,max_lin_bend_dim))
      
        do ilin=1,lin_num
           read(io_cbmc_bend,*)
           do i=1,lin_bend_dim(ilin)
              read(io_cbmc_bend,*) lin_bend_table(ilin,i),lin_bend_prob(ilin,i)
           end do
        end do
     end if

! ---Secondly read single-branched part
     ! first the dimensionality information
     if ( br_num .ge. 1) then
        read(io_cbmc_bend,*)
        do ibr=1,br_num
           read(io_cbmc_bend,*)
           read(io_cbmc_bend,*) br_bend_type(ibr),br_bend_dim1(ibr),&
              br_bend_dim2(ibr),br_bend_dim3(ibr)
           if ( ibr .eq. 1) then
              max_br_bend_dim1=br_bend_dim1(ibr)
              max_br_bend_dim2=br_bend_dim2(ibr)
              max_br_bend_dim3=br_bend_dim3(ibr)
           else
              max_br_bend_dim1=max(br_bend_dim1(ibr),max_br_bend_dim1)
              max_br_bend_dim2=max(br_bend_dim2(ibr),max_br_bend_dim2)
              max_br_bend_dim3=max(br_bend_dim3(ibr),max_br_bend_dim3)
           end if 
        end do
        
        allocate(br_bend_theta1(br_num,max_br_bend_dim1+1))
        allocate(br_bend_theta2(br_num,max_br_bend_dim1+1,max_br_bend_dim2+1))
        allocate(br_bend_phi12(br_num,max_br_bend_dim1+1,max_br_bend_dim2+1,max_br_bend_dim3+1))
        allocate(br_bend_prob(br_num,max_br_bend_dim1*max_br_bend_dim2*max_br_bend_dim3+1))
        
        ! secondly the probability part
        do ibr=1,br_num
           read(io_cbmc_bend,*)
           i=2
           j=2
           k=2
           br_bend_prob(ibr,1)=0
           read(io_cbmc_bend,*) br_bend_theta1(ibr,1),br_bend_theta1(ibr,2),br_bend_theta2(ibr,1,1),&
              br_bend_theta2(ibr,1,2),br_bend_phi12(ibr,1,1,1),&
              br_bend_phi12(ibr,1,1,2),br_bend_prob(ibr,2)
                
           do count=2,br_bend_dim1(ibr)*br_bend_dim2(ibr)*br_bend_dim3(ibr)
              if (mod(count,br_bend_dim2(ibr)*br_bend_dim3(ibr)) .eq. 1) then
                 i=i+1
                 j=2
                 k=2
                 read(io_cbmc_bend,*) trash,br_bend_theta1(ibr,i),br_bend_theta2(ibr,i-1,j-1),br_bend_theta2(ibr,i-1,j),&
                    br_bend_phi12(ibr,i-1,j-1,k-1),br_bend_phi12(ibr,i-1,j-1,k),br_bend_prob(ibr,count+1)
              else if (mod(count,br_bend_dim3(ibr)) .eq. 1) then
                 j=j+1
                 k=2
                 read(io_cbmc_bend,*) trash,trash,trash,br_bend_theta2(ibr,i-1,j),&
                    br_bend_phi12(ibr,i-1,j-1,k-1),br_bend_phi12(ibr,i-1,j-1,k),br_bend_prob(ibr,count+1)
              else
                 k=k+1
                 read(io_cbmc_bend,*) trash,trash,trash,trash,trash,&
                    br_bend_phi12(ibr,i-1,j-1,k),br_bend_prob(ibr,count+1)
              end if
           end do
        end do
     end if
    
     close(io_cbmc_bend)
        
end subroutine setup_cbmc_bend
! --- Q.Paul C. ---

END MODULE sim_initia
