!> \file Calculate properties using the Widom test particle insertion method.
!>
!> \warning Widom insertion is known to fail at high densities.

module prop_widom
  use var_type,only:dp,RealAllocArray3D,IntegerAllocArray3D
  implicit none
  private
  save
  public::read_prop_widom,inc_prop_widom_counter,calc_prop_widom,blk_avg_prop_widom,write_deltaG_map,write_prop_widom,write_prop_widom_with_stats,read_checkpoint_prop_widom,write_checkpoint_prop_widom

  type(RealAllocArray3D),allocatable::deltaG(:,:)
  type(RealAllocArray3D),allocatable::deltaG_count(:,:)
  real,allocatable::blk_avg_setedist(:,:,:),setedist(:,:),blk_avg_setedist_prev(:,:),blk_avg_Uads(:,:),Uads(:),blk_avg_Uads_prev(:),Wrosen(:),blk_avg_Wrosen_prev(:)
  logical,allocatable::ldeltaG_map(:)
  integer::deltaG_ngrid(3),idx(3)
  logical::initialized=.false.
  namelist /widom/ ldeltaG_map,deltaG_ngrid
contains
  subroutine read_prop_widom(io_input,lprint,blockm)
    use util_string,only:format_n
    use util_runtime,only:err_exit
    use util_mp,only:mp_bcast
    use sim_system,only:lgrand,nmolty,myid,rootid,groupid,io_output,nunit
    INTEGER,INTENT(IN)::io_input,blockm
    LOGICAL,INTENT(IN)::lprint
    integer::i,j,jerr

    lgrand=.true.

    if (allocated(Wrosen)) then
       do i=1,nmolty
          do j=0,1
             deallocate(deltaG(j,i)%val,deltaG_count(j,i)%val,stat=jerr)
          end do
       end do
       deallocate(Wrosen,blk_avg_Wrosen_prev,Uads,blk_avg_Uads_prev,blk_avg_Uads,setedist,blk_avg_setedist_prev,blk_avg_setedist,ldeltaG_map,deltaG,deltaG_count,stat=jerr)
    end if
    allocate(Wrosen(nmolty),blk_avg_Wrosen_prev(nmolty),Uads(nmolty),blk_avg_Uads_prev(nmolty),blk_avg_Uads(nmolty,blockm),setedist(4,nmolty),blk_avg_setedist_prev(4,nmolty),blk_avg_setedist(4,nmolty,blockm),ldeltaG_map(nmolty),deltaG(0:1,nmolty),deltaG_count(0:1,nmolty),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'read_prop_widom: memory allocation',jerr)

    Wrosen=0._dp
    blk_avg_Wrosen_prev=0._dp
    Uads=0._dp
    blk_avg_Uads_prev=0._dp
    blk_avg_Uads=0._dp
    setedist=0._dp
    blk_avg_setedist_prev=0._dp
    blk_avg_setedist=0._dp
    ldeltaG_map=.false.
    deltaG_ngrid=0

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=widom,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: widom',jerr)
    end if

    call mp_bcast(ldeltaG_map,nmolty,rootid,groupid)
    call mp_bcast(deltaG_ngrid,3,rootid,groupid)

    if (.not.initialized) call init_prop_widom()

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST WIDOM','------------------------------------------'
       write(io_output,'(A,'//format_n(nmolty,'L2')//')') 'ldeltaG_map: ',ldeltaG_map(1:nmolty)
       write(io_output,'(A,3(1X,I0))') 'deltaG_ngrid: ',deltaG_ngrid
    end if
  end subroutine read_prop_widom

  subroutine init_prop_widom()
    use util_runtime,only:err_exit
    use sim_system,only:nmolty
    integer,parameter::ibox=1
    integer::jerr,j,i
    do i=1,nmolty
       if (ldeltaG_map(i)) then
          do j=0,1
             allocate(deltaG(j,i)%val(0:deltaG_ngrid(1)-1,0:deltaG_ngrid(2)-1,0:deltaG_ngrid(3)-1),deltaG_count(j,i)%val(0:deltaG_ngrid(1)-1,0:deltaG_ngrid(2)-1,0:deltaG_ngrid(3)-1),stat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_prop_widom: memory allocation',jerr)
             deltaG(j,i)%val=0._dp
             deltaG_count(j,i)%val=0._dp
          end do
       end if
    end do

    initialized = .true.
  end subroutine init_prop_widom

  subroutine inc_prop_widom_counter(imolty,tag_bead,r,weight)
    use sim_cell,only:hmati
    use sim_system,only:lsolid,lrect,boxlx,boxly,boxlz
    integer,intent(in)::imolty,tag_bead
    real,intent(in)::r(3),weight
    integer,parameter::ibox=1
    real::scoord(3)

    if (ldeltaG_map(imolty)) then
       ! convert to fractional coordinates, fold back to central simulation box, and find the grid indices
       if (lsolid(ibox) .and. .not.lrect(ibox)) then
          scoord(1) = r(1)*hmati(ibox,1)+r(2)*hmati(ibox,4)+r(3)*hmati(ibox,7)
          scoord(2) = r(1)*hmati(ibox,2)+r(2)*hmati(ibox,5)+r(3)*hmati(ibox,8)
          scoord(3) = r(1)*hmati(ibox,3)+r(2)*hmati(ibox,6)+r(3)*hmati(ibox,9)
       else
          scoord(1) = r(1)/boxlx(ibox)
          scoord(2) = r(2)/boxly(ibox)
          scoord(3) = r(3)/boxlz(ibox)
       end if
       scoord = scoord - floor(scoord)
       idx = scoord * deltaG_ngrid
       deltaG_count(tag_bead,imolty)%val(idx(1),idx(2),idx(3)) = deltaG_count(tag_bead,imolty)%val(idx(1),idx(2),idx(3)) + weight
    end if
  end subroutine inc_prop_widom_counter

  subroutine calc_prop_widom(imolty,weight)
    use moves_cbmc,only:first_bead_to_swap
    use sim_system,only:vnew,ivTot,rxnew,rynew,rznew,xcm,ycm,zcm,nugrow,nchain,beta
    integer,intent(in)::imolty
    real,intent(in)::weight
    real::r(3),s,arg
    integer::i,j

    Wrosen(imolty) = Wrosen(imolty) + weight
    Uads(imolty) = Uads(imolty) + vnew(ivTot)*weight
    s=0.0_dp
    do i=1,3
       if (i.eq.1) then
          r(1)=rxnew(1)
          r(2)=rxnew(nugrow(imolty))
       else if (i.eq.2) then
          r(1)=rynew(1)
          r(2)=rynew(nugrow(imolty))
       else if (i.eq.3) then
          r(1)=rznew(1)
          r(2)=rznew(nugrow(imolty))
       end if
       arg = (r(1)-r(2))**2
       s = s + arg
       setedist(i,imolty) = setedist(i,imolty) + arg*weight
    end do
    setedist(0,imolty) = setedist(0,imolty) + s*weight

    if (ldeltaG_map(imolty)) then
       ! Must start with j=1 because we need idx from earlier call to inc_prop_widom_counter()
       do j=1,0,-1
          if (j.gt.0) then
             r(1)=rxnew(first_bead_to_swap(imolty))
             r(2)=rynew(first_bead_to_swap(imolty))
             r(3)=rznew(first_bead_to_swap(imolty))
          else
             i=nchain+2
             r(1)=xcm(i)
             r(2)=ycm(i)
             r(3)=zcm(i)
             call inc_prop_widom_counter(imolty,j,r,weight/exp(-beta*vnew(ivTot)))
          end if
          deltaG(j,imolty)%val(idx(1),idx(2),idx(3)) = deltaG(j,imolty)%val(idx(1),idx(2),idx(3)) + weight
       end do
    end if
  end subroutine calc_prop_widom

  subroutine blk_avg_prop_widom(nblock)
    use util_math,only:store_block_average
    use sim_system,only:nmolty
    integer,intent(in)::nblock
    real::Wrosen_tmp
    integer::itype,j

    do itype = 1, nmolty
       do j = 0, 3
          Wrosen_tmp = blk_avg_Wrosen_prev(itype)
          call store_block_average(blk_avg_setedist(j,itype,nblock),setedist(j,itype)/Wrosen(itype),Wrosen(itype),blk_avg_setedist_prev(j,itype),Wrosen_tmp)
       end do
       call store_block_average(blk_avg_Uads(itype,nblock),Uads(itype)/Wrosen(itype),Wrosen(itype),blk_avg_Uads_prev(itype),blk_avg_Wrosen_prev(itype))
    end do
  end subroutine blk_avg_prop_widom

  subroutine write_deltaG_map(io_output,imolty,prefactor)
    use var_type,only:default_path_length
    use util_files,only:get_iounit
    use util_runtime,only:err_exit
    use moves_cbmc,only:first_bead_to_swap
    integer,intent(in)::io_output,imolty
    real,intent(in)::prefactor
    character(len=default_path_length)::fname
    integer::iunit,i,j,k,io_deltaG,jerr

    if (ldeltaG_map(imolty)) then
       io_deltaG=get_iounit()

       do iunit=0,1
          if (iunit.eq.1) then
             i=first_bead_to_swap(imolty)
          else
             i=iunit
          end if
          write(fname,'("deltaG-",I2.2,"-",I2.2,".dat")') imolty,i
          open(unit=io_deltaG,access='sequential',action='write',file=fname,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) then
             call err_exit(__FILE__,__LINE__,'cannot open file for writing (deltaG)',-1)
          end if

          where (deltaG_count(iunit,imolty)%val.gt.0) deltaG(iunit,imolty)%val = deltaG(iunit,imolty)%val / real(deltaG_count(iunit,imolty)%val,dp) / prefactor

          do k=0,deltaG_ngrid(3)-1
             do j=0,deltaG_ngrid(2)-1
                do i=0,deltaG_ngrid(1)-1
                   write(io_deltaG,'(G16.9,1X,G16.9)') deltaG(iunit,imolty)%val(i,j,k),deltaG_count(iunit,imolty)%val(i,j,k)
                end do
             end do
          end do

          close(io_deltaG)
       end do
    end if
  end subroutine write_deltaG_map

  subroutine write_prop_widom(io_output)
    use sim_system,only:nmolty
    integer,intent(in)::io_output
    integer::i,j

    Uads(1:nmolty)=Uads(1:nmolty)/Wrosen(1:nmolty)
    do i=0,3
       setedist(i,1:nmolty)=setedist(i,1:nmolty)/Wrosen(1:nmolty)
    end do

    write(io_output,"(/,A)") '  ** Properties from Rosenbluth Sampling **'
    do i=1,nmolty
       write(io_output,"(A,I2,A,1X,F12.3)") ' potential energy of type ',i,'          [K] =',Uads(i)
    end do
    do i=1,nmolty
       write(io_output,"(A,I2,A,4(1X,F12.3))") ' sete & x,y,z len of type ',i,'        [A^2] =',(setedist(j,i),j=0,3)
    end do
  end subroutine write_prop_widom

  subroutine write_prop_widom_with_stats(io_output,nblock)
    use util_math,only:calculate_statistics
    use sim_system,only:nmolty
    integer,intent(in)::io_output,nblock
    real::avg(0:3),sd(0:3),se(0:3)
    integer::itype,j

    write(io_output,"(/,A)") '  ** Properties with Uncertainties from Rosenbluth Sampling **'
    do itype = 1, nmolty
       call calculate_statistics(blk_avg_Uads(itype,1:nblock),avg(0),sd(0),se(0))
       write(io_output,"(2(A,I2),A,3(1X,F12.3))") ' potential energy    itype ',itype,' box ',1,' = ',avg(0),sd(0),se(0)
    end do
    do itype = 1, nmolty
       do j = 0, 3
          call calculate_statistics(blk_avg_setedist(j,itype,1:nblock),avg(j),sd(j),se(j))
       end do
       write(io_output,"(2(A,I2),A,12(1X,F7.3))") ' sete & x,y,z length itype ',itype,' box ',1,' = ',(avg(j),sd(j),se(j),j=0,3)
    end do

  end subroutine write_prop_widom_with_stats

  subroutine read_checkpoint_prop_widom(io_chkpt,blockm)
    use util_mp,only:mp_bcast
    use sim_system,only:myid,rootid,groupid,nmolty
    integer,intent(in)::io_chkpt,blockm
    if (myid.eq.rootid) read(io_chkpt) Wrosen,blk_avg_Wrosen_prev,Uads,blk_avg_Uads_prev,blk_avg_Uads,setedist,blk_avg_setedist_prev,blk_avg_setedist
    call mp_bcast(Wrosen,nmolty,rootid,groupid)
    call mp_bcast(blk_avg_Wrosen_prev,nmolty,rootid,groupid)
    call mp_bcast(Uads,nmolty,rootid,groupid)
    call mp_bcast(blk_avg_Uads_prev,nmolty,rootid,groupid)
    call mp_bcast(blk_avg_Uads,nmolty*blockm,rootid,groupid)
    call mp_bcast(setedist,4*nmolty,rootid,groupid)
    call mp_bcast(blk_avg_setedist_prev,4*nmolty,rootid,groupid)
    call mp_bcast(blk_avg_setedist,4*nmolty*blockm,rootid,groupid)
  end subroutine read_checkpoint_prop_widom

  subroutine write_checkpoint_prop_widom(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) Wrosen,blk_avg_Wrosen_prev,Uads,blk_avg_Uads_prev,blk_avg_Uads,setedist,blk_avg_setedist_prev,blk_avg_setedist
  end subroutine write_checkpoint_prop_widom
end module prop_widom
