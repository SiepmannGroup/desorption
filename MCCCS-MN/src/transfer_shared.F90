module transfer_shared
  use sim_system
  implicit none
  private
  save
  public::read_transfer,update_bias,opt_bias,lopt_bias,freq_opt_bias,read_checkpoint_transfer_shared&
   ,write_checkpoint_transfer_shared,gcmc_setup,gcmc_cleanup,gcmc_exchange

  real,allocatable::u_bias_diff(:,:) !<u_bias_diff(i,j) is the bias potential difference that should be applied to box i to achieve equal distribution of particles of type j
  integer,allocatable::num_update_bias(:,:)
  logical,allocatable::lopt_bias(:)
  integer::freq_opt_bias=500
  namelist /transfer/ lopt_bias,freq_opt_bias
contains
  subroutine read_transfer(io_input,lprint)
    use util_string,only:format_n
    use util_runtime,only:err_exit
    use util_mp,only:mp_bcast
    INTEGER,INTENT(IN)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::jerr

    if (allocated(u_bias_diff)) deallocate(u_bias_diff,num_update_bias,lopt_bias,stat=jerr)
    allocate(u_bias_diff(nbox,nmolty),num_update_bias(nbox,nmolty),lopt_bias(nmolty),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'read_transfer: memory allocation',jerr)

    u_bias_diff=0.0_dp
    num_update_bias=0
    lopt_bias=.false.

    if (myid.eq.rootid) then
       rewind(io_input)
       read(UNIT=io_input,NML=transfer,iostat=jerr)
       if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: transfer',jerr)
    end if

    call mp_bcast(lopt_bias,nmolty,rootid,groupid)
    call mp_bcast(freq_opt_bias,1,rootid,groupid)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST TRANSFER','------------------------------------------'
       write(io_output,'(A,'//format_n(nmolty,'L2')//')') 'lopt_bias: ',lopt_bias(1:nmolty)
       write(io_output,'(A,I0)') 'freq_opt_bias: ',freq_opt_bias
    end if
  end subroutine read_transfer

  subroutine update_bias(u_diff,boxrem,boxins,imolty)
    use util_math,only:update_average
    real,intent(in)::u_diff
    integer,intent(in)::boxrem,boxins,imolty

    num_update_bias(boxins,imolty)=num_update_bias(boxins,imolty)+1
    num_update_bias(boxrem,imolty)=num_update_bias(boxrem,imolty)+1
    call update_average(u_bias_diff(boxins,imolty),u_diff/2.0_dp,num_update_bias(boxins,imolty))
    call update_average(u_bias_diff(boxrem,imolty),-u_diff/2.0_dp,num_update_bias(boxrem,imolty))
  end subroutine update_bias

  subroutine opt_bias
    integer::imolty,ibox

    do imolty=1,nmolty
       if (.not.lopt_bias(imolty)) cycle
       do ibox=1,nbox
          eta2(ibox,imolty)=eta2(ibox,imolty)+u_bias_diff(ibox,imolty)
       end do
    end do
    u_bias_diff=0.0_dp
    num_update_bias=0
  end subroutine opt_bias

  subroutine read_checkpoint_transfer_shared(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) num_update_bias,u_bias_diff
    call mp_bcast(num_update_bias,nbox*nmolty,rootid,groupid)
    call mp_bcast(u_bias_diff,nbox*nmolty,rootid,groupid)
  end subroutine read_checkpoint_transfer_shared

  subroutine write_checkpoint_transfer_shared(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) num_update_bias,u_bias_diff
  end subroutine write_checkpoint_transfer_shared

  subroutine gcmc_setup(imolty,boxrem,irem,iparbox)
    use sim_initia,only:setup_molecule_config
    integer,intent(in)::imolty,boxrem
    integer,intent(out)::irem,iparbox
    nchain=nchain+1
    irem=nchain
    if (lrigid(imolty)) call setup_molecule_config(imolty,irem)
    moltyp(irem)=imolty
    nboxi(irem)=boxrem
    temtyp(imolty)=temtyp(imolty)+1
    nchbox(boxrem)=nchbox(boxrem)+1
    ncmt(boxrem,imolty)=ncmt(boxrem,imolty)+1
    iparbox=ncmt(boxrem,imolty)
  end subroutine gcmc_setup

  subroutine gcmc_cleanup(imolty,boxrem)
    integer,intent(in)::imolty,boxrem
    nchain=nchain-1
    temtyp(imolty)=temtyp(imolty)-1
    nchbox(boxrem)=nchbox(boxrem)-1
    ncmt(boxrem,imolty)=ncmt(boxrem,imolty)-1
  end subroutine gcmc_cleanup

  subroutine gcmc_exchange(irem,iremparbox)
    integer,intent(in)::irem,iremparbox
    integer::imolty,ibox,i,jmolty,jbox

    imolty=moltyp(irem)
    ibox=nboxi(irem)

    ! remove reference to irem in parbox & parall
    parbox(ncmt(ibox,imolty),ibox,imolty)=0
    do i=1,temtyp(imolty)
       if (parall(imolty,i).eq.irem) then
          parall(imolty,i)=parall(imolty,temtyp(imolty))
          parall(imolty,temtyp(imolty))=0
          exit
       end if
    end do

    ! copy nchain to irem
    jmolty=moltyp(nchain)
    moltyp(irem)=jmolty
    jbox=nboxi(nchain)
    nboxi(irem)=jbox
    do i=1,nunit(jmolty)
       rcmu(irem)=rcmu(nchain)
       xcm(irem)=xcm(nchain)
       ycm(irem)=ycm(nchain)
       zcm(irem)=zcm(nchain)
       sxcm(irem)=sxcm(nchain)
       sycm(irem)=sycm(nchain)
       szcm(irem)=szcm(nchain)
       rxu(irem,i)=rxu(nchain,i)
       ryu(irem,i)=ryu(nchain,i)
       rzu(irem,i)=rzu(nchain,i)
       qqu(irem,i)=qqu(nchain,i)
    end do

    ! change reference of nchain to irem in parall & parbox
    do i=1,temtyp(jmolty)
       if (parall(jmolty,i).eq.nchain) then
          parall(jmolty,i)=irem
          exit
       end if
    end do
    if (iremparbox.eq.0) then
       do i=1,ncmt(jbox,jmolty)
          if (parbox(i,jbox,jmolty).eq.nchain) then
             parbox(i,jbox,jmolty)=irem
             exit
          end if
       end do
    else
       parbox(iremparbox,jbox,jmolty)=irem
    end if

    nchain = nchain - 1
    temtyp(imolty)=temtyp(imolty)-1
    nchbox(ibox)=nchbox(ibox)-1
    ncmt(ibox,imolty)=ncmt(ibox,imolty)-1
  end subroutine gcmc_exchange
end module transfer_shared
