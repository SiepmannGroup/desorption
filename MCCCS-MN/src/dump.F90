!> \brief Dumps the final configuration before stopping the program.
subroutine dump(file_cfg)
  use sim_system
  use sim_cell
  implicit none
  character(LEN=*),intent(in)::file_cfg

  integer::io_cfg,jerr,i,j,imolty,ibox

  io_cfg=get_iounit()
  open(unit=io_cfg,access='sequential',action='write',file=file_cfg,form='formatted',iostat=jerr,status='unknown')
  if (jerr.ne.0) then
     call err_exit(__FILE__,__LINE__,'cannot open file '//trim(file_cfg),myid+1)
  end if

  write(io_cfg,*) tmcc
  write(io_cfg,*) Armtrax, Armtray, Armtraz
  do ibox=1,nbox
     do imolty=1,nmolty
        write(io_cfg,*) rmtrax(imolty,ibox), rmtray(imolty,ibox) , rmtraz(imolty,ibox)
        write(io_cfg,*) rmrotx(imolty,ibox), rmroty(imolty,ibox) , rmrotz(imolty,ibox)
     end do
  end do
  do ibox=1, nbox
     write (io_cfg,*) (rmflcq(i,ibox),i=1,nmolty)
  end do
  write(io_cfg,*) (rmvol(ibox),ibox=1,nbox)
  do ibox = 1,nbox
     if (lsolid(ibox) .and. .not. lrect(ibox)) then
        write(io_cfg,*) (rmhmat(ibox,i),i=1,9)
        write(io_cfg,*) (hmat(ibox,i),i=1,9)
     else
        write(io_cfg,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
     end if
  end do
  write(io_cfg,*) nchain
  write(io_cfg,*) nmolty
  write(io_cfg,*) (nunit(i),i=1,nmolty)
  write(io_cfg,*) (moltyp(i),i=1,nchain)
  write(io_cfg,*) (nboxi(i),i=1,nchain)
  do i = 1, nmolty
     if ( lexpand(i) ) write(io_cfg,*) eetype(i)
  end do
  do i = 1, nmolty
     if ( lexpand(i) ) write(io_cfg,*) rmexpc(i)
  end do
  do  i = 1, nchain
     imolty = moltyp(i)
     do j = 1, nunit(imolty)
        write(io_cfg,*) rxu(i,j), ryu(i,j), rzu(i,j), qqu(i,j)
     end do
  end do
  close(io_cfg)
  return
end subroutine dump
