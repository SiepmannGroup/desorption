module util_files
  use util_string,only:is_blank_line
  implicit none
  private
  public::get_iounit,readLine,readNthLine,flush_force

  INTEGER,PARAMETER::max_unit_number=999,reserved_unit_numbers(3)=(/0,5,6/)
CONTAINS

! *****************************************************************************
!> \brief returns the first fortran unit that is not preconnected
!> \return -1 if no free unit exists
!> \author taken from CP2K
! *****************************************************************************
  INTEGER FUNCTION get_iounit()
    INTEGER                                  :: istat
    LOGICAL                                  :: lexist, lopen

    DO get_iounit=1,max_unit_number
       IF (ANY(get_iounit == reserved_unit_numbers)) CYCLE
       INQUIRE(UNIT=get_iounit,EXIST=lexist,OPENED=lopen,IOSTAT=istat)
       IF (lexist.AND.(.NOT.lopen).AND.(istat == 0)) RETURN
    END DO

    get_iounit = -1

  END FUNCTION get_iounit

!> \brief read the following line, skipping comment or blank lines
  subroutine readLine(iounit,line,skipComment,iostat)
    integer,intent(in)::iounit
    logical,intent(in)::skipComment
    character(len=*),intent(out)::line
    integer,intent(out)::iostat

    do
       line=''
       read(unit=iounit,FMT='(A)',iostat=iostat) line
       if  (iostat.eq.0) then
          if (.not.is_blank_line(line,skipComment)) exit
       else
          exit
       end if
    end do

  end subroutine readLine

!> \brief read the n-th line (i.e., skipping n-1 lines and read the following line), skipping
!> comment or blank lines
  subroutine readNthLine(iounit,n,line,skipComment,iostat)
    integer,intent(in)::iounit,n
    logical,intent(in)::skipComment
    character(len=*),intent(out)::line
    integer,intent(out)::iostat
    integer::i

    do i=1,n
       call readLine(iounit,line,skipComment,iostat)
       if (iostat.ne.0) exit
    end do

  end subroutine readNthLine

!> \brief Force flushing IO
!>
!> \note Fortran flush seems able to only transfer data from program-specific cache
!> to system cache, which isn't guaranteed to be written to file system if
!> program crashes. The current work-around isn't guaranteed to work either,
!> but in most cases it works. Only the low-level Unix system call (f)sync can
!> indeed force committing data to file system.
  function flush_force(io_unit) result(pos)
    use var_type,only:default_string_length,default_path_length
    use util_runtime,only:err_exit
    integer::pos
    integer,intent(in)::io_unit
    character(LEN=default_string_length)::mode,action,form
    character(LEN=default_path_length)::fname
    integer::jerr
    logical::lopen

    inquire(unit=io_unit,access=mode,action=action,name=fname,form=form,iostat=jerr,opened=lopen,pos=pos)
    if (jerr.ne.0) then
       write(*,'("Error ",I0," in flush_force file: ",A)') jerr,fname
    else if (lopen.and.ALL(io_unit.ne.reserved_unit_numbers)) then
       close(io_unit)
       open(unit=io_unit,access=mode,action=action,file=fname,form=form,iostat=jerr,position="append",status="old")
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Error when opening file '//fname,jerr)
       inquire(unit=io_unit,pos=pos)
    end if
  end function flush_force
end module util_files
