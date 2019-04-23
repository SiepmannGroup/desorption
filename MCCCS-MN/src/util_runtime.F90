module util_runtime
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  private
  public::err_exit
contains
  subroutine err_exit(file,lineno,msg,code)
    character(LEN=*),intent(in)::file,msg
    integer,intent(in)::lineno,code
    integer::ierr

    write(*,FMT='("ERROR in ",A,": line ",I0)') TRIM(file),lineno
    write(*,FMT='("code ",I0,": ",A)') code,msg

#ifdef __MPI__
    call MPI_ABORT(MPI_COMM_WORLD,code,ierr)
#endif

    stop
  end subroutine err_exit
end module util_runtime
