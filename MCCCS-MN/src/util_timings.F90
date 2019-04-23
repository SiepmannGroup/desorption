module util_timings
#ifdef __OPENMP__
  use omp_lib
#endif
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  private
  public::time_date_str,time_init,time_now
  real,save::time_start

contains
!> \brief Returns a datum in human readable format using a standard Fortran routine
!> \since 10.2009 created [Joost VandeVondele]
  function time_date_str()
    CHARACTER::time_date_str*23,date*8,time*10
    CALL DATE_AND_TIME(date=date, time=time)
    time_date_str=date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "//time(1:2)//":"//time(3:4)//":"//time(5:10)
  end function time_date_str

!> \note subroutine cpu_time(real time) is available in Fortran 95
!> and later, while double precision function MPI_WTIME() is
!> available using MPI, double precision function
!> OMP_GET_WTIME() is available with OpenMP, whether or not
!> the program is running with multiple processes or threads
  subroutine time_init()
#ifdef __MPI__
    time_start=MPI_WTIME()
#elif defined __OPENMP__
!$   time_start=omp_get_wtime()
#else
    call cpu_time(time_start)
#endif
  end subroutine time_init

  function time_now()
    real::time_now
#ifdef __MPI__
    time_now=MPI_WTIME()
#elif defined __OPENMP__
!$   time_now=omp_get_wtime()
#else
    call cpu_time(time_now)
#endif
    time_now=time_now-time_start
  end function time_now
end module util_timings
