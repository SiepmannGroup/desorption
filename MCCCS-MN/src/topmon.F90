program topmon
  use var_type,only:default_path_length
  use util_mp,only:mp_start,mp_end
  use sim_system,only:myid,numprocs,groupid,thread_num,thread_num_max
  use topmon_main,only:output_version,monola
#ifdef __OPENMP__
  use omp_lib
#endif
  implicit none
  integer::narg,iarg,jerr
  character(LEN=default_path_length)::sarg,file_in='topmon.inp'
  logical::lrun=.true.,lsetinput=.false.,lversion=.false.,lusage =.false.
! ----------------------------------------------------------------
! Parse the command line arguments
  narg=command_argument_count()
  iarg=1
  do while (iarg.le.narg)
     call get_command_argument(iarg,sarg)
     select case(sarg)
     case('--version','-v')
        lversion=.true.
        lrun=.false.
     case('--help','-h')
        lusage=.true.
        lrun=.false.
     case('--thread','-t')
        iarg=iarg+1
        if (iarg.gt.narg) then
           lusage=.true.
           lrun=.false.
           exit
        end if
        call get_command_argument(iarg,sarg)
        read(sarg,*,iostat=jerr) thread_num_max
        if (jerr.ne.0) then
           lusage=.true.
           lrun=.false.
           exit
        end if
!$       thread_num=thread_num_max
!$       thread_num_max=omp_get_max_threads()
!$       if (thread_num.gt.thread_num_max) thread_num=thread_num_max
!$       call omp_set_num_threads(thread_num)
     case('--input','-i')
        iarg=iarg+1
        if (iarg.gt.narg) then
           lusage=.true.
           lrun=.false.
           exit
        end if
        call get_command_argument(iarg,sarg)
        lsetinput=.true.
        file_in=sarg
     case default
        if (.not.lsetinput.and.iarg.eq.narg) then
           lsetinput=.true.
           file_in=sarg
        end if
     end select
     iarg=iarg+1
  end do
! ----------------------------------------------------------------
  if (lversion) call output_version(6)

  if (lusage) then
     call get_command_argument(0,sarg)
     write(*,'(A,/,T4,A,/,T4,A,/)') 'Usage: '//trim(sarg)// ' [--version|-v] [--help|-h] [(--threads|-t) number_of_threads_per_processor]','[(--input|-i) name_of_master_input]','If the input file is the last argument, --input or -i can be omitted'
  end if

  if (lrun) then
     call mp_start(numprocs,myid,groupid)

! --- call main program
     call monola(file_in)

     call mp_end()
  end if
! ----------------------------------------------------------------
end program topmon
