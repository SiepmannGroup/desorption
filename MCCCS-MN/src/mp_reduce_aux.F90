!> \file
!> Reduce a distributed variable ps(dim) over the processors.
!> Save the result on either ps or the output variable psout.
!>
!> This version uses a fixed-length buffer of appropriate (?) dim
    USE var_type,ONLY:DP
    USE util_runtime,ONLY:err_exit
    IMPLICIT NONE
#ifdef __MPI__
  include 'mpif.h'
#endif
    INTEGER,  INTENT(IN)    :: dim     !< size of the array
    DATA_TYPE               :: ps(dim) !< array whose elements have to be reduced
    INTEGER,  INTENT(IN)    :: comm    !< communicator
    INTEGER,  INTENT(IN)    :: root    !< if root <  0 perform a reduction to all procs
                                       !< if root >= 0 perform a reduce only to root proc.
    INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
#ifdef COMMON_BLOCK_NAME
    DATA_TYPE :: buff(maxb) !< the use of the common here could help the transfer of data to the network device
    COMMON / COMMON_BLOCK_NAME/ buff
#else
    DATA_TYPE               :: psout(dim)
#endif

#ifdef __MPI__
    INTEGER            :: ierr, ibuf, nbuf, remaining, istart, nproc, myid

#ifdef __DEBUG_MPI__
    write(*,*) 'mp_reduce_aux IN'
#endif
    CALL MPI_comm_size( comm, nproc, ierr )
    IF( ierr /= 0 ) call err_exit(__FILE__,__LINE__, 'mp_reduce_aux error in MPI_comm_size', ierr)

    CALL MPI_comm_rank( comm, myid, ierr )
    IF( ierr /= 0 ) call err_exit(__FILE__,__LINE__, 'mp_reduce_aux error in MPI_comm_rank', ierr)

#ifndef COMMON_BLOCK_NAME
    IF ( dim > 0 .AND. nproc <= 1 ) THEN
       psout = ps
    END IF
#endif

    IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1  ! go to the end of the subroutine

#ifdef __USE_BARRIER
    ! synchronize processes
    CALL mp_synchronize( comm )
#endif

    nbuf = dim / maxb
    DO ibuf = 1, nbuf
       istart = 1 + ( ibuf - 1 ) * maxb
       IF( root >= 0 ) THEN
#ifdef COMMON_BLOCK_NAME
          CALL MPI_REDUCE( ps(istart), buff, maxb, MP_TYPE, MP_OP, root, comm, ierr )
#else
          CALL MPI_REDUCE( ps(istart), psout(istart), maxb, MP_TYPE, MP_OP, root, comm, ierr )
#endif
       ELSE
#ifdef COMMON_BLOCK_NAME
          CALL MPI_ALLREDUCE( ps(istart), buff, maxb, MP_TYPE, MP_OP, comm, ierr )
#else
          CALL MPI_ALLREDUCE( ps(istart), psout(istart), maxb, MP_TYPE, MP_OP, comm, ierr )
#endif
       END IF
       IF( ierr /= 0 ) call err_exit(__FILE__,__LINE__, 'mp_reduce_aux error in mpi_(all)reduce 1', ierr)

#ifdef COMMON_BLOCK_NAME
       IF( root < 0 ) THEN
          ps(istart:(ibuf*maxb)) = buff(1:maxb)
       ELSE IF( root == myid ) THEN
          ps(istart:(ibuf*maxb)) = buff(1:maxb)
       END IF
#endif
    END DO

    ! possible remaining elements < maxb
    remaining = dim - nbuf * maxb
    IF ( remaining > 0 ) THEN
       istart = 1 + nbuf * maxb
       IF( root >= 0 ) THEN
#ifdef COMMON_BLOCK_NAME
          CALL MPI_REDUCE( ps(istart), buff, remaining, MP_TYPE, MP_OP, root, comm, ierr )
#else
          CALL MPI_REDUCE( ps(istart), psout(istart), remaining, MP_TYPE, MP_OP, root, comm, ierr )
#endif
       ELSE
#ifdef COMMON_BLOCK_NAME
          CALL MPI_ALLREDUCE( ps(istart), buff, remaining, MP_TYPE, MP_OP, comm, ierr )
#else
          CALL MPI_ALLREDUCE( ps(istart), psout(istart), remaining, MP_TYPE, MP_OP, comm, ierr )
#endif
       END IF
       IF( ierr /= 0 ) call err_exit(__FILE__,__LINE__, 'mp_reduce_aux error in mpi_(all)reduce 2', ierr)

#ifdef COMMON_BLOCK_NAME
       IF( root < 0 ) THEN
          ps(istart:dim) = buff(1:remaining)
       ELSE IF( root == myid ) THEN
          ps(istart:dim) = buff(1:remaining)
       END IF
#endif
    END IF

1   CONTINUE

#ifdef __DEBUG_MPI__
    write(*,*) 'mp_reduce_aux OUT'
#endif

#endif
    RETURN
