!> \file
!> Wrapper for MPI implementations that have problems with large messages
!>
!> \note In some MPI implementation the communication subsystem
!> crashes when message exceeds a given size, so we need
!> to break down MPI communications in smaller pieces

#define __MSGSIZ_MAX 100000
#define __BCAST_MSGSIZ_MAX 100000

!>  Some implementation of MPI (OpenMPI) if it is not well tuned for the given
!>  network hardware (InfiniBand) tend to lose performance or get stuck inside
!>  collective routines if processors are not well synchronized
!>  A barrier fixes the problem
!#define __USE_BARRIER

#ifdef __DOUBLE_PRECISION__
#define MP_REAL MPI_DOUBLE_PRECISION
#else
#define MP_REAL MPI_REAL
#endif

  SUBROUTINE mp_synchronize( comm )
    USE util_runtime,ONLY:err_exit
    IMPLICIT NONE
#ifdef __MPI__
  include 'mpif.h'
#endif
    INTEGER, INTENT(IN) :: comm
#if defined __MPI__ && defined __USE_BARRIER
    INTEGER :: ierr
    CALL MPI_barrier( comm, ierr )
    IF( ierr /= 0 ) call err_exit(__FILE__,__LINE__, 'mp_synchronize error in MPI_barrier ', ierr)
#endif
    RETURN
  END SUBROUTINE mp_synchronize

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE LOGICAL

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MPI_LOGICAL
  SUBROUTINE BCAST_LOGICAL( array, n, root, comm )
#include "mp_bcast_aux.F90"
  END SUBROUTINE BCAST_LOGICAL

#ifdef COMMON_BLOCK_NAME
#undef COMMON_BLOCK_NAME
#endif
#define COMMON_BLOCK_NAME mp_aux_logical

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_LAND
  SUBROUTINE mp_reduce_land( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_land

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_LOR
  SUBROUTINE mp_reduce_lor( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_lor

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE INTEGER

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MPI_INTEGER
  SUBROUTINE BCAST_INTEGER( array, n, root, comm )
#include "mp_bcast_aux.F90"
  END SUBROUTINE BCAST_INTEGER

#ifdef COMMON_BLOCK_NAME
#undef COMMON_BLOCK_NAME
#endif
#define COMMON_BLOCK_NAME mp_aux_integer

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_MIN
  SUBROUTINE mp_reduce_min_integer( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_min_integer

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_MAX
  SUBROUTINE mp_reduce_max_integer( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_max_integer

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_SUM
  SUBROUTINE mp_reduce_sum_integer( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_sum_integer

#ifdef COMMON_BLOCK_NAME
#undef COMMON_BLOCK_NAME
#endif
  SUBROUTINE mp_reduce_sum_to_integer( dim, ps, psout, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_sum_to_integer

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE REAL(DP)

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MP_REAL
  SUBROUTINE BCAST_REAL( array, n, root, comm )
#include "mp_bcast_aux.F90"
  END SUBROUTINE BCAST_REAL

#ifdef COMMON_BLOCK_NAME
#undef COMMON_BLOCK_NAME
#endif
#define COMMON_BLOCK_NAME mp_aux_real

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_MIN
  SUBROUTINE mp_reduce_min_real( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_min_real

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_MAX
  SUBROUTINE mp_reduce_max_real( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_max_real

#ifdef MP_OP
#undef MP_OP
#endif
#define MP_OP MPI_SUM
  SUBROUTINE mp_reduce_sum_real( dim, ps, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_sum_real

#ifdef COMMON_BLOCK_NAME
#undef COMMON_BLOCK_NAME
#endif
  SUBROUTINE mp_reduce_sum_to_real( dim, ps, psout, comm, root )
#include "mp_reduce_aux.F90"
  END SUBROUTINE mp_reduce_sum_to_real
