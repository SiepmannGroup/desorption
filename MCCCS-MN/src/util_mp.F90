module util_mp
  USE var_type,ONLY:DP, long_int
  IMPLICIT NONE
#ifdef __MPI__
  include 'mpif.h'
#endif
  PRIVATE
  PUBLIC::mp_start, mp_end, mp_rank, mp_size, mp_barrier, mp_sendrecv&
   , mp_bcast, mp_gather, mp_allgather, mp_alltoall&
   , mp_min, mp_max, mp_sum, mp_root_sum, mp_land, mp_lor&
   , mp_circular_shift_left&
   , mp_comm_group, mp_comm_split, mp_comm_create, mp_comm_free&
   , mp_group_create, mp_group_free, mp_set_displs

  INTERFACE mp_sendrecv
     MODULE PROCEDURE mp_sendrecv_i0, mp_sendrecv_i1, mp_sendrecv_r0, mp_sendrecv_r1, mp_sendrecv_r2
  END INTERFACE

  INTERFACE mp_bcast
     MODULE PROCEDURE mp_bcast_l0, mp_bcast_l1, mp_bcast_l2, mp_bcast_i0, mp_bcast_i0_long,&
         mp_bcast_i1, mp_bcast_i2, mp_bcast_i3&
      , mp_bcast_r0, mp_bcast_r1, mp_bcast_r2, mp_bcast_r3, mp_bcast_r4, mp_bcast_r5, mp_bcast_c0, mp_bcast_c1
  END INTERFACE

  INTERFACE mp_gather
     MODULE PROCEDURE mp_gather_i0, mp_gather_i1, mp_gatherv_i1, mp_gatherv_i2, mp_gather_r1, mp_gatherv_r1, mp_gatherv_r2
  END INTERFACE

  INTERFACE mp_allgather
     MODULE PROCEDURE mp_allgather_l1, mp_allgatherv_l1, mp_allgather_i0, mp_allgather_i1, mp_allgatherv_i1 &
      ,  mp_allgather_r0, mp_allgather_r1, mp_allgatherv_r1, mp_allgatherv_r2, mp_allgatherv_r3
  END INTERFACE

  INTERFACE mp_alltoall
     MODULE PROCEDURE mp_alltoall_i3
  END INTERFACE

  INTERFACE mp_min
     MODULE PROCEDURE mp_min_i0, mp_min_i1, mp_min_r0, mp_min_r1
  END INTERFACE

  INTERFACE mp_max
     MODULE PROCEDURE mp_max_i0, mp_max_i1, mp_max_r0, mp_max_r1
  END INTERFACE

  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_i0, mp_sum_i1, mp_sum_i2, mp_sum_i3, mp_sum_r0, mp_sum_r1, mp_sum_r2, mp_sum_r3&
      , mp_sum_r4, mp_sum_r5, mp_sum_r22
  END INTERFACE

  INTERFACE mp_root_sum
     MODULE PROCEDURE mp_root_sum_r2
  END INTERFACE

  INTERFACE mp_land
     MODULE PROCEDURE mp_land_l0, mp_land_l1
  END INTERFACE

  INTERFACE mp_lor
     MODULE PROCEDURE mp_lor_l0, mp_lor_l1
  END INTERFACE

  INTERFACE mp_circular_shift_left
     MODULE PROCEDURE mp_circular_shift_left_d2d_int,mp_circular_shift_left_d2d_double
  END INTERFACE

#ifdef __MPI__
  INTEGER,PUBLIC::mp_comm_null=MPI_COMM_NULL,mp_comm_self=MPI_COMM_SELF
#else
  INTEGER,PUBLIC::mp_comm_null=0,mp_comm_self=0
#endif

#ifdef __DOUBLE_PRECISION__
#define MP_REAL MPI_DOUBLE_PRECISION
#else
#define MP_REAL MPI_REAL
#endif

CONTAINS
  SUBROUTINE mp_stop(code)
    INTEGER, INTENT (IN) :: code
    INTEGER :: ierr
#ifdef __MPI__
    CALL MPI_abort(MPI_COMM_WORLD,code,ierr)
#endif
    STOP
  END SUBROUTINE mp_stop

  SUBROUTINE mp_start(numtask, taskid, groupid)
    INTEGER, INTENT (OUT) :: numtask, taskid, groupid
    INTEGER :: ierr

    ierr = 0
    numtask = 1
    taskid = 0
    groupid = 0

#ifdef __MPI__
    CALL MPI_init(ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
    CALL MPI_comm_rank(MPI_COMM_WORLD,taskid,ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
    CALL MPI_comm_size(MPI_COMM_WORLD,numtask,ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
    groupid = MPI_COMM_WORLD
#endif

    RETURN
  END SUBROUTINE mp_start

  SUBROUTINE mp_end()
    INTEGER :: ierr

    ierr = 0

#ifdef __MPI__
    CALL MPI_finalize(ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
    RETURN
  END SUBROUTINE mp_end

  FUNCTION mp_rank( comm )
    INTEGER :: mp_rank
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr

    ierr = 0
    mp_rank = 0
#ifdef __MPI__
    CALL MPI_comm_rank(comm,mp_rank,ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END FUNCTION mp_rank

  FUNCTION mp_size( comm )
    INTEGER :: mp_size
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr

    ierr = 0
    mp_size = 1
#ifdef __MPI__
    CALL MPI_comm_size(comm,mp_size,ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END FUNCTION mp_size

  SUBROUTINE mp_barrier(comm)
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr
#ifdef __MPI__
    CALL MPI_BARRIER(comm,IERR)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END SUBROUTINE mp_barrier

! ===BEGIN MP_SENDRECV===
#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MPI_INTEGER
  SUBROUTINE mp_sendrecv_i0(msg_src,msg_dest,msglen,src,dest,tag,comm)
    INTEGER :: msg_dest, msg_src
#include "mp_sendrecv.F90"
  END SUBROUTINE mp_sendrecv_i0

  SUBROUTINE mp_sendrecv_i1(msg_src,msg_dest,msglen,src,dest,tag,comm)
    INTEGER :: msg_dest(:), msg_src(:)
#include "mp_sendrecv.F90"
  END SUBROUTINE mp_sendrecv_i1

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MP_REAL
  SUBROUTINE mp_sendrecv_r0(msg_src,msg_dest,msglen,src,dest,tag,comm)
    REAL (DP) :: msg_dest, msg_src
#include "mp_sendrecv.F90"
  END SUBROUTINE mp_sendrecv_r0

  SUBROUTINE mp_sendrecv_r1(msg_src,msg_dest,msglen,src,dest,tag,comm)
    REAL (DP) :: msg_dest(:), msg_src(:)
#include "mp_sendrecv.F90"
  END SUBROUTINE mp_sendrecv_r1

  SUBROUTINE mp_sendrecv_r2(msg_src,msg_dest,msglen,src,dest,tag,comm)
    REAL (DP) :: msg_dest(:,:), msg_src(:,:)
#include "mp_sendrecv.F90"
  END SUBROUTINE mp_sendrecv_r2
! ===END MP_SENDRECV===

! ===BEGIN MP_BCAST===
#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC bcast_logical
  SUBROUTINE mp_bcast_l0(msg,msglen,source,comm)
    LOGICAL :: msg
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_l0

  SUBROUTINE mp_bcast_l1(msg,msglen,source,comm)
    LOGICAL :: msg(:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_l1

  SUBROUTINE mp_bcast_l2(msg,msglen,source,comm)
    LOGICAL :: msg(:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_l2

#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC bcast_integer
  SUBROUTINE mp_bcast_i0(msg,msglen,source,comm)
    INTEGER :: msg
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_i0

  SUBROUTINE mp_bcast_i0_long(msg,msglen,source,comm)
    INTEGER(long_int) :: msg
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_i0_long

  SUBROUTINE mp_bcast_i1(msg,msglen,source,comm)
    INTEGER :: msg(:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_i1

  SUBROUTINE mp_bcast_i2(msg,msglen,source,comm)
    INTEGER :: msg(:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_i2

  SUBROUTINE mp_bcast_i3(msg,msglen,source,comm)
    INTEGER :: msg(:,:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_i3

#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC bcast_real
  SUBROUTINE mp_bcast_r0(msg,msglen,source,comm)
    REAL (DP) :: msg
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_r0

  SUBROUTINE mp_bcast_r1(msg,msglen,source,comm)
    REAL (DP) :: msg(:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_r1

  SUBROUTINE mp_bcast_r2(msg,msglen,source,comm)
    REAL (DP) :: msg(:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_r2

  SUBROUTINE mp_bcast_r3(msg,msglen,source,comm)
    REAL (DP) :: msg(:,:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_r3

  SUBROUTINE mp_bcast_r4(msg,msglen,source,comm)
    REAL (DP) :: msg(:,:,:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_r4

  SUBROUTINE mp_bcast_r5(msg,msglen,source,comm)
    REAL (DP) :: msg(:,:,:,:,:)
#include "mp_bcast.F90"
  END SUBROUTINE mp_bcast_r5

  SUBROUTINE mp_bcast_c0(msg,source,comm)
    CHARACTER (len=*) :: msg
    INTEGER :: source
    INTEGER, INTENT(IN) :: comm
    INTEGER :: msglen, ierr, i
    INTEGER, ALLOCATABLE :: imsg(:)
#ifdef __MPI__
    ierr = 0
    msglen = len(msg)
    ALLOCATE (imsg(1:msglen), STAT=ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
    DO i = 1, msglen
       imsg(i) = ichar(msg(i:i))
    END DO
    CALL bcast_integer( imsg, msglen, source, comm )
    DO i = 1, msglen
       msg(i:i) = char(imsg(i))
    END DO
    DEALLOCATE (imsg, STAT=ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END SUBROUTINE mp_bcast_c0

  SUBROUTINE mp_bcast_c1(msg,source,comm)
    CHARACTER (len=*) :: msg(:)
    INTEGER :: source
    INTEGER, INTENT(IN) :: comm
    INTEGER :: msglen, m1, m2, ierr, i, j
    INTEGER, ALLOCATABLE :: imsg(:,:)
#ifdef __MPI__
    ierr = 0
    m1 = LEN(msg)
    m2 = SIZE(msg)
    msglen = LEN(msg)*SIZE(msg)
    ALLOCATE (imsg(1:m1,1:m2), STAT=ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
    DO j = 1, m2
       DO i = 1, m1
          imsg(i,j) = ichar(msg(j)(i:i))
       END DO
    END DO
    CALL bcast_integer( imsg, msglen, source, comm )
    DO j = 1, m2
       DO i = 1, m1
          msg(j)(i:i) = char(imsg(i,j))
       END DO
    END DO
    DEALLOCATE (imsg, STAT=ierr)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END SUBROUTINE mp_bcast_c1
! ===END MP_BCAST===

! ===BEGIN MP_GATHER===
#ifdef ALLGATHER
#undef ALLGATHER
#endif

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE INTEGER

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MPI_INTEGER
  SUBROUTINE mp_gather_i0(mydata, alldata, root, comm)
#include "mp_gather_0.F90"
  END SUBROUTINE mp_gather_i0

  SUBROUTINE mp_gather_i1(mydata, alldata, root, comm)
#include "mp_gather_1.F90"
  END SUBROUTINE mp_gather_i1

  SUBROUTINE mp_gatherv_i1( mydata, alldata, recvcount, displs, root, comm)
#include "mp_gatherv_1.F90"
  END SUBROUTINE mp_gatherv_i1

  SUBROUTINE mp_gatherv_i2( mydata, alldata, recvcount, displs, root, comm)
#include "mp_gatherv_2.F90"
  END SUBROUTINE mp_gatherv_i2

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE REAL(DP)

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MP_REAL
  SUBROUTINE mp_gather_r1(mydata, alldata, root, comm)
#include "mp_gather_1.F90"
  END SUBROUTINE mp_gather_r1

  SUBROUTINE mp_gatherv_r1( mydata, alldata, recvcount, displs, root, comm)
#include "mp_gatherv_1.F90"
  END SUBROUTINE mp_gatherv_r1

  SUBROUTINE mp_gatherv_r2( mydata, alldata, recvcount, displs, root, comm)
#include "mp_gatherv_2.F90"
  END SUBROUTINE mp_gatherv_r2
! ===END MP_GATHER===

! ===BEGIN MP_ALLGATHER===
#ifndef ALLGATHER
#define ALLGATHER
#endif

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE LOGICAL

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MPI_LOGICAL
  SUBROUTINE mp_allgather_l1(mydata, alldata, comm)
#include "mp_gather_1.F90"
  END SUBROUTINE mp_allgather_l1

  SUBROUTINE mp_allgatherv_l1(mydata, alldata, recvcount, displs, comm)
#include "mp_gatherv_1.F90"
  END SUBROUTINE mp_allgatherv_l1

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE INTEGER

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MPI_INTEGER
  SUBROUTINE mp_allgather_i0(mydata, alldata, comm)
#include "mp_gather_0.F90"
  END SUBROUTINE mp_allgather_i0

  SUBROUTINE mp_allgather_i1(mydata, alldata, comm)
#include "mp_gather_1.F90"
  END SUBROUTINE mp_allgather_i1

  SUBROUTINE mp_allgatherv_i1(mydata, alldata, recvcount, displs, comm)
#include "mp_gatherv_1.F90"
  END SUBROUTINE mp_allgatherv_i1

#ifdef DATA_TYPE
#undef DATA_TYPE
#endif
#define DATA_TYPE REAL(DP)

#ifdef MP_TYPE
#undef MP_TYPE
#endif
#define MP_TYPE MP_REAL
  SUBROUTINE mp_allgather_r0(mydata, alldata, comm)
#include "mp_gather_0.F90"
  END SUBROUTINE mp_allgather_r0

  SUBROUTINE mp_allgather_r1(mydata, alldata, comm)
#include "mp_gather_1.F90"
  END SUBROUTINE mp_allgather_r1

  SUBROUTINE mp_allgatherv_r1(mydata, alldata, recvcount, displs, comm)
#include "mp_gatherv_1.F90"
  END SUBROUTINE mp_allgatherv_r1

  SUBROUTINE mp_allgatherv_r2(mydata, alldata, recvcount, displs, comm)
#include "mp_gatherv_2.F90"
  END SUBROUTINE mp_allgatherv_r2

  SUBROUTINE mp_allgatherv_r3(mydata, alldata, recvcount, displs, comm)
#include "mp_gatherv_3.F90"
  END SUBROUTINE mp_allgatherv_r3
! ===END MP_ALLGATHER===

! ===BEGIN MP_ALLTOALL===
  SUBROUTINE mp_alltoall_i3( sndbuf, rcvbuf, comm )
    INTEGER :: sndbuf( :, :, : )
    INTEGER :: rcvbuf( :, :, : )
    INTEGER, INTENT(IN) :: comm
    INTEGER :: nsiz, ierr, npe

#ifdef __MPI__
    CALL MPI_comm_size( comm, npe, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    IF ( SIZE( sndbuf, 3 ) < npe ) CALL mp_stop(__LINE__)
    IF ( SIZE( rcvbuf, 3 ) < npe ) CALL mp_stop(__LINE__)

    nsiz = SIZE( sndbuf, 1 ) * SIZE( sndbuf, 2 )

    CALL MPI_ALLTOALL( sndbuf, nsiz, MPI_INTEGER, &
     rcvbuf, nsiz, MPI_INTEGER, comm, ierr )

    IF (ierr/=0) CALL mp_stop(__LINE__)
#else
    rcvbuf = sndbuf
#endif

    RETURN
  END SUBROUTINE mp_alltoall_i3
! ===END MP_ALLTOALL===

! ===BEGIN MP_MIN===
#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_min_integer
  SUBROUTINE mp_min_i0(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_min_i0

  SUBROUTINE mp_min_i1(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_min_i1

#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_min_real
  SUBROUTINE mp_min_r0(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_min_r0

  SUBROUTINE mp_min_r1(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_min_r1
! ===END MP_MIN===

! ===BEGIN MP_MAX===
#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_max_integer
  SUBROUTINE mp_max_i0(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_max_i0

  SUBROUTINE mp_max_i1(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_max_i1

#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_max_real
  SUBROUTINE mp_max_r0(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_max_r0

  SUBROUTINE mp_max_r1(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_max_r1
! ===END MP_MAX===

! ===BEGIN MP_SUM===
#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_sum_integer
  SUBROUTINE mp_sum_i0(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_i0

  SUBROUTINE mp_sum_i1(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_i1

  SUBROUTINE mp_sum_i2(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg(:,:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_i2

  SUBROUTINE mp_sum_i3(msg,msglen,comm)
    INTEGER, INTENT (INOUT) :: msg(:,:,:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_i3

#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_sum_real
  SUBROUTINE mp_sum_r0(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_r0

  SUBROUTINE mp_sum_r1(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_r1

  SUBROUTINE mp_sum_r2(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:,:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_r2

  SUBROUTINE mp_sum_r3(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:,:,:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_r3

  SUBROUTINE mp_sum_r4(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:,:,:,:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_r4

  SUBROUTINE mp_sum_r5(msg,msglen,comm)
    REAL (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_sum_r5

  SUBROUTINE mp_sum_r22( msg, res, root, comm )
    REAL (DP), INTENT (IN) :: msg(:,:)
    REAL (DP), INTENT (OUT) :: res(:,:)
    INTEGER, OPTIONAL, INTENT (IN) :: root
    INTEGER, INTENT (IN) :: comm
    INTEGER :: msglen
    INTEGER :: myid, ierr

    msglen = size(msg)

#ifdef __MPI__
    IF( PRESENT( root ) ) THEN
       !
       CALL MPI_comm_rank( comm, myid, ierr)
       IF( ierr /= 0 ) CALL mp_stop(__LINE__)

       IF( myid == root ) THEN
          IF( msglen > size(res) ) CALL mp_stop(__LINE__)
       END IF
       !
       CALL mp_reduce_sum_to_real( msglen, msg, res, comm, root )
       !
    ELSE
       !
       IF( msglen > size(res) ) CALL mp_stop(__LINE__)
       !
       CALL mp_reduce_sum_to_real( msglen, msg, res, comm, -1 )
       !
    END IF
#else
    res = msg
#endif
  END SUBROUTINE mp_sum_r22
! ===END MP_SUM===

! ===BEGIN MP_ROOT_SUM===
  SUBROUTINE mp_root_sum_r2( msg, res, root, comm )
    REAL (DP), INTENT (IN)  :: msg(:,:)
    REAL (DP), INTENT (OUT) :: res(:,:)
    INTEGER,   INTENT (IN)  :: root
    INTEGER, INTENT (IN) :: comm
    INTEGER :: msglen, ierr, myid

#ifdef __MPI__
    msglen = size(msg)

    CALL MPI_comm_rank( comm, myid, ierr)
    IF( ierr /= 0 ) CALL mp_stop(__LINE__)
    !
    IF( myid == root ) THEN
       IF( msglen > size(res) ) CALL mp_stop(__LINE__)
    END IF

    CALL mp_reduce_sum_to_real( msglen, msg, res, comm, root )
#else
    res = msg
#endif
  END SUBROUTINE mp_root_sum_r2
! ===END MP_ROOT_SUM===

! ===BEGIN MP_LAND===
#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_land
  SUBROUTINE mp_land_l0(msg,msglen,comm)
    LOGICAL, INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_land_l0

  SUBROUTINE mp_land_l1(msg,msglen,comm)
    LOGICAL, INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_land_l1
! ===END MP_LAND===

! ===BEGIN MP_LOR===
#ifdef MP_AUX_FUNC
#undef MP_AUX_FUNC
#endif
#define MP_AUX_FUNC mp_reduce_lor
  SUBROUTINE mp_lor_l0(msg,msglen,comm)
    LOGICAL, INTENT (INOUT) :: msg
#include "mp_reduce.F90"
  END SUBROUTINE mp_lor_l0

  SUBROUTINE mp_lor_l1(msg,msglen,comm)
    LOGICAL, INTENT (INOUT) :: msg(:)
#include "mp_reduce.F90"
  END SUBROUTINE mp_lor_l1
! ===END MP_LOR===

! ===BEGIN MP_CIRCULAR_SHIFT_LEFT===
  SUBROUTINE mp_circular_shift_left_d2d_int( buf, itag, comm )
    INTEGER :: buf
    INTEGER, INTENT(IN) :: itag
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr, npe, sour, dest, myid

#ifdef __MPI__
    INTEGER :: istatus( MPI_status_size )

    CALL MPI_comm_size( comm, npe, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
    CALL MPI_comm_rank( comm, myid, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    sour = myid + 1
    IF( sour == npe ) sour = 0
    dest = myid - 1
    IF( dest == -1 ) dest = npe - 1

    CALL MPI_Sendrecv_replace( buf, 1, MPI_INTEGER, &
     dest, itag, sour, itag, comm, istatus, ierr)

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    ! do nothing
#endif
    RETURN
  END SUBROUTINE mp_circular_shift_left_d2d_int

  SUBROUTINE mp_circular_shift_left_d2d_double( buf, itag, comm )
    REAL(DP) :: buf( :, : )
    INTEGER, INTENT(IN) :: itag
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr, npe, sour, dest, myid

#ifdef __MPI__
    INTEGER :: istatus( MPI_status_size )

    CALL MPI_comm_size( comm, npe, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
    CALL MPI_comm_rank( comm, myid, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    sour = myid + 1
    IF( sour == npe ) sour = 0
    dest = myid - 1
    IF( dest == -1 ) dest = npe - 1

    CALL MPI_Sendrecv_replace( buf, SIZE(buf), MP_REAL, &
     dest, itag, sour, itag, comm, istatus, ierr)

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    ! do nothing
#endif
    RETURN
  END SUBROUTINE mp_circular_shift_left_d2d_double
! ===END MP_CIRCULAR_SHIFT_LEFT===

  SUBROUTINE mp_comm_group( comm, group )
    INTEGER, INTENT (IN) :: comm
    INTEGER, INTENT (OUT) :: group
    INTEGER :: ierr
    ierr = 0
#ifdef __MPI__
    CALL MPI_comm_group( comm, group, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
#else
    group = 0
#endif
  END SUBROUTINE mp_comm_group

  SUBROUTINE mp_comm_split( old_comm, color, key, new_comm )
    INTEGER, INTENT (IN) :: old_comm
    INTEGER, INTENT (IN) :: color, key
    INTEGER, INTENT (OUT) :: new_comm
    INTEGER :: ierr
    ierr = 0
#ifdef __MPI__
    CALL MPI_COMM_SPLIT( old_comm, color, key, new_comm, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
#else
    new_comm = old_comm
#endif
  END SUBROUTINE mp_comm_split

  SUBROUTINE mp_comm_create( old_comm, new_grp, new_comm )
    INTEGER, INTENT (IN) :: old_comm
    INTEGER, INTENT (IN) :: new_grp
    INTEGER, INTENT (OUT) :: new_comm
    INTEGER :: ierr

    ierr = 0
    new_comm = old_comm
#ifdef __MPI__
    CALL MPI_comm_create( old_comm, new_grp, new_comm, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END SUBROUTINE mp_comm_create

  SUBROUTINE mp_comm_free( comm )
    INTEGER, INTENT (INOUT) :: comm
    INTEGER :: ierr
    ierr = 0
#ifdef __MPI__
    IF( comm /= MPI_COMM_NULL ) THEN
       CALL MPI_comm_free( comm, ierr )
       IF (ierr/=0) CALL mp_stop(__LINE__)
    END IF
#endif
    RETURN
  END SUBROUTINE mp_comm_free

  SUBROUTINE mp_group_create( group_list, group_size, old_grp, new_grp )
    INTEGER, INTENT (IN) :: group_list(:), group_size, old_grp
    INTEGER, INTENT (OUT) :: new_grp
    INTEGER :: ierr

    ierr = 0
    new_grp = old_grp
#ifdef __MPI__
    CALL MPI_group_incl( old_grp, group_size, group_list, new_grp, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END SUBROUTINE mp_group_create

  SUBROUTINE mp_group_free( comm )
    INTEGER, INTENT (INOUT) :: comm
    INTEGER :: ierr
    ierr = 0
#ifdef __MPI__
    CALL MPI_group_free( comm, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif
  END SUBROUTINE mp_group_free

!> \brief  Given the number of elements on each processor (recvcount), this subroutine
!>  sets the correct offsets (displs) to collect them on a single
!>  array with contiguous elemets
  SUBROUTINE mp_set_displs( recvcount, displs, ntot, nproc )
    INTEGER, INTENT(IN) :: recvcount(:) !< number of elements on each processor
    INTEGER, INTENT(OUT) :: displs(:)   !< offsets/displacements
    INTEGER, INTENT(OUT) :: ntot
    INTEGER, INTENT(IN) :: nproc
    INTEGER :: i

    displs( 1 ) = 0
    !
#ifdef __MPI__
    IF( nproc < 1 ) CALL mp_stop(__LINE__)
    DO i = 2, nproc
       displs( i ) = displs( i - 1 ) + recvcount( i - 1 )
    END DO
    ntot = displs( nproc ) + recvcount( nproc )
#else
    ntot = recvcount( 1 )
#endif
    RETURN
  END SUBROUTINE mp_set_displs
END MODULE util_mp
