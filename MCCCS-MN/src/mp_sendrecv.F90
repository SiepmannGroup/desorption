    INTEGER, INTENT(IN) :: dest, src, tag
    INTEGER, INTENT(INOUT) :: msglen
    INTEGER, INTENT(IN) :: comm
    INTEGER :: myid
#ifdef __MPI__
    INTEGER :: istatus(MPI_STATUS_SIZE)
    INTEGER :: ierr, nrcv

    CALL MPI_comm_rank( comm, myid, ierr )
    IF( ierr /= 0 ) CALL mp_stop(__LINE__)
#else
    myid = 0
#endif

    IF (src .NE. dest) THEN
#ifdef __MPI__
       IF(myid .EQ. src) THEN
          CALL MPI_SEND( msg_src, msglen, MP_TYPE, dest, tag, comm, ierr)
          IF (ierr/=0) CALL mp_stop(__LINE__)
       ELSE IF(myid .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, msglen, MP_TYPE, src, tag, comm, istatus, ierr )
          IF (ierr/=0) CALL mp_stop(__LINE__)
          CALL MPI_GET_COUNT(istatus, MP_TYPE, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(__LINE__)
          msglen=nrcv
       ELSE
          ! processors not taking part in the communication have 0 length message
          msglen = 0
       END IF
#endif
    ELSE IF (myid .EQ. src)THEN
       msg_dest = msg_src
    END IF

#ifdef __USE_BARRIER
    CALL mp_barrier(comm)
#endif

    RETURN
