    DATA_TYPE, INTENT(IN) :: mydata(:)
    DATA_TYPE, INTENT(OUT) :: alldata(:)
    INTEGER, INTENT(IN) :: comm
#ifndef ALLGATHER
    INTEGER, INTENT(IN)::root
#endif
    INTEGER :: msglen, ierr

    msglen = SIZE(mydata)
    !IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop(__LINE__)

#ifdef __MPI__

#ifdef ALLGATHER
    CALL MPI_ALLGATHER(mydata, msglen, MP_TYPE, alldata, msglen, MP_TYPE, comm, IERR)
#else
    CALL MPI_GATHER(mydata, msglen, MP_TYPE, alldata, msglen, MP_TYPE, root, comm, IERR)
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    alldata(1:msglen) = mydata
#endif
    RETURN
