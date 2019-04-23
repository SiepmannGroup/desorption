    DATA_TYPE, INTENT(IN) :: mydata
    DATA_TYPE, INTENT(OUT) :: alldata(:)
    INTEGER, INTENT(IN) :: comm
#ifndef ALLGATHER
    INTEGER, INTENT(IN) :: root
#endif
    INTEGER :: ierr

#ifdef __MPI__

#ifdef ALLGATHER
    CALL MPI_ALLGATHER(mydata, 1, MP_TYPE, alldata, 1, MP_TYPE, comm, IERR)
#else
    CALL MPI_GATHER(mydata, 1, MP_TYPE, alldata, 1, MP_TYPE, root, comm, IERR)
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    alldata(1) = mydata
#endif
    RETURN
