    INTEGER, INTENT(IN) :: msglen
    INTEGER, INTENT(IN) :: comm
#ifdef __MPI__
    CALL MP_AUX_FUNC( msglen, msg, comm, -1 )
#endif
