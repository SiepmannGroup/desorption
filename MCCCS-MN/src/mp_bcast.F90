    INTEGER, INTENT(IN) :: source, msglen
    INTEGER, INTENT(IN) :: comm
#ifdef __MPI__
    CALL MP_AUX_FUNC( msg, msglen, source, comm )
#endif
