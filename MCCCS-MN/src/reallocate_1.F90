    INTEGER :: istat, lb1, lb1_old, ub1, ub1_old

    IF (allocated(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)

       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)

       work(lb1:ub1) = p(lb1:ub1)

       DEALLOCATE(p,STAT=istat)
       IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)
    END IF

    ALLOCATE(p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)

    p(:) = p0

    IF (ALLOCATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1) = work(lb1:ub1)
       DEALLOCATE (work,STAT=istat)
       IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)
    END IF
