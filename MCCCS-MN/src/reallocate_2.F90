    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, ub1, ub1_old, ub2, ub2_old

    IF (allocated(p)) THEN
       lb1_old = LBOUND(p,1)
       ub1_old = UBOUND(p,1)
       lb2_old = LBOUND(p,2)
       ub2_old = UBOUND(p,2)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       lb2 = MAX(lb2_new,lb2_old)
       ub2 = MIN(ub2_new,ub2_old)

       ALLOCATE (work(lb1:ub1,lb2:ub2),STAT=istat)
       IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)

       work(lb1:ub1,lb2:ub2) = p(lb1:ub1,lb2:ub2)

       DEALLOCATE(p,STAT=istat)
       IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)
    END IF

    ALLOCATE(p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)

    p(:,:) = p0

    IF (ALLOCATED(p).AND.ALLOCATED(work)) THEN
       p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
       DEALLOCATE (work,STAT=istat)
       IF (istat /= 0) call err_exit(__FILE__,__LINE__,'',istat)
    END IF
