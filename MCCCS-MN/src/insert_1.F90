  if (.not.allocated(p)) then
     allocate(p(pos:pos+999))
  else if (pos.gt.ubound(p,1)) then
     call reallocate(p,lbound(p,1),pos+999)
  else if (pos.lt.lbound(p,1)) then
     call reallocate(p,pos-999,ubound(p,1))
  end if
  p(pos)=val
