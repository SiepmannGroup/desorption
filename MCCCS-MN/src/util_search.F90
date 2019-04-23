module util_search
  use util_runtime,only:err_exit
  use util_memory,only:reallocate
  use util_string,only:integer_to_string
  implicit none
  private
  public::LookupTable,initiateTable,destroyTable,addToTable,tightenTable,indexOf,hunt,locate

  type LookupTable
     integer::size !< number of element in the table
     integer,allocatable::list(:) !< correspondence table
  end type LookupTable

CONTAINS
  !> \brief initiate table with specified initial size
  !>
  !> The table size will increase automatically when adding a new item will cause an out-of-bound error
  subroutine initiateTable(table,initialSize)
    type(LookupTable),intent(inout)::table
    integer,intent(in)::initialSize

    integer::jerr

    table%size=0
    allocate(table%list(1:initialSize),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'initiateTable: allocation error',jerr)
  end subroutine initiateTable

  !> \brief Destroy specified table
  subroutine destroyTable(table)
    type(LookupTable),intent(inout)::table

    integer::jerr

    table%size=0
    deallocate(table%list,stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'destroyTable: deallocation error',jerr)
  end subroutine destroyTable

  !> \brief Add item to table
  !>
  !> Ignore if exist, increase table size by 2 if table is full and \a expand is set to .true., signal an error otherwise
  integer function addToTable(table,item,expand)
    type(LookupTable),intent(inout)::table
    integer,intent(in)::item
    logical,intent(in),optional::expand

    integer::i
    logical::lexpand

    if (present(expand)) then
       lexpand=expand
    else
       lexpand=.false.
    end if

    do i=1,table%size
       if (table%list(i).eq.item) exit
    end do

    if (i.gt.table%size) then
       ! item not found
       if (i.gt.size(table%list)) then
          if (lexpand) then
             call reallocate(table%list,1,2*size(table%list))
          else
             call err_exit(__FILE__,__LINE__,'addToTable: table full',-1)
          end if
       end if
       table%size=i
       table%list(i)=item
    end if
    addToTable=i
  end function addToTable

  subroutine tightenTable(table)
    type(LookupTable),intent(inout)::table

    if (table%size.gt.0) then
       call reallocate(table%list,1,table%size)
    end if
  end subroutine tightenTable

  !> \brief look up item in table
  !>
  !> \return index of item in table; 0 if not found
  integer function indexOf(table,item)
    type(LookupTable),intent(in)::table
    integer,intent(in)::item

    integer::i

    do i=1,table%size
       if (table%list(i).eq.item) exit
    end do
    if (i.gt.table%size) i=0
    indexOf=i
  end function indexOf

!> \brief Hunt the index of a value in a monotonic array using hunting algorithm
!>
!> Given an array \a xx(1:n) and a value \a x, return index
!> \a j such that \a x is (insofar as possible) centered in
!> the subrange \a xx(j:j+mm-1). The values in \a xx must be
!> monotonic, either increasing or decreasing. The returned
!> value is not less than 1, nor greater than \a n.
  FUNCTION HUNT(XX,N,X,JO,MM) RESULT(J)
    INTEGER::J
    INTEGER,INTENT(IN)::N,JO,MM
    REAL,INTENT(IN)::XX(N),X
    LOGICAL::ASCND
    INTEGER::JLO,JHI,JM,INC

    IF (N.LT.2.OR.MM.LT.2.OR.MM.GT.N) CALL ERR_EXIT(__FILE__,__LINE__,'HUNT SIZE ERROR',-1)

    ASCND=XX(N).GT.XX(1)

    IF(JLO.LT.1.OR.JLO.GT.N)THEN
       JLO=1
       JHI=N
    ELSE
       INC=1
       JLO=JO
       IF (X.GT.XX(JLO).EQV.ASCND) THEN
          DO
             JHI=JLO+INC
             IF (JHI.GE.N) THEN
                JHI=N
             ELSE IF (X.GT.XX(JHI).EQV.ASCND) THEN
                JLO=JHI
                INC=INC+INC
             ELSE
                EXIT
             END IF
          END DO
       ELSE
          JHI=JLO
          DO
             JLO=JHI-INC
             IF (JLO.LE.1) THEN
                JLO=1
             ELSE IF (X.LE.XX(JLO).EQV.ASCND) THEN
                JHI=JLO
                INC=INC+INC
             ELSE
                EXIT
             END IF
          END DO
       END IF
    END IF
    DO WHILE (JHI-JLO.GT.1)
       JM=(JHI+JLO)/2
       IF (X.GT.XX(JM).EQV.ASCND) THEN
          JLO=JM
       ELSE
          JHI=JM
       END IF
    END DO
    J=MAX(1,MIN(N-MM+1,JLO-(MM-2)/2))
    RETURN
  END FUNCTION HUNT

!> \brief Locate the index of a value in a monotonic array using bisection
!>
!> Given an array \a xx(1:n) and a value \a x, return index
!> \a j such that \a x is (insofar as possible) centered in
!> the subrange \a xx[j:j+mm-1]. The values in \a xx must be
!> monotonic, either increasing or decreasing. The returned
!> value is not less than 1, nor greater than \a n.
  FUNCTION LOCATE(XX,N,X,MM) RESULT(J)
    INTEGER::J
    INTEGER,INTENT(IN)::N,MM
    REAL,INTENT(IN)::XX(N),X
    INTEGER::JL,JM,JU
    LOGICAL::ASCND

    IF (N.LT.2.OR.MM.LT.2.OR.MM.GT.N) CALL ERR_EXIT(__FILE__,__LINE__,'LOCATE SIZE ERROR',-1)

    ASCND=(XX(N).GT.XX(1))

    JL=1
    JU=N
    DO WHILE (JU-JL.GT.1)
       JM=(JU+JL)/2
       IF(ASCND.EQV.(X.GT.XX(JM))) THEN
          JL=JM
       ELSE
          JU=JM
       END IF
    END DO
    J=MAX(1,MIN(N-MM+1,JL-(MM-2)/2))
    RETURN
  END FUNCTION LOCATE
end module util_search
