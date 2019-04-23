module util_memory
  use util_string,only:integer_to_string
  use util_runtime,only:err_exit
  implicit none
  private
  public::reallocate,insert

  interface reallocate
     module procedure reallocate_l1,reallocate_l2,reallocate_l3,reallocate_l4,reallocate_l5&
      ,reallocate_i1,reallocate_i2,reallocate_i3,reallocate_i4,reallocate_i5&
      ,reallocate_r1,reallocate_r2,reallocate_r3,reallocate_r4,reallocate_r5&
      ,reallocate_c1,reallocate_c2
  end interface

  interface insert
     module procedure insert_to_array_l1,insert_to_array_c1
  end interface

contains
 SUBROUTINE reallocate_l1(p,lb1_new,ub1_new)
    LOGICAL,ALLOCATABLE,INTENT(INOUT)::p(:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new

    LOGICAL,ALLOCATABLE::work(:)
    LOGICAL,PARAMETER::p0=.FALSE.
#include "reallocate_1.F90"
  END subroutine reallocate_l1

  SUBROUTINE reallocate_l2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    LOGICAL,ALLOCATABLE,INTENT(INOUT)::p(:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new

    LOGICAL,ALLOCATABLE::work(:,:)
    LOGICAL,PARAMETER::p0=.FALSE.
#include "reallocate_2.F90"
  END subroutine reallocate_l2

  SUBROUTINE reallocate_l3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
    LOGICAL,ALLOCATABLE,INTENT(INOUT)::p(:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new

    LOGICAL,ALLOCATABLE::work(:,:,:)
    LOGICAL,PARAMETER::p0=.FALSE.
#include "reallocate_3.F90"
  END subroutine reallocate_l3

  SUBROUTINE reallocate_l4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new)
    LOGICAL,ALLOCATABLE,INTENT(INOUT)::p(:,:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new

    LOGICAL,ALLOCATABLE::work(:,:,:,:)
    LOGICAL,PARAMETER::p0=.FALSE.
#include "reallocate_4.F90"
  END subroutine reallocate_l4

  SUBROUTINE reallocate_l5(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new)
    LOGICAL,ALLOCATABLE,INTENT(INOUT)::p(:,:,:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new

    LOGICAL,ALLOCATABLE::work(:,:,:,:,:)
    LOGICAL,PARAMETER::p0=.FALSE.
#include "reallocate_5.F90"
  END subroutine reallocate_l5

  SUBROUTINE reallocate_i1(p,lb1_new,ub1_new)
    INTEGER,ALLOCATABLE,INTENT(INOUT)::p(:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new

    INTEGER,ALLOCATABLE::work(:)
    INTEGER,PARAMETER::p0=0
#include "reallocate_1.F90"
  END subroutine reallocate_i1

  SUBROUTINE reallocate_i2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    INTEGER,ALLOCATABLE,INTENT(INOUT)::p(:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new

    INTEGER,ALLOCATABLE::work(:,:)
    INTEGER,PARAMETER::p0=0
#include "reallocate_2.F90"
  END subroutine reallocate_i2

  SUBROUTINE reallocate_i3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
    INTEGER,ALLOCATABLE,INTENT(INOUT)::p(:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new

    INTEGER,ALLOCATABLE::work(:,:,:)
    INTEGER,PARAMETER::p0=0
#include "reallocate_3.F90"
  END subroutine reallocate_i3

  SUBROUTINE reallocate_i4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new)
    INTEGER,ALLOCATABLE,INTENT(INOUT)::p(:,:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new

    INTEGER,ALLOCATABLE::work(:,:,:,:)
    INTEGER,PARAMETER::p0=0
#include "reallocate_4.F90"
  END subroutine reallocate_i4

  SUBROUTINE reallocate_i5(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new)
    INTEGER,ALLOCATABLE,INTENT(INOUT)::p(:,:,:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new

    INTEGER,ALLOCATABLE::work(:,:,:,:,:)
    INTEGER,PARAMETER::p0=0
#include "reallocate_5.F90"
  END subroutine reallocate_i5

  SUBROUTINE reallocate_r1(p,lb1_new,ub1_new)
    real,ALLOCATABLE,INTENT(INOUT)::p(:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new

    REAL,ALLOCATABLE::work(:)
    REAL,PARAMETER::p0=0.
#include "reallocate_1.F90"
  END subroutine reallocate_r1

  SUBROUTINE reallocate_r2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    real,ALLOCATABLE,INTENT(INOUT)::p(:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new

    REAL,ALLOCATABLE::work(:,:)
    REAL,PARAMETER::p0=0.
#include "reallocate_2.F90"
  END subroutine reallocate_r2

  SUBROUTINE reallocate_r3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
    real,ALLOCATABLE,INTENT(INOUT)::p(:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new

    REAL,ALLOCATABLE::work(:,:,:)
    REAL,PARAMETER::p0=0.
#include "reallocate_3.F90"
  END subroutine reallocate_r3

  SUBROUTINE reallocate_r4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new)
    real,ALLOCATABLE,INTENT(INOUT)::p(:,:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new

    REAL,ALLOCATABLE::work(:,:,:,:)
    REAL,PARAMETER::p0=0.
#include "reallocate_4.F90"
  END subroutine reallocate_r4

  SUBROUTINE reallocate_r5(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new)
    real,ALLOCATABLE,INTENT(INOUT)::p(:,:,:,:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new

    REAL,ALLOCATABLE::work(:,:,:,:,:)
    REAL,PARAMETER::p0=0.
#include "reallocate_5.F90"
  END subroutine reallocate_r5

  SUBROUTINE reallocate_c1(p,lb1_new,ub1_new)
    CHARACTER(LEN=*),ALLOCATABLE,INTENT(INOUT)::p(:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new

    CHARACTER(LEN=len(p(1))),ALLOCATABLE::work(:)
    CHARACTER(LEN=*),PARAMETER::p0=''
#include "reallocate_1.F90"
  END subroutine reallocate_c1

  SUBROUTINE reallocate_c2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    CHARACTER(LEN=*),ALLOCATABLE,INTENT(INOUT)::p(:,:)
    INTEGER,INTENT(IN)::lb1_new,ub1_new,lb2_new,ub2_new

    CHARACTER(LEN=len(p(1,1))),ALLOCATABLE::work(:,:)
    CHARACTER(LEN=*),PARAMETER::p0=''
#include "reallocate_2.F90"
  END subroutine reallocate_c2

  SUBROUTINE insert_to_array_l1(p,val,pos)
    logical,ALLOCATABLE,INTENT(INOUT)::p(:)
    logical,INTENT(IN)::val
    integer,intent(in)::pos
#include "insert_1.F90"
  END SUBROUTINE insert_to_array_l1

  SUBROUTINE insert_to_array_c1(p,val,pos)
    character(LEN=*),ALLOCATABLE,INTENT(INOUT)::p(:)
    character(LEN=*),INTENT(IN)::val
    integer,intent(in)::pos
#include "insert_1.F90"
  END SUBROUTINE insert_to_array_c1
end module util_memory
