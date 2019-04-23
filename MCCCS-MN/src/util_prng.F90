! *****************************************************************************
!> \brief Parallel (pseudo)random number generator (PRNG) for multiple streams
!>      and substreams of random numbers.
!>
!>      In detail, this RNG provides 2**64 random number streams each with a
!>      length of 2**127 resulting in a length of 2**191 for the total RNG.
!>      Moreover, each stream is divided in 2**51 substream of length 2**76.
!>      The stream lengths refer to the default precision of 32 bit random
!>      number, but also an extended precision of 53 bit per random number
!>      can be requested. In this case, two 32 bit random numbers are used
!>      to generate a 53 bit random number and therefore the stream length
!>      is halved when extended precision are requested.
!>
!>      Usage hint:
!>
!>      error=create_rng_stream(rng_stream,name,...)
!>
!>      to generate the first stream. Optionally, you may define a different
!>      seed or create a stream of extended precision (53 bits). Then
!>
!>      error=create_rng_stream(next_rng_stream,name,last_rng_stream=rng_stream)
!>
!>      to create all the following RNG streams w.r.t. the previous stream.
!>      The command line
!>
!>      x = next_random_number(rng_stream)
!>
!>      will provide the next real random number x between 0 and 1 and
!>
!>      ix = next_random_number(rng_stream,low,high)
!>
!>      the next integer random number ix between low and high from stream
!>      rng_stream. The default distribution type is a uniform distribution
!>      [0,1], but also other distribution types are available (see below).
!>
!>      Finally, do not forget to delete each created RNG stream when it is
!>      no longer needed by
!>
!>      CALL delete_rng_stream(rng_stream)
!>
!>      to avoid memory leaks!
!> \par Literature
!>      P. L'Ecuyer, R. Simard, E. J. Chen, and W. D. Kelton,
!>      "An object-oriented random-number package with many long streams
!>       and substreams", Operations Research 50(6), 1073-1075 (2002)
!> \author C++ code converted to Fortran 90/95 (18.05.2005, Matthias Krack)
! *****************************************************************************
MODULE util_prng
  USE var_type,only:dp=>double_precision
  IMPLICIT NONE
  PRIVATE
  SAVE
  ! Global parameters in this module

  CHARACTER(len=*), PARAMETER :: rng_record_format = "(A40,I2,3L2,ES25.16,18F20.1)"
  INTEGER, PARAMETER          :: rng_record_length = 433

  ! Distribution types:
  INTEGER, PARAMETER          :: GAUSSIAN = 1,& !< Gaussian distribution with zero mean and unit variance
                                 UNIFORM  = 2 !< Uniform distribution [0,1] with 1/2 mean (default)

  REAL(KIND=dp), PARAMETER :: norm  = 2.328306549295727688e-10_dp,&
                              m1    = 4294967087.0_dp,&
                              m2    = 4294944443.0_dp,&
                              a12   = 1403580.0_dp,&
                              a13n  = 810728.0_dp,&
                              a21   = 527612.0_dp,&
                              a23n  = 1370589.0_dp,&
                              two17 = 131072.0_dp,&            !< 2**17
                              two53 = 9007199254740992.0_dp,&  !< 2**53
                              fact  = 5.9604644775390625e-8_dp !< 1/2**24
  REAL(KIND=dp), DIMENSION(3,2) :: nextSeed = 12345.0_dp

! *****************************************************************************
  ! Data type definitions

  !> \brief Information on a stream.
  !>
  !> The arrays \a bg, \a cg, and \a ig contain the current
  !> state of the stream, the starting state of the current substream, and the
  !> starting state of the stream. This stream generates antithetic variates
  !> if \a antithetic = .TRUE.. It also generates numbers with extended precision
  !> (53 bits, if machine follows IEEE 754 standard), if \a extended_precision =
  !> .TRUE., otherwise, numbers with 32 bits precision are generated.
  TYPE rng_stream_type
    PRIVATE
    CHARACTER(LEN=40)             :: name
    INTEGER                       :: distribution_type
    REAL(KIND=dp), DIMENSION(3,2) :: bg,cg,ig
    LOGICAL                       :: antithetic,extended_precision
    ! only used for distribution type GAUSSIAN
    REAL(KIND=dp)                 :: buffer
    LOGICAL                       :: buffer_filled
  END TYPE rng_stream_type

  TYPE rng_stream_type_ptr
     TYPE(rng_stream_type),POINTER::val
  END TYPE rng_stream_type_ptr

! *****************************************************************************
  ! The following are the transition matrices of the two MRG components
  ! (in matrix form), raised to the powers -1, 1, 2**76, and 2**127, resp.

  REAL(KIND=dp), DIMENSION(3,3), parameter :: inv_a1 = reshape((/184888585.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp&
   ,1945170933.0_dp,0.0_dp,0.0_dp/),shape(inv_a1)) !< Inverse of a1p0
  REAL(KIND=dp), DIMENSION(3,3), parameter :: inv_a2 = reshape((/0.0_dp,1.0_dp,0.0_dp,360363334.0_dp,0.0_dp,1.0_dp&
   ,4225571728.0_dp,0.0_dp,0.0_dp/),shape(inv_a2)) !< Inverse of a2p0
  REAL(KIND=dp), DIMENSION(3,3), parameter :: a1p0 = reshape((/0.0_dp,0.0_dp,-810728.0_dp,1.0_dp,0.0_dp,1403580.0_dp&
   ,0.0_dp,1.0_dp,0.0_dp/),shape(a1p0)) !< Transition matrix for the first component raised to the power 2**0
  REAL(KIND=dp), DIMENSION(3,3), parameter :: a2p0 = reshape((/0.0_dp,0.0_dp,-1370589.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp&
   ,1.0_dp,527612.0_dp/),shape(a2p0)) !< Transition matrix for the second component raised to the power 2**0
  REAL(KIND=dp), DIMENSION(3,3), parameter :: a1p76 = reshape((/82758667.0_dp,3672831523.0_dp,3672091415.0_dp,1871391091.0_dp&
   ,69195019.0_dp,3528743235.0_dp,4127413238.0_dp,1871391091.0_dp,69195019.0_dp/),shape(a1p76)) !< Transition matrix for the first component raised to the power 2**76
  REAL(KIND=dp), DIMENSION(3,3), parameter :: a2p76 = reshape((/1511326704.0_dp,4292754251.0_dp,3859662829.0_dp,3759209742.0_dp&
   ,1511326704.0_dp,4292754251.0_dp,1610795712.0_dp,3889917532.0_dp,3708466080.0_dp/),shape(a2p76)) !< Transition matrix for the second component raised to the power 2**76
  REAL(KIND=dp), DIMENSION(3,3), parameter :: a1p127 = reshape((/2427906178.0_dp,226153695.0_dp,1988835001.0_dp,3580155704.0_dp&
   ,1230515664.0_dp,986791581.0_dp,949770784.0_dp,3580155704.0_dp,1230515664.0_dp/),shape(a1p127)) !< Transition matrix for the first component raised to the power 2**127
  REAL(KIND=dp), DIMENSION(3,3), parameter :: a2p127 = reshape((/1464411153.0_dp,32183930.0_dp,2824425944.0_dp,277697599.0_dp&
   ,1464411153.0_dp,32183930.0_dp,1610723613.0_dp,1022607788.0_dp,2093834863.0_dp/),shape(a2p127)) !< Transition matrix for the second component raised to the power 2**127

  ! Interfaces

  INTERFACE next_random_number
    MODULE PROCEDURE next_integer_random_number,&
                     next_real_random_number
  END INTERFACE

  INTERFACE random_numbers
    MODULE PROCEDURE random_numbers_1,&
                     random_numbers_2,&
                     random_numbers_3
  END INTERFACE

  ! Public parameters

  PUBLIC :: rng_record_length,&
            GAUSSIAN,&
            UNIFORM

  ! Public data types

  PUBLIC :: rng_stream_type, rng_stream_type_ptr

  ! Public subroutines

  PUBLIC :: advance_rng_state,&
            check_rng,&
            create_rng_stream,&
            delete_rng_stream,&
            dump_rng_stream,&
            get_rng_stream,&
            next_random_number,&
            random_numbers,&
            read_rng_stream,&
            reset_rng_stream,&
            reset_rng_substream,&
            reset_to_next_rng_substream,&
            set_rng_stream,&
            write_rng_matrices,&
            write_rng_stream

CONTAINS
  !> \brief Advance the state by n steps
  !>
  !> i.e. jump n steps forward, if n > 0, or
  !> backward if n < 0.
  !> IF e > 0, let n = 2**e + c
  !> IF e < 0, let n = -2**(-e) + c
  !> IF e = 0, let n = c
  !> \attention The use of this method is discouraged.
  SUBROUTINE advance_rng_state(rng_stream,e,c)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    INTEGER, INTENT(IN)                      :: e, c

    REAL(KIND=dp), DIMENSION(3, 2)           :: x
    REAL(KIND=dp), DIMENSION(3, 3)           :: u1, u2, v1, v2, w1, w2
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    u1 = 0.0_dp
    u2 = 0.0_dp
    v1 = 0.0_dp
    v2 = 0.0_dp
    w1 = 0.0_dp
    w2 = 0.0_dp

    IF (e > 0) THEN
       CALL mat_two_pow_mod_m(a1p0,u1,m1,e)
       CALL mat_two_pow_mod_m(a2p0,u2,m2,e)
    ELSE IF (e < 0) THEN
       CALL mat_two_pow_mod_m(inv_a1,u1,m1,-e)
       CALL mat_two_pow_mod_m(inv_a2,u2,m2,-e)
    END IF

    IF (c >= 0) THEN
       CALL mat_pow_mod_m(a1p0,v1,m1,c)
       CALL mat_pow_mod_m(a2p0,v2,m2,c)
    ELSE
       CALL mat_pow_mod_m(inv_a1,v1,m1,-c)
       CALL mat_pow_mod_m(inv_a2,v2,m2,-c)
    END IF

    IF (e == 0) THEN
       w1 = v1
       w2 = v2
    ELSE
       CALL mat_mat_mod_m(u1,v1,w1,m1)
       CALL mat_mat_mod_m(u2,v2,w2,m2)
    END IF

    x = 0.0_dp

    CALL mat_vec_mod_m(w1,rng_stream%cg(:,1),x(:,1),m1)
    CALL mat_vec_mod_m(w2,rng_stream%cg(:,2),x(:,2),m2)

    rng_stream%cg = x

  END SUBROUTINE advance_rng_state

  !> \brief Check the parallel (pseudo)random number generator (RNG).
  FUNCTION check_rng(output_unit) RESULT(error)
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL                                  :: error

    INTEGER                                  :: i, sumi
    REAL(KIND=dp)                            :: sum, sum3
    REAL(KIND=dp), DIMENSION(3, 2)           :: germe
    TYPE(rng_stream_type), POINTER           :: cantor, g1, g2, g3, galois, laplace, poisson
! -------------------------------------------------------------------------
! Test 1

    NULLIFY (g1)
    NULLIFY (g2)
    NULLIFY (g3)

    ! Create RNG test streams

    error= create_rng_stream(g1,"g1")
    error= create_rng_stream(g2,"g2")
    error= create_rng_stream(g3,"g3")

    WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
     "RESULTS OF THE (PSEUDO)RANDOM NUMBER GENERATOR TEST RUNS",&
     "Initial states of the (pseudo)random number streams (test 1):"
    CALL write_rng_stream(g1,output_unit)
    CALL write_rng_stream(g2,output_unit)
    CALL write_rng_stream(g3,output_unit)

    sum = next_random_number(g2) + next_random_number(g3)

    CALL advance_rng_state(g1,5,3)
    sum = sum + next_random_number(g1)

    CALL reset_rng_stream(g1)
    DO i=1,35
       CALL advance_rng_state(g1,0,1)
    END DO
    sum = sum + next_random_number(g1)

    CALL reset_rng_stream(g1)

    sumi = 0
    DO i=1,35
      sumi = sumi + next_random_number(g1,1,10)
    END DO
    sum = sum + sumi/100.0_dp

    sum3 = 0.0_dp
    DO i=1,100
      sum3 = sum3 + next_random_number(g3)
    END DO
    sum = sum + sum3/10.0_dp

    CALL reset_rng_stream(g3)
    DO i=1,5
      sum = sum + next_random_number(g3)
    END DO

    CALL reset_rng_stream(g3)
    DO i=1,4
      CALL reset_to_next_rng_substream(g3)
    END DO
    DO i=1,5
      sum = sum + next_random_number(g3)
    END DO

    CALL reset_rng_substream(g3)
    DO i=1,5
      sum = sum + next_random_number(g3)
    END DO

    CALL reset_to_next_rng_substream(g2)
    sum3 = 0.0_dp
    DO i=1,100000
      sum3 = sum3 + next_random_number(g2)
    END DO
    sum = sum + sum3/10000.0_dp

    error=set_rng_stream(g3,antithetic=.TRUE.)
    sum3 = 0.0_dp
    DO i=1,100000
      sum3 = sum3 + next_random_number(g3)
    END DO
    sum = sum + sum3/10000.0_dp

    WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
     "Final states of the (pseudo)random number streams (test 1):"
    CALL write_rng_stream(g1,output_unit)
    CALL write_rng_stream(g2,output_unit)
    CALL write_rng_stream(g3,output_unit)
    WRITE (UNIT=output_unit,FMT="(/,(T2,A))")&
     "This test routine should print for test 1 the number 25.342059"
    WRITE (UNIT=output_unit,FMT="(T2,A,F10.6)")&
     "The actual result of test 1 is                      ",sum

    CALL delete_rng_stream(g1)
    CALL delete_rng_stream(g2)
    CALL delete_rng_stream(g3)

    ! -------------------------------------------------------------------------
    ! Test 2

    NULLIFY (cantor)
    NULLIFY (galois)
    NULLIFY (laplace)
    NULLIFY (poisson)

    germe(:,:) = 1

    error=create_rng_stream(poisson,"Poisson",seed=germe)
    error=create_rng_stream(laplace,"Laplace")
    error=create_rng_stream(galois,"Galois")
    error=create_rng_stream(cantor,"Cantor")

    WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
     "Initial states of the (pseudo)random number streams (test 2):"
    CALL write_rng_stream(poisson,output_unit)
    CALL write_rng_stream(laplace,output_unit)
    CALL write_rng_stream(galois,output_unit)
    CALL write_rng_stream(cantor,output_unit)

    sum = sum + next_random_number(poisson) +&
                next_random_number(laplace) +&
                next_random_number(galois) +&
                next_random_number(cantor)

    CALL advance_rng_state(galois,-127,0)
    sum = sum + next_random_number(galois)

    CALL reset_to_next_rng_substream(galois)
    error=set_rng_stream(galois,extended_precision=.TRUE.)
    sum3 = 0.0_dp
    DO i=1,100000
      sum3 = sum3 + next_random_number(galois)
    END DO
    sum = sum + sum3/10000.0_dp

    error=set_rng_stream(galois,antithetic=.TRUE.)
    sum3 = 0.0_dp
    DO i=1,100000
      sum3 = sum3 + next_random_number(galois)
    END DO
    sum = sum + sum3/10000.0_dp
    error=set_rng_stream(galois,antithetic=.FALSE.)

    error=set_rng_stream(galois,extended_precision=.FALSE.)
    sum = sum + next_random_number(poisson) +&
                next_random_number(laplace) +&
                next_random_number(galois) +&
                next_random_number(cantor)

    WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
     "Final states of the (pseudo)random number streams (test 2):"
    CALL write_rng_stream(poisson,output_unit)
    CALL write_rng_stream(laplace,output_unit)
    CALL write_rng_stream(galois,output_unit)
    CALL write_rng_stream(cantor,output_unit)
    WRITE (UNIT=output_unit,FMT="(/,(T2,A))")&
     "This test routine should print for test 2 the number 39.697547"
    WRITE (UNIT=output_unit,FMT="(T2,A,F10.6)")&
     "The actual result of test 2 is                      ",sum

    CALL delete_rng_stream(cantor)
    CALL delete_rng_stream(galois)
    CALL delete_rng_stream(laplace)
    CALL delete_rng_stream(poisson)

  END FUNCTION check_rng

  !> \brief Check that the seeds are legitimate values.
  !>
  !> \return .FALSE.
  !> if legal seeds, .TRUE. otherwise.
  FUNCTION check_seed(seed) RESULT(error)
    REAL(KIND=dp), DIMENSION(3, 2), &
      INTENT(IN)                             :: seed
    LOGICAL                                  :: error

    CHARACTER(LEN=*), PARAMETER :: fmtstr = "(A,I1,A,ES23.14,A,ES23.14)"

    INTEGER                                  :: i
! -------------------------------------------------------------------------
    error=.false.
    DO i=1,3
      ! Check condition: 0 <= seed(:,1) < m1
      IF (seed(i,1) < 0.0_dp) THEN
        WRITE (UNIT=*,FMT=fmtstr) "ERROR: Seed(",i,",1) = ",seed(i,1)," < ",0.0_dp
        error=.true.
      END IF
      IF (seed(i,1) >= m1) THEN
        WRITE (UNIT=*,FMT=fmtstr) "ERROR: Seed(",i,",1) = ",seed(i,1)," >= ",m1
        error=.true.
      END IF

      ! Check condition: 0 <= seed(:,2) < m2
      IF (seed(i,2) < 0.0_dp) THEN
        WRITE (UNIT=*,FMT=fmtstr) "ERROR: Seed(",i,",2) = ",seed(i,2)," < ",0.0_dp
        error=.true.
      END IF
      IF (seed(i,2) >= m2) THEN
        WRITE (UNIT=*,FMT=fmtstr) "ERROR: Seed(",i,",2) = ",seed(i,2)," >= ",m2
        error=.true.
      END IF
    END DO

    ! Check condition: first or second seed is 0
    IF (ALL(seed(:,1) < 1.0_dp)) THEN
       WRITE(*,*) "ERROR: First 3 seeds = 0"
        error=.true.
    END IF
    IF (ALL(seed(:,2) < 1.0_dp)) THEN
       WRITE(*,*) "ERROR: Last 3 seeds = 0"
        error=.true.
    END IF

  END FUNCTION check_seed

  !> \brief Create a new RNG stream.
  !>
  !> \param last_rng_stream is used as a reference stream, if it is specified.
  !>
  !> \par Usage hint:
  !> CALL create_rng_stream(rng_stream,name,...) to generate the first stream.
  !>
  !> Then CALL create_rng_stream(next_rng_stream,name,rng_stream) to create
  !> all the following RNG streams.
  FUNCTION create_rng_stream(rng_stream,name,last_rng_stream,distribution_type,seed,antithetic,extended_precision) RESULT(error)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(rng_stream_type), OPTIONAL, POINTER :: last_rng_stream
    INTEGER, INTENT(IN), OPTIONAL            :: distribution_type
    REAL(KIND=dp), DIMENSION(3, 2), &
      INTENT(IN), OPTIONAL                   :: seed
    LOGICAL, INTENT(IN), OPTIONAL            :: antithetic, extended_precision
    LOGICAL                                  :: error
    INTEGER                                  :: istat
! -------------------------------------------------------------------------
    IF (ASSOCIATED(rng_stream)) CALL delete_rng_stream(rng_stream)

    ALLOCATE (rng_stream,STAT=istat)
    !assert(istat==0)

    rng_stream%name = name

    IF (PRESENT(seed)) THEN
       error=check_seed(seed)
       rng_stream%ig = seed
    ELSE
       error=.false.
       rng_stream%ig = nextSeed
    END IF

    CALL mat_vec_mod_m(a1p127,rng_stream%ig(:,1),nextSeed(:,1),m1)
    CALL mat_vec_mod_m(a2p127,rng_stream%ig(:,2),nextSeed(:,2),m2)

    rng_stream%cg = rng_stream%ig
    rng_stream%bg = rng_stream%ig

    IF (PRESENT(distribution_type)) THEN
       SELECT CASE (distribution_type)
       CASE (GAUSSIAN)
          rng_stream%distribution_type = GAUSSIAN
       CASE (UNIFORM)
          rng_stream%distribution_type = UNIFORM
       CASE DEFAULT
          WRITE(*,*) "Invalid distribution type specified"
          error=.TRUE.
       END SELECT
    ELSE IF (PRESENT(last_rng_stream)) THEN
       rng_stream%distribution_type = last_rng_stream%distribution_type
    ELSE
       rng_stream%distribution_type = UNIFORM
    END IF

    IF (PRESENT(antithetic)) THEN
       rng_stream%antithetic = antithetic
    ELSE IF (PRESENT(last_rng_stream)) THEN
       rng_stream%antithetic = last_rng_stream%antithetic
    ELSE
       rng_stream%antithetic = .FALSE.
    END IF

    IF (PRESENT(extended_precision)) THEN
       rng_stream%extended_precision = extended_precision
    ELSE IF (PRESENT(last_rng_stream)) THEN
       rng_stream%extended_precision = last_rng_stream%extended_precision
    ELSE
       rng_stream%extended_precision = .FALSE.
    END IF

    ! Initialize buffer for distribution type GAUSSIAN

    rng_stream%buffer = 0.0_dp
    rng_stream%buffer_filled = .FALSE.

  END FUNCTION create_rng_stream

  !> \brief Delete a random number stream.
   SUBROUTINE delete_rng_stream(rng_stream)
    TYPE(rng_stream_type), POINTER           :: rng_stream

    INTEGER                                  :: istat
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))
    DEALLOCATE (rng_stream,STAT=istat)
    !assert(istat==0)

  END SUBROUTINE delete_rng_stream

  !> \brief Dump a RNG stream as a record given as an internal file (string).
  SUBROUTINE dump_rng_stream(rng_stream,rng_record)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    CHARACTER(LEN=rng_record_length), &
      INTENT(OUT)                            :: rng_record
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    WRITE (UNIT=rng_record,FMT=rng_record_format)&
     rng_stream%name,&
     rng_stream%distribution_type,&
     rng_stream%antithetic,&
     rng_stream%extended_precision,&
     rng_stream%buffer_filled,&
     rng_stream%buffer,&
     rng_stream%cg,&
     rng_stream%bg,&
     rng_stream%ig

  END SUBROUTINE dump_rng_stream

  !> \brief Get the components of a RNG stream.
  !> \since 2009-11-04 lwalewsk: changed type of bg, cg and ig
  !>            from INTEGER, DIMENSION(3, 2) to REAL(KIND=dp), DIMENSION(3, 2)
  !> \since 2009-11-09 lwalewsk: getting of buffer and buffer_filled components added
  SUBROUTINE get_rng_stream(rng_stream,name,distribution_type,bg,cg,ig,&
                            antithetic,extended_precision,&
                            buffer,buffer_filled)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: name
    INTEGER, INTENT(OUT), OPTIONAL           :: distribution_type
    REAL(KIND=dp), DIMENSION(3, 2), &
      INTENT(OUT), OPTIONAL                  :: bg, cg, ig
    LOGICAL, INTENT(OUT), OPTIONAL           :: antithetic, extended_precision
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: buffer
    LOGICAL, INTENT(OUT), OPTIONAL           :: buffer_filled
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    IF (PRESENT(name)) name = rng_stream%name
    IF (PRESENT(distribution_type)) distribution_type = rng_stream%distribution_type
    IF (PRESENT(bg)) bg = rng_stream%bg
    IF (PRESENT(cg)) cg = rng_stream%cg
    IF (PRESENT(ig)) ig = rng_stream%ig
    IF (PRESENT(antithetic)) antithetic = rng_stream%antithetic
    IF (PRESENT(extended_precision)) extended_precision = rng_stream%extended_precision
    IF (PRESENT(buffer)) buffer = rng_stream%buffer
    IF (PRESENT(buffer_filled)) buffer_filled = rng_stream%buffer_filled

  END SUBROUTINE get_rng_stream

  !> \brief Returns c = MODULO(a*b,m).
  SUBROUTINE mat_mat_mod_m(a,b,c,m)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a, b
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: c
    REAL(KIND=dp), INTENT(IN)                :: m

    INTEGER                                  :: i
! -------------------------------------------------------------------------
    DO i=1,3
      CALL mat_vec_mod_m(a,b(:,i),c(:,i),m)
    END DO

  END SUBROUTINE mat_mat_mod_m

  !> \brief Compute matrix b = MODULO(a**n,m)
  SUBROUTINE mat_pow_mod_m(a,b,m,n)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: b
    REAL(KIND=dp), INTENT(IN)                :: m
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3, 3)           :: u, v, w
! -------------------------------------------------------------------------
! Initialize: u = v = a; b = I

    w = a

    b(1,1) = 1.0_dp
    b(2,1) = 0.0_dp
    b(3,1) = 0.0_dp
    b(1,2) = 0.0_dp
    b(2,2) = 1.0_dp
    b(3,2) = 0.0_dp
    b(1,3) = 0.0_dp
    b(2,3) = 0.0_dp
    b(3,3) = 1.0_dp

    ! Compute b = MODULO(a**n,m) using the binary decomposition of n

    i = n

    DO
      IF (MODULO(i,2) /= 0) THEN
        u = w
        v = b
        CALL mat_mat_mod_m(u,v,b,m)
      END IF
      i = i/2
      IF (i == 0) EXIT
      u = w
      v = w
      CALL mat_mat_mod_m(u,v,w,m)
    END DO

  END SUBROUTINE mat_pow_mod_m

  !> \brief Compute matrix b = MODULO(a**(2**e),m)
  SUBROUTINE mat_two_pow_mod_m(a,b,m,e)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: b
    REAL(KIND=dp), INTENT(IN)                :: m
    INTEGER, INTENT(IN)                      :: e

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3, 3)           :: u, v
! -------------------------------------------------------------------------
    b = a

    ! Compute b = MODULO(a**(2**e),m)

    DO i=1,e
      u = b
      v = b
      CALL mat_mat_mod_m(u,v,b,m)
    END DO

  END SUBROUTINE mat_two_pow_mod_m

  !> \brief Returns v = MODULO(a*s,m). Assumes that -m < s(i) < m.
  SUBROUTINE mat_vec_mod_m(a,s,v,m)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: v
    REAL(KIND=dp), INTENT(IN)                :: m

    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: a1, a2, c
! -------------------------------------------------------------------------
    v = 0.0_dp

    DO i=1,3
      DO j=1,3
        a2 = a(i,j)
        c = v(i)
        v(i) = a2*s(j) + c
        IF ((v(i) >= two53).OR.(v(i) <= -two53)) THEN
          a1 = INT(a2/two17)
          a2 = a2 - a1*two17
          v(i) = a1*s(j)
          a1 = INT(v(i)/m)
          v(i) = v(i) - a1*m
          v(i) = v(i)*two17 + a2*s(j) + c
        END IF
        a1 = INT(v(i)/m)
        v(i) = v(i) - a1*m
        IF (v(i) < 0.0_dp) v(i) = v(i) + m
      END DO
    END DO

  END SUBROUTINE mat_vec_mod_m

  !> \brief Get the next integer random number between low and high from the stream
  !> rng_stream.
  FUNCTION next_integer_random_number(rng_stream,low,high) RESULT(u)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    INTEGER, INTENT(IN)                      :: low, high
    INTEGER                                  :: u

    REAL(KIND=dp)                            :: r
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))
    !assert((rng_stream%distribution_type == UNIFORM))

    r = next_real_random_number(rng_stream)
    u = low + INT(r*REAL(high - low + 1,dp))

  END FUNCTION next_integer_random_number

  !> \brief Get the next real random number from the stream rng_stream.
  !> variance: variance of the Gaussian distribution (defaults to 1)
  FUNCTION next_real_random_number(rng_stream,variance) RESULT(u)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: variance
    REAL(KIND=dp)                            :: u

    REAL(KIND=dp)                            :: f, r, u1, u2, var
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    SELECT CASE (rng_stream%distribution_type)
    CASE (GAUSSIAN)
       IF (PRESENT(variance)) THEN
          var = variance
       ELSE
          var = 1.0_dp
       END IF
       ! take the random number from the buffer, if the buffer is filled
       IF (rng_stream%buffer_filled) THEN
          u = SQRT(var)*rng_stream%buffer
          rng_stream%buffer_filled = .FALSE.
       ELSE
          DO
             IF (rng_stream%extended_precision) THEN
                u1 = 2.0_dp*rn53(rng_stream) - 1.0_dp
                u2 = 2.0_dp*rn53(rng_stream) - 1.0_dp
             ELSE
                u1 = 2.0_dp*rn32(rng_stream) - 1.0_dp
                u2 = 2.0_dp*rn32(rng_stream) - 1.0_dp
             END IF
             r = u1*u1 + u2*u2
             IF ((r > 0.0_dp).AND.(r < 1.0_dp)) EXIT
          END DO
          ! Box-Muller transformation
          f = SQRT(-2.0_dp*LOG(r)/r)
          u = SQRT(var)*f*u1
          ! save the second random number for the next call
          rng_stream%buffer = f*u2
          rng_stream%buffer_filled = .TRUE.
       END IF
    CASE (UNIFORM)
       IF (rng_stream%extended_precision) THEN
          u = rn53(rng_stream)
       ELSE
          u = rn32(rng_stream)
       END IF
    END SELECT

  END FUNCTION next_real_random_number

  !> \brief Fill entity array with random numbers from the RNG stream rng_stream.
  SUBROUTINE random_numbers_1(array,rng_stream)
    REAL(KIND=dp), DIMENSION(:)              :: array
    TYPE(rng_stream_type), POINTER           :: rng_stream

    INTEGER                                  :: i
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    DO i=1,SIZE(array)
       array(i) = next_random_number(rng_stream)
    END DO

  END SUBROUTINE random_numbers_1

  !> \brief Fill entity array with random numbers from the RNG stream rng_stream.
  SUBROUTINE random_numbers_2(array,rng_stream)
    REAL(KIND=dp), DIMENSION(:, :)           :: array
    TYPE(rng_stream_type), POINTER           :: rng_stream

    INTEGER                                  :: i, j
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    DO j=1,SIZE(array,2)
       DO i=1,SIZE(array,1)
          array(i,j) = next_random_number(rng_stream)
       END DO
    END DO

  END SUBROUTINE random_numbers_2

  !> \brief Fill entity array with random numbers from the RNG stream rng_stream.
  SUBROUTINE random_numbers_3(array,rng_stream)
    REAL(KIND=dp), DIMENSION(:, :, :)        :: array
    TYPE(rng_stream_type), POINTER           :: rng_stream

    INTEGER                                  :: i, j, k
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    DO k=1,SIZE(array,3)
       DO j=1,SIZE(array,2)
          DO i=1,SIZE(array,1)
             array(i,j,k) = next_random_number(rng_stream)
          END DO
       END DO
    END DO

  END SUBROUTINE random_numbers_3

  !> \brief Read a RNG stream from a record given as an internal file (string).
  SUBROUTINE read_rng_stream(rng_stream,rng_record)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    CHARACTER(LEN=rng_record_length), &
      INTENT(IN)                             :: rng_record

    INTEGER                                  :: istat
! -------------------------------------------------------------------------
    IF (ASSOCIATED(rng_stream)) CALL delete_rng_stream(rng_stream)

    ALLOCATE (rng_stream,STAT=istat)
    !assert(istat==0)

    READ (UNIT=rng_record,FMT=rng_record_format)&
     rng_stream%name,&
     rng_stream%distribution_type,&
     rng_stream%antithetic,&
     rng_stream%extended_precision,&
     rng_stream%buffer_filled,&
     rng_stream%buffer,&
     rng_stream%cg,&
     rng_stream%bg,&
     rng_stream%ig

  END SUBROUTINE read_rng_stream

  !> \brief Reset a random number stream to its initial state.
  SUBROUTINE reset_rng_stream(rng_stream)
    TYPE(rng_stream_type), POINTER           :: rng_stream
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    rng_stream%cg = rng_stream%ig
    rng_stream%bg = rng_stream%ig

  END SUBROUTINE reset_rng_stream

  !> \brief Reset a random number stream to the beginning of its current substream.
  SUBROUTINE reset_rng_substream(rng_stream)
    TYPE(rng_stream_type), POINTER           :: rng_stream
! -------------------------------------------------------------------------

    rng_stream%cg = rng_stream%bg

  END SUBROUTINE reset_rng_substream

  !> \brief Reset a random number stream to the beginning of its next substream.
  SUBROUTINE reset_to_next_rng_substream(rng_stream)
    TYPE(rng_stream_type), POINTER           :: rng_stream

    REAL(KIND=dp), DIMENSION(3, 2)           :: u
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    u = 0.0_dp

    CALL mat_vec_mod_m(a1p76,rng_stream%bg(:,1),u(:,1),m1)
    CALL mat_vec_mod_m(a2p76,rng_stream%bg(:,2),u(:,2),m2)

    rng_stream%bg = u
    rng_stream%cg = u

  END SUBROUTINE reset_to_next_rng_substream

  !> \brief Generate the next random number with standard precision (32 bits).
  FUNCTION rn32(rng_stream) RESULT(u)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp)                            :: u

    INTEGER                                  :: k
    REAL(KIND=dp)                            :: p1, p2
! -------------------------------------------------------------------------
    ! Component 1

    p1 = a12*rng_stream%cg(2,1) - a13n*rng_stream%cg(1,1)
    k = INT(p1/m1)
    p1 = p1 - k*m1
    IF (p1 < 0.0_dp) p1 = p1 + m1
    rng_stream%cg(1,1) = rng_stream%cg(2,1)
    rng_stream%cg(2,1) = rng_stream%cg(3,1)
    rng_stream%cg(3,1) = p1

    ! Component 2

    p2 = a21*rng_stream%cg(3,2) - a23n*rng_stream%cg(1,2)
    k = INT(p2/m2)
    p2 = p2 - k*m2
    IF (p2 < 0.0_dp) p2 = p2 + m2
    rng_stream%cg(1,2) = rng_stream%cg(2,2)
    rng_stream%cg(2,2) = rng_stream%cg(3,2)
    rng_stream%cg(3,2) = p2

    ! Combination

    IF (p1 > p2) THEN
      u = (p1 - p2)*norm
    ELSE
      u = (p1 - p2 + m1)*norm
    END IF

    IF (rng_stream%antithetic) u = 1.0_dp - u

  END FUNCTION rn32

  !> \brief Generate the next random number with extended precision (53 bits).
  FUNCTION rn53(rng_stream) RESULT(u)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp)                            :: u
! -------------------------------------------------------------------------
    u = rn32(rng_stream)

    ! Note: rn32 returns 1 - u in the antithetic case

    IF (rng_stream%antithetic) THEN
      u = u + (rn32(rng_stream) - 1.0_dp)*fact
      IF (u < 0.0_dp) u = u + 1.0_dp
    ELSE
      u = u + rn32(rng_stream)*fact
      IF (u >= 1.0_dp) u = u - 1.0_dp
    END IF

  END FUNCTION rn53

  !> \brief Set the components of a RNG stream.
  !>
  !> \attention The manipulation of an active RNG stream is discouraged.
  !> \since 2009-11-09 lwalewsk: setting of buffer and buffer_filled components added
  FUNCTION set_rng_stream(rng_stream,name,distribution_type,bg,cg,ig,&
                            seed,antithetic,extended_precision,&
                            buffer,buffer_filled) RESULT(error)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    INTEGER, INTENT(IN), OPTIONAL            :: distribution_type
    REAL(KIND=dp), DIMENSION(3, 2), &
      INTENT(IN), OPTIONAL                   :: bg, cg, ig, seed
    LOGICAL, INTENT(IN), OPTIONAL            :: antithetic, extended_precision
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: buffer
    LOGICAL, INTENT(IN), OPTIONAL            :: buffer_filled
    LOGICAL                                  :: error
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))
    error=.FALSE.
    IF (PRESENT(name)) rng_stream%name = name
    IF (PRESENT(distribution_type)) rng_stream%distribution_type = distribution_type
    IF (PRESENT(bg)) rng_stream%bg = bg
    IF (PRESENT(cg)) rng_stream%cg = cg
    IF (PRESENT(ig)) rng_stream%ig = ig
    IF (PRESENT(seed)) THEN
       ! Sets the initial seed of the stream to seed
       ! NOTE: The use of this method is discouraged
       error=check_seed(seed)
       rng_stream%ig = seed
       rng_stream%cg = seed
       rng_stream%bg = seed
    END IF
    IF (PRESENT(antithetic)) rng_stream%antithetic = antithetic
    IF (PRESENT(extended_precision)) rng_stream%extended_precision = extended_precision
    IF (PRESENT(buffer)) rng_stream%buffer = buffer
    IF (PRESENT(buffer_filled)) rng_stream%buffer_filled = buffer_filled

  END FUNCTION set_rng_stream

  !> \brief Write the transformation matrices of the two MRG components (raised to
  !> the powers -1, 1, 2**76, and 2**127) the a logical output unit.
  SUBROUTINE write_rng_matrices(output_unit)
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(LEN=*), PARAMETER :: fmtstr = "(/,T4,A,/,/,(2X,3F14.1))"
    INTEGER                                  :: i, j
! -------------------------------------------------------------------------
! Print the transformation matrices for both components

    WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
      "TRANSFORMATION MATRICES FOR THE PARALLEL (PSEUDO)RANDOM NUMBER GENERATOR"

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A1",((a1p0(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A2",((a2p0(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A1**(2**76)",((a1p76(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A2**(2**76)",((a2p76(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A1**(2**127)",((a1p127(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A2**(2**127)",((a2p127(i,j),j=1,3),i=1,3)

  END SUBROUTINE write_rng_matrices

  !> \param write_all: if .TRUE., then print all stream informations.
  !>            (the default is .FALSE.)
  SUBROUTINE write_rng_stream(rng_stream,output_unit,write_all)
    TYPE(rng_stream_type), POINTER           :: rng_stream
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN), OPTIONAL            :: write_all

    LOGICAL                                  :: my_write_all
! -------------------------------------------------------------------------
    !assert(ASSOCIATED(rng_stream))

    IF (PRESENT(write_all)) THEN
      my_write_all = write_all
    ELSE
      my_write_all = .FALSE.
    END IF

    WRITE (UNIT=output_unit,FMT="(/,T2,A,/)")&
     "Random number stream <"//TRIM(rng_stream%name)//">:"

    SELECT CASE (rng_stream%distribution_type)
    CASE (GAUSSIAN)
       WRITE (UNIT=output_unit,FMT="(T4,A)")&
        "Distribution type: "//&
        "Normal Gaussian distribution with zero mean"
    CASE (UNIFORM)
       WRITE (UNIT=output_unit,FMT="(T4,A)")&
        "Distribution type: "//&
        "Uniform distribution [0,1] with 1/2 mean"
    END SELECT

    IF (rng_stream%antithetic) THEN
       WRITE (UNIT=output_unit,FMT="(T4,A)") "Antithetic:        yes"
    ELSE
       WRITE (UNIT=output_unit,FMT="(T4,A)") "Antithetic:        no"
    END IF

    IF (rng_stream%extended_precision) THEN
       WRITE (UNIT=output_unit,FMT="(T4,A)") "Precision:         53 Bit"
    ELSE
       WRITE (UNIT=output_unit,FMT="(T4,A)") "Precision:         32 Bit"
    END IF

    IF (my_write_all) THEN
       WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,(T4,A,3F20.1))")&
        "Initial state of the stream:",&
        "Component 1:",rng_stream%ig(:,1),&
        "Component 2:",rng_stream%ig(:,2)
       WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,(T4,A,3F20.1))")&
        "Initial state of the current substream:",&
        "Component 1:",rng_stream%bg(:,1),&
        "Component 2:",rng_stream%bg(:,2)
    END IF

    WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,(T4,A,3F20.1))")&
     "Current state of the stream:",&
     "Component 1:",rng_stream%cg(:,1),&
     "Component 2:",rng_stream%cg(:,2)

  END SUBROUTINE write_rng_stream
END MODULE util_prng
