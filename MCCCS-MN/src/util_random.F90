!> \brief (Pseudo)random number generator (RNG)
!>
!> -- Mersenne-Twister (MT19937) with initialization improved
!> \author Originally coded in C by Takuji Nishimura and Makoto Matsumoto (2002-01-26)
!> \author Translated to Fortran 77 by Tsuyoshi TADA. (2005-12-19)
!> \author Modified for use in topmon by Peng Bai. (2012-02-11)
module util_random
  use var_type,only:dp
  use const_math,only:twopi
  use util_runtime,only:err_exit
  use util_prng,only:rng_stream_type_ptr,create_rng_stream,next_random_number
  implicit none
  private
  save
  public::ranset,random,sphere

  TYPE(rng_stream_type_ptr),ALLOCATABLE::rng(:)
contains
!-----------------------------------------------------------------------
!> \brief Initialize large number (over 32-bit constant number)
!-----------------------------------------------------------------------
  subroutine mt_initln
    integer::ALLBIT_MASK
    integer::TOPBIT_MASK
    integer::UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
    integer::mag01(0:1)
    common /mt_mask1/ ALLBIT_MASK
    common /mt_mask2/ TOPBIT_MASK
    common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
    common /mt_mag01/ mag01
!C    TOPBIT_MASK = Z'80000000'
!C    ALLBIT_MASK = Z'ffffffff'
!C    UPPER_MASK  = Z'80000000'
!C    LOWER_MASK  = Z'7fffffff'
!C    MATRIX_A    = Z'9908b0df'
!C    T1_MASK     = Z'9d2c5680'
!C    T2_MASK     = Z'efc60000'
    TOPBIT_MASK=1073741824
    TOPBIT_MASK=ishft(TOPBIT_MASK,1)
    ALLBIT_MASK=2147483647
    ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
    UPPER_MASK=TOPBIT_MASK
    LOWER_MASK=2147483647
    MATRIX_A=419999967
    MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
    T1_MASK=489444992
    T1_MASK=ior(T1_MASK,TOPBIT_MASK)
    T2_MASK=1875247104
    T2_MASK=ior(T2_MASK,TOPBIT_MASK)
    mag01(0)=0
    mag01(1)=MATRIX_A
    return
  end subroutine mt_initln

!-----------------------------------------------------------------------
!> \brief Initialize mt(0:N-1) with a seed
!>
!> This subroutine should be called once before using the RNG, otherwise
!> a default seed will be used.
!>
!> \param s seed, must be non-negative
!-----------------------------------------------------------------------
  subroutine RANSET(s,nStream)
    use util_string,only:integer_to_string
    integer,intent(in)::s,nStream
    integer,parameter::N=624,DONE=123456789
    integer::ALLBIT_MASK
    integer::mti,initialized
    integer::mt(0:N-1)
    common /mt_mask1/ ALLBIT_MASK
    common /mt_state1/ mti,initialized
    common /mt_state2/ mt
    real(8)::iseeds(3,2)
    integer::i,ierr
    logical::err

    call mt_initln
    mt(0)=iand(s,ALLBIT_MASK)
    do mti=1,N-1
       mt(mti)=1812433253*ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
       mt(mti)=iand(mt(mti),ALLBIT_MASK)
    end do
    initialized=DONE

    if (nStream.gt.1) then
       if (allocated(rng)) deallocate(rng,stat=ierr)
       allocate(rng(0:nStream-1),stat=ierr)
       if (ierr.ne.0) call err_exit(__FILE__,__LINE__,'RANSET: allocation error',ierr)

       do i=0,nStream-1
          rng(i)%val=>null()
       end do

       iseeds=s
       err=create_rng_stream(rng(0)%val,integer_to_string(0),seed=iseeds)
       if (err) call err_exit(__FILE__,__LINE__,'RANSET: error creating prng stream 0',-1)

       do i=1,nStream-1
          err=create_rng_stream(rng(i)%val,integer_to_string(i),last_rng_stream=rng(i-1)%val)
          if (err) call err_exit(__FILE__,__LINE__,'RANSET: error creating prng stream '//integer_to_string(i),-1)
       end do
    end if

    return
  end subroutine RANSET

!-----------------------------------------------------------------------
!> \brief Generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
  function random(iStream)
    integer,intent(in)::iStream
    real::random
    integer,parameter::N=624,M=397,DONE=123456789
    integer::UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
    integer::mti,initialized
    integer::mt(0:N-1)
    integer::mag01(0:1)
    common /mt_state1/ mti,initialized
    common /mt_state2/ mt
    common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
    common /mt_mag01/ mag01
    integer::y,kk

    if (iStream.ge.0) then
       RANDOM=next_random_number(rng(iStream)%val)
       return
    end if

    ! iStream.lt.0

    if(initialized.ne.DONE)then
       call RANSET(21641,1)
    end if
!   First generates a random number on [0,0xffffffff]-interval
    if(mti.ge.N)then
       do kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
       end do
       do kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
       end do
       y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
       mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
       mti=0
    end if

    y=mt(mti)
    mti=mti+1

    y=ieor(y,ishft(y,-11))
    y=ieor(y,iand(ishft(y,7),T1_MASK))
    y=ieor(y,iand(ishft(y,15),T2_MASK))
    y=ieor(y,ishft(y,-18))

    if(y.lt.0.E0_dp) then
       random=dble(y)+2.E0_dp**32
    else
       random=dble(y)
    end if
! Divide by 2^32-1, for [0,1]-real-interval
!    random=random/4294967295.E0_dp
! Divide by 2^32, for [0,1)-real-interval
    random=random/4294967296.E0_dp
    return
  end function random

!> \brief Generate random vector of length 1 on unit spheres
!>
!> \see http://mathworld.wolfram.com/SpherePointPicking.html for further explanations
!> \see P349 of Allen & Tildesley: G.4 Random vectors on the surface of a sphere
  subroutine sphere(x,y,z,iStream)
    real,intent(out)::x,y,z
    integer,intent(in)::iStream

#ifdef RANDOM_SPHERE_DIRECT
    real::theta,phi
    theta=twopi*random(iStream)
    z=2.0_dp*random(iStream)-1.0_dp
    phi=acos(z)
    x=sin(phi)*cos(theta)
    y=sin(phi)*sin(theta)
#else
    integer::ii
    real::xi1,xi2,xisq
    do ii = 1,100
       xi1 = ( 2.0E0_dp * random(iStream) ) - 1.0E0_dp
       xi2 = ( 2.0E0_dp * random(iStream) ) - 1.0E0_dp
       xisq = xi1**2 + xi2**2
       if ( xisq .lt. 1.0E0_dp ) then
          x = 2.0E0_dp * xi1 * sqrt( 1.0E0_dp - xisq )
          y = 2.0E0_dp * xi2 * sqrt( 1.0E0_dp - xisq )
          z = ( 1.0E0_dp - 2.0E0_dp * xisq )
          return
       end if
    end do
    call err_exit(__FILE__,__LINE__,'exceeded 100 tries to get a vector in sphere',-1)
#endif
  end subroutine sphere

!> \brief Generate a normal deviate
!>
!> \param mu mean
!> \param sig standard deviation
!>
!> Leva's ratio-of-uniforms method; \see $7.3.9 of Numerical Recipes (3rd ed.)
!> \note The algorithms in P347 of Allen & Tildesley, G.4 Generating non-uniform
!> distributions, are inferior in either speed, accuracy, or both
  function gaussian(mu,sig,iStream) result(gau)
    real,intent(in)::mu,sig
    integer,intent(in)::iStream
    real::gau

    real::u,v,x,y,q

    do
       u=random(iStream)
       v=1.7156_dp*(random(iStream)-0.5_dp)
       x=u-0.449871_dp
       y=abs(v)+0.386595_dp
       q=x*x+y*(0.19600_dp*y-0.25472_dp*x)
       if (q.le.0.27597) exit
       if (q.le.0.27846 .and. v*v.le.-4.0_dp*log(u)*u*u) exit
    end do
    gau=mu+sig*v/u
  end function gaussian
end module util_random
