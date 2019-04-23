module util_math
  use var_type,only:dp
  use const_math,only:onepi
  implicit none
  private
  public::cone_angle,erfunc,mbessel,polint,spline,splint,update_average,store_block_average,calculate_statistics,cross_product,normalize_vector

  interface store_block_average
     module procedure store_block_average_i,store_block_average_r
  end interface

contains
!****************************************************************
!> \brief takes two unit vectors in spherical coordinates and computes
!> the angle between them.
!>
!> \param thetaone is the angle between vector one and the z-axis
!> \param phione is the angle between vector one and the x-axis
!> \param thetatwo is the angle between vector two and the z-axis
!> \param phitwo is the angle between vector two and the x-axis
!> \param angle is the angle between the two vectors
!> \remark x = r sin (theta) cos (phi)
!> y = r sin (theta) sin (phi)
!> z = r cos (theta)
!> \author M.G. Martin
!> \date 2-4-98
!****************************************************************
  function cone_angle( thetaone, phione, thetatwo, phitwo ) result(angle)
    real,intent(in)::thetaone,thetatwo,phione,phitwo
    real::angle
    real::sintheone,costheone,sinthetwo,costhetwo,sinphione,cosphione,sinphitwo,cosphitwo,cosangle

    sintheone = sin(thetaone)
    costheone = cos(thetaone)
    sinthetwo = sin(thetatwo)
    costhetwo = cos(thetatwo)
    sinphione = sin(phione)
    cosphione = cos(phione)
    sinphitwo = sin(phitwo)
    cosphitwo = cos(phitwo)

    cosangle = sintheone*cosphione*sinthetwo*cosphitwo + sintheone*sinphione*sinthetwo*sinphitwo + costheone*costhetwo

    angle = acos(cosangle)

    return
  end function cone_angle

!> \brief complementary error function
  elemental function erfunc(x)
    real::erfunc
    real,intent(in)::x

#ifdef __USEOWN__
    real,parameter::p=0.3275911E0_dp,a1=0.254829592E0_dp,a2=-0.284496736E0_dp,a3=1.421413741E0_dp,a4=-1.453152027E0_dp,a5=1.061405429E0_dp
    real::tt,eee

    eee = exp(-x*x)
    tt = 1.0E0_dp/(1.0E0_dp + p*x)
    erfunc = ((((a5*tt+a4)*tt+a3)*tt+a2)*tt+a1)*tt*eee
#else
    erfunc = erfc(x)
#endif
    return
  end function erfunc

  elemental function mbessel(z,nu)
    real::mbessel
    real,intent(in)::z,nu

! simple form
    mbessel = sqrt(onepi/(2.0E0_dp*z))*exp(-z)
!     &         *(1.0E0_dp + (4.0E0_dp*nu**2-1)/(8.0E0_dp*z) +
!     &         (4.0E0_dp*nu**2-1)*(4.0E0_dp*nu**2-9.0E0_dp)/(2.0E0_dp*64.0E0_dp*z**2))
  end function mbessel
  
  !< \Brief normalize vector a
  !< \param a: input vector
  function normalize_vector(a)
      real, dimension(3), intent(in) :: a
      real, dimension(3) :: normalize_vector
      real :: norm_factor
      norm_factor = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      normalize_vector(1) = a(1) / norm_factor
      normalize_vector(2) = a(2) / norm_factor
      normalize_vector(3) = a(3) / norm_factor
  end function normalize_vector

  !< \Brief calculate the normalized cross product of vector a and b
  !< \param a, b input vectors a and b, with dimension of 3
  function cross_product(a, b)
      real, dimension(3) :: cross_product
      real, dimension(3), intent(in) :: a, b
      cross_product(1) = a(2) * b(3) - a(3) * b(2)
      cross_product(2) = a(3) * b(1) - a(1) * b(3)
      cross_product(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

!> \copyright (C) Copr. 1986-92 Numerical Recipes Software +3Y.
  pure subroutine polint(xa,ya,n,x,y)
    integer,intent(in)::n
    real,intent(in)::xa(n),ya(n),x
    real,intent(out)::y

    integer,parameter::nmax=10
    integer::i,m,ns
    real::dy,den,dif,dift,ho,hp,w,c(nmax),d(nmax)
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       end if
       c(i)=ya(i)
       d(i)=ya(i)
    end do
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
    return
  end subroutine polint

!> \brief Set up cubic spline derivative array
!>
!>   Given arrays \a x(1:n) and \a y(1:n) containing a tabulated function, i.e.,
!> y(i) = f(x(i)) with x ascending in order, and given values \a yp1 and \a ypn
!> for the first derivative of the interpolating function at the point 1
!> and \a n, respectively, this routine returns an array \a y2(1:n) of length \a n
!> which contains the second derivatives of the interpolating function at
!> the tabulated points \a x(i). If \a yp1 and/or \a ypn are equal to 1E30_dp or larger,
!> the routine is signaled to set the corresponding boundary condition for
!> a natural spline, with zero second derivative on that boundary
!>
!> \copyright  Numerical recipes, 2nd ed, 1992.
  subroutine spline(x,y,n,yp1,ypn,y2)
    integer::n,NMAX !< NMAX is the largest anticipated value of n
    real::yp1,ypn,x(n),y(n),y2(n)
    parameter (NMAX=500)

    integer::i,k
    real::p,qn,sig,un,u(NMAX)

    if (yp1.gt.0.99E30_dp) then
       y2(1) = 0.0E0_dp
       u(1) = 0.0E0_dp
    else
       y2(1) = -0.5E0_dp
       u(1) = (3.0E0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if
    do i = 2,n-1
       sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = sig*y2(i-1)+2.0E0_dp
       y2(i) = (sig-1.0E0_dp)/p
       u(i) = (6.0E0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn .gt. 0.99E30_dp) then
       qn = 0.0E0_dp
       un = 0.0E0_dp
    else
       qn = 0.5E0_dp
       un = (3.0E0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0E0_dp)
    do k = n-1, 1, -1
       y2(k) = y2(k)*y2(k+1)+u(k)
    end do
    return
  end subroutine spline

!> \brief Spline interpolation
!>
!>   Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
!> function (with the xa(i) in order), and given the array y2a(1:n), which
!> is from the ouput from spline, and given a value of x, this routine
!> returns a cubic-spline interpolated value y.
  subroutine splint(xa,ya,y2a,n,x,y)
    use util_runtime,only:err_exit
    use util_search,only:locate
    integer::n
    real::x,y,xa(n),y2a(n),ya(n)
    integer::khi,klo
    real::a,b,h

    klo = locate(xa,n,x,2)
    khi = klo+1
    h = xa(khi)-xa(klo)
    if (h.eq.0.) call err_exit(__FILE__,__LINE__,'bad xa input in splint',-1)
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0E0_dp
    return
  end subroutine splint

  elemental subroutine update_average(aval,val,count)
    real,intent(inout)::aval
    real,intent(in)::val
    integer,intent(in)::count
    real::inv_count

    if (count.gt.0) then
       inv_count=1.0_dp/real(count,dp)
       ! The sequence of calculations is for avoiding overflow problems
       aval=inv_count*(count-1)*aval+inv_count*val
    end if
  end subroutine update_average

!> \brief Calculate and store block averages
  elemental subroutine store_block_average_i(block_average,new_value,new_count,last_value,last_count)
    real,intent(out)::block_average
    real,intent(in)::new_value
    integer,intent(in)::new_count
    real,intent(inout)::last_value
    integer,intent(inout)::last_count
#include "store_block_average.F90"
  end subroutine store_block_average_i

  elemental subroutine store_block_average_r(block_average,new_value,new_count,last_value,last_count)
    real,intent(out)::block_average
    real,intent(in)::new_value,new_count
    real,intent(inout)::last_value,last_count
#include "store_block_average.F90"
  end subroutine store_block_average_r

  pure subroutine calculate_statistics(block_values,mean,stdev,sterr)
    real,intent(in)::block_values(:)
    real,intent(out)::mean,stdev,sterr
    integer::nblock

    nblock=size(block_values)
    mean=sum(block_values)/nblock
    stdev=sqrt(sum((block_values-mean)**2)/nblock)
    if (nblock.gt.1) sterr=stdev/sqrt(real(nblock-1,dp))
  end subroutine calculate_statistics



!      subroutine coordinate_transform(x,y,z,invh,sx,sy,sz)
!      real,intent(in)::x,y,z,invhmat
!      real,intent(out)::sx,sy,sz

!      end subroutine coordinate_transform
end module util_math
