!> \brief Sami's parameters: to be used with lsami = .true. AND lmuir = .true.
MODULE energy_sami
  use sim_system
  implicit none
  private
  save
  public::susami,ljsami,exsami,sumuir,ljmuir,exmuir

! EXTERNALMUIR.INC
  real::sigpri,c9ch2,c3ch2,c9ch3,c3ch3,zprmin,v2prmin,v3prmin,beta1,beta2,betac2,betac3
! LJSAMIPARA.INC
  real::alpha1,alpha2
  real::sij(9),eij(9),vsh(9),vsha(9)
  integer::tau1,tau2
  parameter ( alpha1=21.162E0_dp, alpha2=-21.162E0_dp, beta1=661.6E0_dp, beta2=6616.0E0_dp, tau1=-32, tau2=-16 )
contains
!> \brief Set-up SAMI's potentials
  subroutine susami()
    real::hsig,heps,tsig,teps
    parameter ( hsig=4.220E0_dp, heps=110.68816E0_dp,  tsig=3.527E0_dp, teps=79.982210E0_dp )
    real::rcsami
    integer::ij
! --------------------------------------------------------------------
    rcsami = 2.5E0_dp*tsig
    rcut(1) = rcsami
    sij(1)=hsig
    sij(2)=0.5E0_dp*(hsig+tsig)
    sij(3)=sij(2)
    sij(4)=sij(2)
    sij(5)=tsig
    sij(6)=tsig
    sij(7)=sij(2)
    sij(8)=tsig
    sij(9)=tsig

    eij(1)=heps
    eij(2)=sqrt(heps*teps)
    eij(3)=eij(2)
    eij(4)=eij(2)
    eij(5)=teps
    eij(6)=teps
    eij(7)=eij(2)
    eij(8)=teps
    eij(9)=teps

    vsh(1)  = eij(1) * ( ( 13.0E0_dp * (sij(1)/rcsami)**12 ) + (  4.0E0_dp * (sij(1)/rcsami)**3  ) )
    vsha(1) = eij(1) * ( ( 12.0E0_dp * sij(1)**12 / rcsami**13 ) + (  3.0E0_dp * sij(1)**3  / rcsami**4  ) )

    do ij = 2, 9
       vsh(ij)  = 4.0E0_dp * eij(ij) *  ( ( 13.0E0_dp * (sij(ij)/rcsami)**12 ) - (  7.0E0_dp * (sij(ij)/rcsami)**6  ) )
       vsha(ij) = 4.0E0_dp * eij(ij) * ( ( 12.0E0_dp * sij(ij)**12 / rcsami**13 ) - (  6.0E0_dp * sij(ij)**6  / rcsami**7  ) )
    end do

    ! do ij = 1,9
    !    write(io_output,*) 'ij',ij,'vsh',(vsh(ij)/80.0E0_dp),'vsha',(vsha(ij)/80.0E0_dp),'eij',(eij(ij)/80.0E0_dp)
    ! end do
  end subroutine susami

!> \brief Calculates SAMI's LJ and 12+3 energy for a bead
  function ljsami( rijsq, ntij )
    real::ljsami, rijsq, rij, sr
    integer::ntij
! --------------------------------------------------------------------
    rij = sqrt( rijsq )
    sr = sij(ntij) / rij

    if ( ntij .eq. 1 ) then
       ! head-head interaction ( repulsive 12+3 interaction )
       ljsami = ( eij(1) * sr**3 * ( 1.0E0_dp + sr**9 ) ) - vsh(1) + ( rij * vsha(1) )
    else
       ! head-tail or tail-tail interaction ( LJ 12-6 interaction )
       ljsami = ( 4.0E0_dp * eij(ntij) * sr**6 * ( sr**6 - 1.0E0_dp ) ) - vsh(ntij) + ( rij * vsha(ntij) )
    end if
  end function ljsami

!> \brief Calculates the SAMI's external energy for a bead
  function exsami( z, ntj )
    real::exsami, z
    integer::ntj
! --------------------------------------------------------------------
    if ( ntj .eq. 1 ) then
       ! HEADgroup potential ---
       if ( z .le. alpha2 ) then
          exsami = 0.0E0_dp
       else
          exsami = beta2 / ( 1.0E0_dp + ( (z/alpha2) - 1.0E0_dp )**tau2 )
       end if
    else
       ! TAILgroup potential ---
       if ( z .ge. alpha1 ) then
          exsami = 0.0E0_dp
       else
          exsami = beta1 / ( 1.0E0_dp + ( 1.0E0_dp - (z/alpha1) )**tau1 )
       end if
    end if
  end function exsami

!> \brief Calculate constants for lmuir external potential
  subroutine sumuir()
    sigpri = 0.715E0_dp * sqrt( 3.8E0_dp * 3.93E0_dp )
    c9ch2 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*47.0E0_dp) ) * sigpri**9
    c3ch2 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*47.0E0_dp) ) * sigpri**3
    c9ch3 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*114.0E0_dp) ) * sigpri**9
    c3ch3 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*114.0E0_dp) ) * sigpri**3
    zprmin = ( 3.0E0_dp**(1/6.0E0_dp) ) * sigpri
    v2prmin = c9ch2 / zprmin**9 - c3ch2 / zprmin**3
    v3prmin = c9ch3 / zprmin**9 - c3ch3 / zprmin**3
    betac2 = beta1 - v2prmin
    betac3 = beta1 - v3prmin
  end subroutine sumuir

!> \brief Calculates SAMI's 12+3 energy for headgroup
!> and normal LJ energy for tail
  function ljmuir ( rijsq, ntij )
    real::ljmuir, rijsq, sr, sr2, sr6, epshead, sighead
    integer::ntij
    ! attention: eps_hh / 4 used, since later multiplied by 4 ---
    ! parameter (epshead=27.67204E0_dp,sighead=4.22E0_dp)
    parameter (epshead=27.7204E0_dp,sighead=6.5E0_dp)
! --------------------------------------------------------------------
    if ( ntij .eq. 1 ) then
       sr = sighead / sqrt( rijsq )
       ! write(io_output,*) 'sr',sr,'v',4.0E0_dp*epshead*sr**3*(sr**9+1.0E0_dp)
       ljmuir = 4.0_dp*epshead * sr**3 * ( sr**9 + 1.0E0_dp )
    else
       sr2 = vvdW(3,ntij) / rijsq
       sr6 = sr2 * sr2 * sr2
       ljmuir = vvdW(1,ntij) * sr6 * ( sr6 - 1.0E0_dp)
    end if
  end function ljmuir

!> \brief calculates the lmuir external energy for a bead
  function exmuir( z, ntj )
    real::exmuir, z
    integer::ntj
! --------------------------------------------------------------------
    if ( ntj .eq. 1 ) then
       ! HEADgroup potential ---
       if ( z .le. alpha2 ) then
          exmuir = 0.0E0_dp
       else
          exmuir = beta2 / ( 1.0E0_dp + ( (z/alpha2) - 1.0E0_dp )**tau2 )
       end if
    else
       ! TAILgroup potential ---
       if ( z .ge. alpha1 ) then
          exmuir = 0.0E0_dp
       else
          if ( z .lt. zprmin ) then
             if ( ntj .eq. 2 ) then
                exmuir = betac2 / (1.0E0_dp+(1.0E0_dp-(z/alpha1))**tau1 ) + v2prmin
             else
                exmuir = betac3 / (1.0E0_dp+(1.0E0_dp-(z/alpha1))**tau1 ) + v3prmin
             end if
          else
             if ( ntj .eq. 2 ) then
                exmuir = betac2 / (1.0E0_dp+(1.0E0_dp-(z/alpha1))**tau1 ) + c9ch2 / z**9 -  c3ch2 / z**3
             else
                exmuir = betac3 / (1.0E0_dp+(1.0E0_dp-(z/alpha1))**tau1 ) + c9ch3 / z**9 -  c3ch3 / z**3
             end if
          end if
       end if
    end if
  end function exmuir
end MODULE energy_sami
