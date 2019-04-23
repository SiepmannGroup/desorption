!> \brief Computes the 2nd virial coefficient
!>
!> B2(T) = -2Pi Int 0toInf [ Exp[-beta*u(r)] -1] r^2 dr \n
!> Using the trapazoid method of numerical integration give \n
!> B2(T) = -2*Pi*stepvir* Sum(i=2,n-1)[ <Exp[-beta*u(r)]-1> ri^2 \n
!> 1/2( <Exp[-beta*u(r1)]-1> r1 + <Exp[-beta*u(rn)]-1> rn
!> \author Marcus Martin 1-15-97
!> revised by Robert Hembree 2-10-2016
!> \brief Computes the 2nd virial coefficient
!>
!> Uses a trapezoidal integration to integrate the virial coefficient of a
!> specific conformation of two molecules. The contribution to B2 is then weighted
!> according to the configurational intramolecular energy.
!> I.e. we compute the boltzmann weighted (via intramolecular energy) B2 virail
!> coefficient
!> subroutine parameters:
!>      bVirialClassical
!>      bVirialQuantum
subroutine virial(bVirialClassical,bVirialQuantum)
  use const_math,only:onepi,twopi
  use util_runtime,only:err_exit
  use sim_system
  use energy_pairwise,only:U2,type_2body,energy
  use energy_intramolecular,only:U_bonded
  implicit none
      integer::i, imolty, ii, j, jmolty, jj, ntii, ntjj , ntij,nnn,ip,itemp,iii
      real::xdiff,ydiff ,zdiff,dvircm
      real::binvir
      real::mass_t, factor,corr,vprev,deri_u
      real::bVirialClassical,bVirialQuantum
      real::Uintramol1,Uintramol2,vvbend,vvib,vvtors,uIntermolecular
      real::xcmsep,ycmsep,zcmsep,deviation,their_distance
      real::mayerterm,boltzfact,smallexpfact,fullexpfact,integralvalue
      real::storefirstval
      real::vEnergy(nEnergy),vmol1(nEnergy),vmol2(nEnergy)
      logical::olp=.false., firstval=.true.
      integer::iunit
#ifdef __DEBUG__
      write(io_output,*) 'start VIRIAL in ',myid
#endif

      firstval=.true.
      ! normally the algorithms describe first computing the intramolecular
      ! energy. This is not needed here as the correct weighting is already
      ! enforced by cbmc. We do however need to calculate the intramolecular
      ! portion of the electrostatic interactions. The energy subroutine makes
      ! no distinction between inter- and intra- molecular electrostatics. By
      ! computing them here in separate boxes we can later subtract out the
      ! contribution of the electrostatic energies.
      call energy(2,moltyp(2),vmol2,1,nboxi(2),1,nunit(2),.true.,olp,.false.,.false.,.false.,.false.)

      if ( nboxi(1) .eq. nboxi(2) ) then
         !write(io_output,*) 'particles found in same box'
         call err_exit(__FILE__,__LINE__,'',myid+1)
      end if
      olp=.false.
      ! calculate the differences in their COM in each direction
      xdiff = xcm(2) - xcm(1)
      ydiff = ycm(2) - ycm(1)
      zdiff = zcm(2) - zcm(1)
      imolty = moltyp(1)
      jmolty = moltyp(2)
      iunit = nunit(imolty)
      ! old code retained for future development
      mass_t = 0.0E0_dp
      do ii = 1, iunit
         mass_t = mass_t + mass(ntype(imolty,ii))
      end do
      mass_t = mass_t/1000E0_dp
      factor = -(6.6260755E-34_dp)**2*6.0221367E23_dp*1E20_dp /  (24.0E0_dp*onepi*mass_t*1.380658E-23_dp*twopi)
      deviation = starvir
      integralvalue=0.0E0_dp

      do while (deviation > rmin)
        !set the "trial" location of the chain2
        do jj=1,nunit(jmolty)
            ! center them ontop of each other in the y- and z- directions.
            ! Slowly translate molecule 2 down the x axis until we reach begin
            ! to overlap. If we overlap then the energy will be infinity.
            rxuion(jj,2) = rxu(2,jj)-xdiff+deviation
            ryuion(jj,2) = ryu(2,jj)-ydiff
            rzuion(jj,2) = rzu(2,jj)-zdiff
        end do
        ! compute the energy of the system as if the second molecule was in hte
        ! same box as the first molecule at a distance of deviation away.
        call energy(2,jmolty,vEnergy,2,nboxi(1),1,nunit(jmolty),.false.,olp,.false.,.false.,.false.,.false.)
        uIntermolecular = vEnergy(ivInterLJ)+vEnergy(ivTail)+vEnergy(ivElect)+vEnergy(iv3body)+vEnergy(ivFlucq)
        ! subtract out the intramolecular contrubutions to the electric
        ! energies that were calculated before leaving only the
        ! intermolecular contributions. Only the intramolecular components
        ! for mol2 were calculated here so only those need to be removed.
        uIntermolecular = uIntermolecular-vmol2(ivElect)
        ! lets do the soft cut stuff
        smallexpfact = -uIntermolecular/virtemp
        if(smallexpfact < -2.3E0_dp*softcut) then
            mayerterm=0.0
        else if (smallexpfact > softcut.or.olp) then
            ! essentially these are states where the energy appears to be
            ! infinite.
            mayerterm = 0.0E0_dp
        else
            mayerterm = exp(-uIntermolecular/virtemp)
        end if
        ! the 2* accounts for the double counting of terms in the trapezoidal
        ! method. These are later removed by dividing by 2.
        integralvalue=integralvalue+2*(1-mayerterm)*deviation**2
        if(firstval) then
          !stores the first term because the first and last terms have half
          !value of the rest of the terms in the trapezoidal integration scheme.
          firstval = .false.
          storefirstval = integralvalue/2.0
        end if
        deviation = deviation-stepvir
      end do

      ! handles the end points of the trapezoidal rule which have half of the
      ! value of all of the other points.
      integralvalue = integralvalue-storefirstval
      integralvalue = integralvalue-(1.0-mayerterm)*(deviation+stepvir)**2
      ! now multiply by the constant terms in the integral... 2pi*dr and divide
      ! by two for the trapezoidal rule
      integralvalue = integralvalue*twopi*stepvir/2.0

      bVirialClassical = bVirialClassical+integralvalue

#ifdef __DEBUG__
      write(io_output,*) 'end VIRIAL in ',myid
#endif

end subroutine virial
