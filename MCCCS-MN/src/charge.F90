!> \brief Calculates the intramolcular polarization energies
!> for fluctuating charge moves
!>
!> \par History
!> written by Marcus G. Martin 9-16-97 \n
!> rewritten by Bin Chen 6-25-99
!> updated to remove vewald calculation 1-26-15
subroutine charge( i, qion, vflucq )
  use sim_system
  use energy_pairwise,only:type_2body
  implicit none

  integer::i,imolty,iunit,ii,jj,ntii,ntjj,ntij,ibox
  real::vflucq,qion(numax),qqii
  vflucq = 0.0E0_dp
  imolty = moltyp(i)
  iunit = nunit(imolty)
  ibox = nboxi(i)

! *************************************
! INTRACHAIN FLUCQ INTERACTIONS ***
! *************************************

  ! calculate intramolecular flucq energy for chain i
  do ii = 1, iunit
     ntii = ntype(imolty,ii)
     qqii = qion(ii)
     
     do jj = ii, iunit
        if ( ii .eq. jj ) then
           vflucq = vflucq + xiq(ntii)*qqii + jayself(ntii)*qqii*qqii
        else
           ntjj = ntype(imolty,jj)
           ntij = type_2body(ntii,ntjj)
           
           vflucq = vflucq + jayq(ntij)*qqii*qion(jj)
        end if
     end do
  end do
  ! remove the ground state gas phase energy
  vflucq = vflucq - fqegp(imolty)
  
  return
end subroutine charge






