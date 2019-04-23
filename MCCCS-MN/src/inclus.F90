subroutine inclus(inclnum,inclmol,inclbead,inclsign,ncarbon,ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)
  use util_runtime,only:err_exit
  use sim_system
  implicit none

  integer::inclnum,inclmol(ntmax*numax*numax),inclbead(ntmax*numax*numax,2),inclsign(ntmax*numax*numax),ncarbon(ntmax),ainclnum,ainclmol(ntmax*numax*numax),ainclbead(ntmax*numax*numax,2),a15t(ntmax*numax*numax)
  ! variables added (3/24/05) for variable 1-4 interactions
  real::ofscale(ntmax*numax*numax),ofscale2(ntmax*numax*numax)
  integer::imolty,m,n,nb,mb,ioffset
! ----------------------------------------------------------------
!> \bug This is modified to work only for TATB NR-2007! (Is this still true?)

  ! triple loop over all types of molecules
  do imolty = 1, nmolty
     do m = 1, nunit(imolty)
        do n = 1, nunit(imolty)
           linclu(imolty,m,n) = .true.
           lqinclu(imolty,m,n) = .true.
           ! by default, don't want any 1-5 r^12 interactions
           lainclu(imolty,m,n) = .false.
           ljscale(imolty,m,n) = 1.0
           qscale2(imolty,m,n) = 1.0
        end do
     end do

     ! double loop over all units
     do m = 1, nunit(imolty)
        ! exclude all self interactions -
        linclu(imolty,m,m) = .false.
        lqinclu(imolty,m,m) = .false.

        ! exclude all directly bonded beads (vibrations) -
        do n = 1, invib(imolty,m)
           nb = ijvib(imolty,m,n)
           linclu(imolty,m,nb) = .false.
           lqinclu(imolty,m,nb) = .false.
        end do

        !> \bug Need to check: exclude carbons around a quaternary center for explct
        if (invib(imolty,m) .eq. 4) then
           do n = 1,4
              do nb = 1,4
                 linclu(imolty,ijvib(imolty,m,n),ijvib(imolty,m,nb))=.false.
                 lqinclu(imolty,ijvib(imolty,m,n),ijvib(imolty,m,nb))=.false.
              end do
           end do
        end if

        ! exclude all next-nearest neighbor bonded beads (bending)
        do n = 1, inben(imolty,m)
           nb = ijben3(imolty,m,n)
           linclu(imolty,m,nb) = .false.
           lqinclu(imolty,m,nb) = .false.
        end do

        ! exclude all third-nearest neighbor bonded beads (torsions)
        do n = 1, intor(imolty,m)
           nb = ijtor4(imolty,m,n)
           linclu(imolty,m,nb) = .false.
           ! don't set lqinclu since we want 1-4 interactions, unless lq14scale is F
           if (.not.lq14scale(imolty)) then
              lqinclu(imolty,m,nb) = .false.
           else
              qscale2(imolty,m,nb) = qscale(imolty)
              qscale2(imolty,nb,m) = qscale(imolty)
           end if
        end do
     end do


     ! add in 1-5 interactions according to aincl
     do m = 1,ainclnum
        if ( ainclmol(m) .eq. imolty ) then
           mb = ainclbead(m,1)
           nb = ainclbead(m,2)
           lainclu(imolty,mb,nb) = .true.
           lainclu(imolty,nb,mb) = .true.
           a15type(imolty,mb,nb) = a15t(m)
           a15type(imolty,nb,mb) = a15t(m)
        end if
     end do

     !> \bug Need to check: exclude all hydrogens that have their carbons excluded
     if ( ncarbon(imolty) .lt. nunit(imolty) ) then
        if (ncarbon(imolty) .eq. 3 .and. nunit(imolty) .eq. 8) then
           ! ethane with bead 3 being hydrogen
           ioffset = 0
        else
           ioffset = 1
        end if
        do m = ncarbon(imolty)+ioffset,nunit(imolty)
           ! hydrogens only have one vibration and that is to the C atom
           mb = ijvib(imolty,m,1)
           do nb = 1, ncarbon(imolty)
              if ( .not. linclu(imolty,mb,nb) ) then
                 linclu(imolty,m,nb) = .false.
                 linclu(imolty,nb,m) = .false.
                 lqinclu(imolty,m,nb) = .false.
                 lqinclu(imolty,nb,m) = .false.
              end if
           end do
           do n=m+1,nunit(imolty)
              nb = ijvib(imolty,n,1)
              linclu(imolty,n,nb) = .false.
              linclu(imolty,nb,n) = .false.
              lqinclu(imolty,n,nb) = .false.
              lqinclu(imolty,nb,n) = .false.
              if ( .not. linclu(imolty,mb,nb) ) then
                 linclu(imolty,m,n) = .false.
                 linclu(imolty,n,m) = .false.
                 linclu(imolty,m,nb) = .false.
                 linclu(imolty,nb,m) = .false.
                 linclu(imolty,n,mb) = .false.
                 linclu(imolty,mb,n) = .false.
                 lqinclu(imolty,m,n) = .false.
                 lqinclu(imolty,n,m) = .false.
                 lqinclu(imolty,m,nb) = .false.
                 lqinclu(imolty,nb,m) = .false.
                 lqinclu(imolty,n,mb) = .false.
                 lqinclu(imolty,mb,n) = .false.
              end if
           end do
        end do
     end if

     ! exclude charge interactions if lqchg is false
     do m = 1,nunit(imolty)
        do n = m+1,nunit(imolty)
           if (.not.lqchg(ntype(imolty,m)).or..not.lqchg(ntype(imolty,n))) then
              lqinclu(imolty,m,n) = .false.
              lqinclu(imolty,n,m) = .false.
           end if
        end do
     end do


     !> \bug Removing this part for TATB (NR-2007) (Does TATB has interactions between rigid beads?)
     if (lrigid(imolty)) then
        ! dont include rigid beads

        ! there will be no intramolecular forces between rigid beads
        ! or beads connected one away from a rigid bead
        do m = 1,invib(imolty,riutry(imolty,1)) ! loop over all beads connected to grow point
           mb = ijvib(imolty,riutry(imolty,1),m)
           do n = riutry(imolty,1), nunit(imolty) ! loop over all beads higher than growpoint
              linclu(imolty,n,mb) = .false.
              linclu(imolty,mb,n) = .false.
              lqinclu(imolty,n,mb) = .false.
              lqinclu(imolty,mb,n) = .false.
           end do
        end do
        do m = riutry(imolty,1), nunit(imolty)
           do n = riutry(imolty,1), nunit(imolty) ! loop over all beads higher than grow point
              linclu(imolty,m,n) = .false.
              linclu(imolty,n,m) = .false.
              lqinclu(imolty,m,n) = .false.
              lqinclu(imolty,n,m) = .false.
           end do
        end do
     end if

     ! include or exclude additional beads accoring to incl
     do n = 1,inclnum
        if ( inclmol(n) .eq. imolty ) then
           m = inclbead(n,1)
           nb = inclbead(n,2)
           if (inclsign(n).gt.0) then
              linclu(imolty,m,nb) = .true.
              linclu(imolty,nb,m) = .true.
              lqinclu(imolty,m,nb) = .true.
              lqinclu(imolty,nb,m) = .true.
              ljscale(imolty,m,nb) = ofscale(n)
              ljscale(imolty,nb,m) = ofscale(n)
              qscale2(imolty,m,nb) = ofscale2(n)
              qscale2(imolty,nb,m) = ofscale2(n)
           else if (inclsign(n).lt.0) then
              linclu(imolty,m,nb) = .false.
              linclu(imolty,nb,m) = .false.
              lqinclu(imolty,m,nb) = .false.
              lqinclu(imolty,nb,m) = .false.
           else
              write(io_output,*) 'INCLUS: n,inclsign(n)',n,inclsign(n)
              call err_exit(__FILE__,__LINE__,'inclusign must be 1 or -1',myid+1)
           end if
        end if
     end do

     ! self consistency check
     do m = 1,nunit(imolty)
        do n = m+1,nunit(imolty)
           if (linclu(imolty,m,n) .neqv. linclu(imolty,n,m)) then
              write(io_output,*) 'LJ interaction between beads',m, n,'was not consistent, turning off'
              linclu(imolty,m,n) = .false.
              linclu(imolty,n,m) = .false.
           end if

           if (lqinclu(imolty,m,n) .neqv. lqinclu(imolty,n,m)) then
              write(io_output,*) 'Charge interaction between beads',m, n,'was not consistent, turning off'
              lqinclu(imolty,m,n) = .false.
              lqinclu(imolty,n,m) = .false.
           end if
        end do
     end do
  end do
end subroutine inclus






