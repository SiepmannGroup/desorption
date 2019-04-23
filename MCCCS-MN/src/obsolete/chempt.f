      subroutine chempt(bsswap,imolty,ntries)

!    ********************************************************************
!    ** Ghost insertion of molecule into the boxes to compute chem pot **
!    ** using CBMC insertion techniques.  Works for linear or branched **
!    ** molecules and for anisotropic and Explicit atom                **
!    ** Written M.G. Martin  9-24-97                                   **
!    ********************************************************************
 
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'system.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'rosen.inc' 
!$$$      include 'boltzmann.inc'
!$$$      include 'poten.inc'
!$$$      include 'external.inc'
!$$$      include 'connect.inc'
!$$$      include 'inputdata.inc'

      logical::lctrl
      parameter (lctrl=.false.)

      logical::ovrlap,lterm,lbo,ltors

      integer(KIND=normal_int)::boxins,ichoi,ip,iwalk,ntries,itry,iincre
      integer(KIND=normal_int)::istt,iett
      integer(KIND=normal_int)::ibranp(20),nbranp

      integer(KIND=normal_int)::iu,iutry,iulast,iut,iub,icbu,islen,iins ,iii,ibr,ibr1,j,jj,jj2,jj3,jj4,invtry,ibox,iunit ,imolty,jmt,igrow
      real(KIND=double_precision)::v,vintra,vinter,vext,velect,vtornew ,vtordum,delen,vewald


      real(KIND=double_precision)::bsswap
      real(KIND=double_precision)::random,rdir,rbf,bsum
      real(KIND=double_precision)::waddnew

      real(KIND=double_precision)::v1insext,v1ins,w1ins ,v1insint ,volins,rho,arg,coru,v1inselc
      dimension bsswap(ntmax,4)

! --------------------------------------------------------------------

!      write(io_output,*) 'start CHEMP'

! *** store number of units in iunit and igrow ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)
! *** give i a phony number ***
      iins = nchain + 1
      moltyp(iins) = imolty

! *** give i a phony number ***
      iins = nchain + 1
      moltyp(iins) = imolty
!     give charges to phony number 
      do icbu = 1, iunit
         qqu(iins,icbu) = qelect(ntype(imolty,icbu))
      end do

      do 500 itry = 1,ntries
! ---    select a box at random
      
         if (random().lt.0.5d00) then
            boxins=1
         else
            boxins=2
         end if

!      X = 2 is # attempts into box 2 X=3 is success into box 1
!      X = 4 is success into box 2
!     bsswap is the same thing but keeps track of successful growths
         
         bsswap(imolty,boxins) = bsswap(imolty,boxins) + 1.0d0

! *** select a position of the first/starting unit at RANDOM ***
! *** and calculate the boltzmann weight                     ***
! *** for the chain to be INSERTED                           ***
         ichoi = nchoi1(imolty)
         do icbu = 1,ichoi
            rxp(1,icbu) = boxlx(boxins) * random()
            ryp(1,icbu) = boxly(boxins) * random()
            if (lpbcz) then
               rzp(1,icbu) = boxlz(boxins) * random()
            else if ( lsami .or. lmuir .or. ljoe ) then
               rzp(1,icbu) = 20*random()-10
            else
               rzp(1,icbu) = 0.0d0
            end if
         end do
         
! *** select starting unit ***
! --- always using bead 1 as the starting unit
         iutry = 1
         lbo = .true.

! --  insert the first atom
         call boltz( lbo,iins,iins,imolty,igrow,ovrlap,boxins,.true. ,0.0d0,0,ichoi)
         bnchem(boxins,imolty) = bnchem(boxins,imolty) + 1.0d0
         if ( ovrlap ) goto 500

! *** perform the walk according to the availibility of the choices ***
! *** and calculate the correct weight for the trial walk           ***

         w1ins = 0.0d0
         do ip = 1, ichoi
            w1ins = w1ins + bfac(ip)
         end do

! --- check for termination of walk ---
         if ( w1ins .lt. softlog ) goto 500

! --- select one position at random ---
         if ( ichoi .gt. 1 ) then
            rbf = w1ins * random()
            bsum = 0.0d0 
            do ip = 1, ichoi
               if ( .not. lovr(ip) ) then
                  bsum = bsum + bfac(ip)
                  if ( rbf .lt. bsum ) then
!                    --- select ip position ---
                     iwalk = ip
                     goto 180
                  end if
               end if
            end do
         else
            iwalk = 1
         end if

 180     v1ins =  vtry(iwalk)  
         v1insext = vtrext(iwalk)
         v1insint = vtrinter(iwalk)
         v1inselc = vtrelect(iwalk)
         
         rxnew(1) = rxp(1,iwalk)
         rynew(1) = ryp(1,iwalk)
         rznew(1) = rzp(1,iwalk)

!         if ( lelect(imolty) ) then
!            call qqcheck(iins,boxins,rxnew(1),rynew(1),rznew(1))
!         end if

! *** set walk conditions ***
         invtry = invib(imolty,iutry)
         if ( invtry .eq. 0 ) then
! --- Bead 1 is the only bead to be grown ---
            islen = 0
            goto 100
         else if ( invtry .eq. 1 ) then
! --- Bead 1 is an endpoint ---
            nbranp = 0
            iincre = 1

            iut = ijvib(imolty,iutry,1)
            lexist(iut) = .false.
            wsched(1) = iut
            islen = 1
         else
! --- Bead 1 is a midpoint or branchpoint ---
            nbranp = 0
! --- all branches will be regrown - starting at the lowest/highest ---
            rdir = random()
            if ( rdir .le. 0.5d0 ) then
               iincre = -1
            else
               iincre = 1
            end if
            if ( iincre .eq. -1 ) then
               write(io_output,*) 'imolty,iutry,invtry',imolty,iutry,invtry
               iut = ijvib(imolty,iutry,invtry)
               do iii = invtry-1, 1, -1
                  nbranp = nbranp + 1
                  ibranp(nbranp) = ijvib(imolty,iutry,iii)
               end do
            else
               iut = ijvib(imolty,iutry,1)
               do iii = 2, invtry
                  nbranp = nbranp + 1
                  ibranp(nbranp) = ijvib(imolty,iutry,iii)
               end do
            end if
            lexist(iut) = .false.
            wsched(1) = iut
            islen = 1
         end if

! *** calculate number of trial segments ***
! *** set LEXIST of trial segments to false ***
! *** generate growing schedule ***
         if ( iincre .eq. 1 ) then
! --- find trial segment with lowest index ---
 21         iulast = wsched(islen)
            iut = 0
            iutry = 0
            do iu = 1, invib(imolty,iulast)
               iub = ijvib(imolty,iulast,iu)
               if ( lexist(iub) ) then
                  iut = iub
                  go to 22
               end if
            end do
 22         if ( iut .ne. 0 ) then
               if ( nbranp .gt. 0 ) then
! --- compare to index of lowest branchpoint ---
                  if ( iut .gt. ibranp(1) ) then
! --- select branchpoint and remove it from list of branchpoints ---
                     iutry = ibranp(1)
                     lexist(iutry) = .false.
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = 0
                     else
                        do ibr = 2, nbranp
                           ibr1 = ibr - 1
                           ibranp(ibr1) = ibranp(ibr)
                        end do
                        ibranp(nbranp) = 0
                     end if
                     nbranp = nbranp - 1
                  else
! --- select IUT and leave branches ---
                     iutry = iut
                     lexist(iutry) = .false.
                  end if
               else
! --- select IUT ---
                  iutry = iut
                  lexist(iutry) = .false.
               end if
! --- add other connections (including non-selected IUT) to branchpoints ---
 23            do iu = 1, invib(imolty,iulast)
                  iub = ijvib(imolty,iulast,iu)
                  if ( lexist(iub) ) then
                     nbranp = nbranp + 1
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = iub
                     else
                        do ibr = nbranp-1, 1, -1
                           ibr1 = ibr + 1
                           if ( iub .lt. ibranp(ibr) ) then
                              ibranp(ibr1) = ibranp(ibr)
                              if ( ibr .eq. 1 ) ibranp(ibr) = iub
                           else
                              ibranp(ibr1) = iub
                              go to 24
                           end if
                        end do
                     end if
                  end if
 24               continue
               end do
            else
               if ( nbranp .gt. 0 ) then
! --- take lowest branchpoint ---
                  iutry = ibranp(1)
                  lexist(iutry) = .false.
                  do ibr = 2, nbranp
                     ibr1 = ibr - 1
                     ibranp(ibr1) = ibranp(ibr)
                  end do
                  nbranp = nbranp - 1
               else
! --- iulast is the last unit in the walkschedule ---
                  go to 100
               end if
            end if

            islen = islen + 1
            wsched(islen) = iutry
            go to 21
         else

! --- find trial segment with highest index ---
 31         iulast = wsched(islen)
            write(io_output,*) '*************'
            write(io_output,*) 'iulast',iulast
            write(io_output,*) 'ibranp', (ibranp(iu),iu=1,nbranp)
            iut = 0
            iutry = 0
            do iu = invib(imolty,iulast), 1, -1
               iub = ijvib(imolty,iulast,iu)
               if ( lexist(iub) ) then
                  iut = iub
                  go to 32
               end if
            end do
 32         if ( iut .ne. 0 ) then
               if ( nbranp .gt. 0 ) then
! --- compare to index of highest branchpoint ---
                  if ( iut .lt. ibranp(1) ) then
! --- select branchpoint and remove it from list of branchpoints ---
                     iutry = ibranp(1)
                     lexist(iutry) = .false.
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = 0
                     else
                        do ibr = 2, nbranp
                           ibr1 = ibr - 1
                           ibranp(ibr1) = ibranp(ibr)
                        end do
                        ibranp(nbranp) = 0
                     end if
                     nbranp = nbranp - 1
                  else
! --- select IUT and leave branches ---
                     iutry = iut
                     lexist(iutry) = .false.
                  end if
               else
! --- select IUT ---
                  iutry = iut
                  lexist(iutry) = .false.
               end if
! --- add other connections (including non-selected IUT) to branchpoints ---
 33            do iu = invib(imolty,iulast), 1, -1
                  iub = ijvib(imolty,iulast,iu)
                  if ( lexist(iub) ) then
                     nbranp = nbranp + 1
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = iub
                     else
                        do ibr = nbranp-1, 1, -1
                           ibr1 = ibr + 1
                           if ( iub .gt. ibranp(ibr) ) then
                              ibranp(ibr1) = ibranp(ibr)
                              if ( ibr .eq. 1 ) ibranp(ibr) = iub
                           else
                              ibranp(ibr1) = iub
                              go to 34
                           end if
                        end do
                     end if
                  end if
 34               continue
               end do
            else
               if ( nbranp .gt. 0 ) then
! --- take highest branchpoint ---
                  iutry = ibranp(1)
                  lexist(iutry) = .false.
                  do ibr = 2, nbranp
                     ibr1 = ibr - 1
                     ibranp(ibr1) = ibranp(ibr)
                  end do
                  nbranp = nbranp - 1
               else
! --- iulast is the last unit in the walkschedule ---
                  go to 100
               end if
            end if

            islen = islen + 1
            wsched(islen) = iutry
            go to 31

         end if

 100     continue

! - create walk schedule for bonded interactions -
         do iii = 1, igrow
            lexist(iii) = .true.
         end do
         do iii = 1, islen
            iut = wsched(iii)
            lexist(iut) = .false.
         end do
         do iii = 1, islen
            iut = wsched(iii)
            do j = 1, invib(imolty,iut)
               jj = ijvib(imolty,iut,j)
               if ( lexist(jj) ) then
                  wschvib(iut,j) = .true.
               else
                  wschvib(iut,j) = .false.
               end if
            end do
            do j = 1, inben(imolty,iut)
               jj2 = ijben2(imolty,iut,j)
               jj3 = ijben3(imolty,iut,j)
               if ( lexist(jj2) .and. lexist(jj3) ) then
                  wschben(iut,j) = .true.
               else
                  wschben(iut,j) = .false.
               end if
            end do
            do j = 1, intor(imolty,iut)
               jj2 = ijtor2(imolty,iut,j)
               jj3 = ijtor3(imolty,iut,j)
               jj4 = ijtor4(imolty,iut,j)
               if ( lexist(jj2) .and. lexist(jj3)  .and. lexist(jj4) ) then
                  wschtor(iut,j) = .true.
               else
                  wschtor(iut,j) = .false.
               end if
            end do
            lexist(iut) = .true.
         end do

! ------------------------------------------------------------------
!     --- not currently working 2-11-98
         lterm = .false.
!         call rosnbr ( iins, imolty, islen, boxins, lterm, igrow  )
! --- termination of cbmc attempt due to walk termination ---
         if ( lterm ) goto 500


! --- Begin DC-CBMC Corrections for NEW configuration
         waddnew = 1.0d0
         if (ldual .or. ((.not. lchgall) .and. lelect(imolty))) then 
!     compute the energy of the inserted molecule
!           calculate the full site-site energy
!           iii=1 new conformation
            do j=1,igrow
               rxuion(j,1) = rxnew(j)
               ryuion(j,1) = rynew(j)
               rzuion(j,1) = rznew(j)
               qquion(j,1) = qqu(iins,j)
            end do

            ibox=boxins
            nboxi(iins) = ibox

            istt=1
            iett = igrow
            call energy (iins,imolty, v, vintra,vinter,vext,velect ,vewald,1, ibox, istt, iett, .false.,ovrlap,.false. ,vtordum,.false.,.false.)
            
            if (ovrlap) goto 500
            delen = v - ( vnewintra + vnewinter + vnewext +vnewelect) - v1ins

            waddnew = waddnew*exp(-(beta*delen))

            vnewt     = vnewt + delen
            vnewinter = vinter - v1insint
            vnewext   = vext - v1insext
            vnewelect = velect - v1inselc
         end if
!     End DC-CBMC Corrections for NEW configuration

!     Begin Explicit Atom Corrections for NEW configuration
         if ( iunit .ne. igrow ) then
!        calculate the true Lennard-Jones energy for the hydrogens
!        iii=1 new conformation
            do j=1,iunit
               rxu(iins,j)=rxnew(j)
               ryu(iins,j)=rynew(j)
               rzu(iins,j)=rznew(j)
            end do
            ibox = boxins
            call explct(iins,vtornew,.false.,.false.)
            ltors = .false.
            
            do j=1,iunit
               rxuion(j,1) = rxu(iins,j)
               ryuion(j,1) = ryu(iins,j)
               rzuion(j,1) = rzu(iins,j)
               qquion(j,1) = qqu(iins,j)
            end do
! Calculate the energy of the non-backbone beads 
            istt = igrow+1
            iett = iunit
            call energy (iins,imolty,v, vintra,vinter,vext,velect ,vewald,1 ,ibox,istt,iett, .true.,ovrlap,ltors,vtordum, .true.,.false.)
            if (ovrlap) goto 500
            delen = v - vnewintra + vtornew

            waddnew = waddnew*exp(-beta*delen)

            vnewt     = vnewt + delen
            vnewintra = vintra
            vnewinter = vnewinter + vinter 
            vnewext   = vnewext + vext
            vnewtg    = vnewtg + vtornew
            vnewelect = vnewelect + velect
         end if

!     End Explicit Atom Corrections for NEW configuration

!     add in the contribution from the first bead
         vnewt     = vnewt + v1ins
         vnewinter = vnewinter + v1insint
         vnewext   = vnewext + v1insext
         vnewelect = vnewelect + v1inselc

         if (lpbcz) then
            volins=boxlx(boxins)*boxly(boxins)*boxlz(boxins)
         else
            volins=boxlx(boxins)*boxly(boxins)
         end if

         arg = w1ins * waddnew * weight * volins  / dble( ncmt(boxins,imolty)+1 )

         if (ltailc) then
            do jmt = 1, nmolty
               if ( jmt .eq. imolty ) then
                  rho = dble(ncmt(boxins,jmt)+1)/volins
               else
                  rho = dble(ncmt(boxins,jmt))/volins
               end if
               arg=arg*exp(-(beta*2.0d0*coru(imolty,jmt,rho,boxins)))
            end do
         end if

         acchem(boxins,imolty) = acchem(boxins,imolty)+arg
         bsswap(imolty,boxins+2) = bsswap(imolty,boxins+2) + 1.0d0

 500  continue


!       write(io_output,*) 'end CHEMPT'

      return
      end

