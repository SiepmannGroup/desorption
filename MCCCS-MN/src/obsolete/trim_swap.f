      subroutine trim_swap(bsswap,bnswap,bnswap_in,bnswap_out,cnt_wf1, cnt_wf2,cnt_wra1,cnt_wra2)

!    ********************************************************************
!    ** removes a molecule from one box and inserts it into the other  **
!    ** using CBMC insertion techniques.  Works for linear or branched **
!    ** molecules and for DC-CBMC and Explicit atom                    **
!    ** Rewritten from old swap and swapbr subroutines by M.G. Martin  **
!    ** 9-18-97                                                        ** 
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
!$$$      include 'conver.inc'
!$$$      include 'coord2.inc'
!$$$      include 'system.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'rosen.inc' 
!$$$      include 'boltzmann.inc'
!$$$      include 'external.inc'
!$$$      include 'connect.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'poten.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'neigh.inc'
!$$$      include 'cell.inc'
      
      logical::ovrlap,lterm,lnew,lempty,ldone,ltors,lovrh,lfavor, laccept,lswapinter,lrem_out,lins_in,lneighij,linsk_in, lremk_in,lrem_clu,lins_clu,lfixnow

     
      integer(KIND=normal_int)::boxins,boxrem,imol,ichoi,ip,iwalk,idum
      integer(KIND=normal_int)::istt,iett,ncount,itype,ipair,ipairb,beg,flagon

      integer(KIND=normal_int)::iutry,icbu,ifrom,irem,iins,glist,findex ,iii,j,ibox,iunit,ic,pointp,imolty,imt,jmt,igrow ,pointp2,jins,jmolty,neighj_num,neighk_num ,joffset,koffset,kmolty,kins,target,cnt_wf1 ,cnt_wf2,neigh_old,cnt_wra1,cnt_wra2,k

      dimension glist(numax),cnt_wf1(0:6,0:6,4),cnt_wf2(0:6,0:6,4), cnt_wra1(1000,4),cnt_wra2(1000,4)

      real(KIND=double_precision)::sx,sy,sz,ddum(27)

      real(KIND=double_precision)::v,vintra,vinter,vext,velect,vtorold,vtornew ,vewald,vflucq,delen,deleo,rpair
      real(KIND=double_precision)::vnewflucq,voldflucq,qion,ctorfo,ctorfn
      dimension qion(numax)
      real(KIND=double_precision)::rxuold,ryuold,rzuold
      dimension rxuold(numax),ryuold(numax),rzuold(numax)
      real(KIND=double_precision)::bsswap,bnswap,bnswap_in,bnswap_out
      real(KIND=double_precision)::random,rmol,rbf,bsum
      real(KIND=double_precision)::waddnew,waddold

      real(KIND=double_precision)::total_NBE,vintran,velectn,vewaldn,vtgn
      real(KIND=double_precision)::vbendn,vvibn 



      real(KIND=double_precision)::v1insext,v1remext,v1ins,w1ins,v1rem,w1rem ,v1insint,v1remint,v1insewd,v1remewd ,wnlog,wolog,wdlog,wratio,vinsta,vremta ,volins,volrem,rho,arg,coru,v1inselc,v1remelc
      real(KIND=double_precision)::rvol,x,y,z,rijsq,wbias_ins,wbias_rem,r ,xi1,xi2,xisq
      dimension bsswap(ntmax,npabmax,nbxmax*2), bnswap(ntmax,npabmax,nbxmax*2),bnswap_in(ntmax,2), bnswap_out(ntmax,2)
      real(KIND=double_precision)::vrecipn,vrecipo,vdum,whins,whrem
      real(KIND=double_precision)::rxuh,ryuh,rzuh,delenh,vtrhext,vtrhintra ,vtrhinter,vtrhelect,vtrhewald,vtrhtg,bfach
     
      dimension bfach(nchmax),delenh(nchmax),vtrhinter(nchmax) ,vtrhext(nchmax),vtrhintra(nchmax),vtrhelect(nchmax) ,vtrhewald(nchmax),vtrhtg(nchmax)
     
      dimension lovrh(nchmax)
      dimension rxuh(numax,nchmax),ryuh(numax,nchmax) ,rzuh(numax,nchmax)

! --------------------------------------------------------------------

!      write(io_output,*) 'START TRIM_SWAP'
!      write(11,*) '1:',neigh_cnt(18)

      lempty = .false.
      lfixnow = .false.
      lins_in = .false.
      linsk_in = .false.
      
   
! --- select a molecule typ with probabilities given in pmswmt
      rmol = random()
      ldone = .false.
      do imol = 1, nmolty
         if ( rmol .lt. pmswmt(imol) ) then
            if ( .not. ldone ) then
               imolty = imol
               ldone = .true.
            end if
         end if
      end do
         
! ---    select a box given in pmswatyp
      if ( nswapb(imolty) .gt. 1 ) then
         rpair = random()
         do 96 ipair = 1, nswapb(imolty)
            if ( rpair .lt. pmswapb(imolty,ipair) ) then
               ipairb = ipair
               goto 97
            end if
 96      continue
      else
         ipairb = 1
      end if

 97   if (random().lt.0.5d0) then
         boxins=box1(imolty,ipairb)
         boxrem=box2(imolty,ipairb)
      else
         boxins=box2(imolty,ipairb)
         boxrem=box1(imolty,ipairb)
      end if
      if ( boxins .eq. boxrem ) then
         lswapinter = .false.
      else
         lswapinter = .true.
      end if
!      write(io_output,*) 'boxins:',boxins,'boxrem:',boxrem
      if ( .not. (lgibbs .or. lgrand) .and. lswapinter )  call cleanup('no interbox swap if not gibbs/grand ensemble!')
         
! *** select a chain in BOXREM at random ***

      if ( ncmt(boxrem,imolty) .eq. 0 ) then
         lempty = .true.
         if ( .not. lswapinter .or. lrigid(imolty) ) return
      else if ( lswapinter .or. lavbmc1(imolty) .and. .not.   (lavbmc2(imolty) .or. lavbmc3(imolty)) ) then

! *** for the advanced AVBMC algorithm, this particle will be selected in
! *** sub-regions defined by Vin

         pointp = int( dble(ncmt(boxrem,imolty))*random() ) + 1
         irem = parbox(pointp,boxrem,imolty)
         if ( moltyp(irem) .ne. imolty )  write(io_output,*) 'screwup swap, irem:',irem,moltyp(irem),imolty
         ibox = nboxi(irem)
         if ( ibox .ne. boxrem ) call cleanup('problem in swap')
      end if

!$$$      write(io_output,*) 'particle ',irem,' is being removed, imolty is:',
!$$$     &     imolty,' and the box is:',boxrem


! ===>  for both gibbs and grand-canonical we have:
! --- insert a chain in box: boxins 
! --- remove one in box: boxrem
!      write(io_output,*) 'boxrem',boxrem,' imolty',imolty,' lempty',lempty
!     bnswap(imolty,X) decoder X = 1 is # attempts into box 1
!      X = 2 is # attempts into box 2 X=3 is success into box 1
!      X = 4 is success into box 2
!     bsswap is the same thing but keeps track of successful growths

      if (.not. lempty) bnswap(imolty,ipairb,boxins)  = bnswap(imolty,ipairb,boxins) + 1.0d0
      bsswap(imolty,ipairb,boxins) = bsswap(imolty,ipairb,boxins)+1.0d0
      
!     *** store number of units in iunit ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)
!     *** give i a phony number ***
      if ( lswapinter ) then
         iins = nchain + 1
         moltyp(iins) = imolty
!     give charges to phony number 
         if ( lempty ) then
            do icbu = 1, iunit
               qqu(iins,icbu) = qelect(ntype(imolty,icbu))
            end do
         else
            do icbu = 1, iunit
               qqu(iins,icbu) = qqu(irem,icbu)
            end do
         end if
      else if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then
         iins = 0
      else
         iins = irem
      end if

! *** select a position of the first/starting unit at RANDOM ***
! *** and calculate the boltzmann weight                     ***
! *** for the chain to be INSERTED                           ***


      if (lrigid(imolty)) then
         beg = riutry(imolty,1)
      else
         beg = 1
      end if

      wbias_ins = 1.0d0
      wbias_rem = 1.0d0

      ichoi = nchoi1(imolty)

      if (lsolid(boxins) .and. .not. lrect(boxins)) then
            ibox = boxins
            do icbu = 1,ichoi
                  sx = random()
                  sy = random()
                  sz = random()
                  rxp(1,icbu) = sx*hmat(ibox,1)+sy*hmat(ibox,4) +sz*hmat(ibox,7)
                  ryp(1,icbu) = sx*hmat(ibox,2)+sy*hmat(ibox,5) +sz*hmat(ibox,8)
                  rzp(1,icbu) = sx*hmat(ibox,3)+sy*hmat(ibox,6) +sz*hmat(ibox,9)
            end do
         else
            do icbu = 1,ichoi
               rxp(1,icbu) = boxlx(boxins) * random()
               ryp(1,icbu) = boxly(boxins) * random()
               if (lpbcz .or. lslit) then
                  rzp(1,icbu) = boxlz(boxins) * random()
               else if ( lsami .or. lmuir .or. ljoe ) then
                  if ( lempty ) then
                     rzp(1,icbu) = 20*random()-10
                  else
                     rzp(1,icbu) = rzu(irem,1)
                  end if
               else
                  rzp(1,icbu) = 0.0d0
               end if
            end do
      end if


!     *** select starting unit ***
!     --- always using bead 1 as the starting unit
      iutry = beg
      lnew = .true.
      glist(1) = beg

!     --  insert the first atom

      call boltz(lnew,.true.,ovrlap,iins,iins,imolty,boxins, ichoi,idum,1,glist,0.0d0)

      if ( ovrlap ) return
      
!     *** perform the walk according to the availibility of the choices ***
!     *** and calculate the correct weight for the trial walk           ***
      
      w1ins = 0.0d0
      do ip = 1, ichoi
         w1ins = w1ins + bfac(ip)
      end do
      
!     --- check for termination of walk ---
      if ( w1ins .lt. softlog ) then
         write(io_output,*) 'caught in swap'
         return
      end if

!     --- select one position at random ---
      if ( ichoi .gt. 1 ) then
         rbf = w1ins * random()
         bsum = 0.0d0 
         do ip = 1, ichoi
            if ( .not. lovr(ip) ) then
               bsum = bsum + bfac(ip)
               if ( rbf .lt. bsum ) then
!     --- select ip position ---
                  iwalk = ip
                  goto 180
               end if
            end if
         end do
         write(io_output,*) 'w1ins:',w1ins,'rbf:',rbf
         call cleanup('big time screwup -- w1ins')
      else
         iwalk = 1
      end if
      
 180  v1ins =  vtry(iwalk)  
      v1insext = vtrext(iwalk)
      v1insint = vtrinter(iwalk)
      v1inselc = vtrelect(iwalk)
      v1insewd = vtrewald(iwalk)


!      write(io_output,*)'vtry(iwalk)vtrext(iwalk)vtrinter(iwalk)vtrelect(iwalk)'
!     & ,vtry(iwalk),vtrext(iwalk),vtrinter(iwalk),vtrelect(iwalk),
!     & vtrewald(iwalk)

      rxnew(beg) = rxp(1,iwalk)
      rynew(beg) = ryp(1,iwalk)
      rznew(beg) = rzp(1,iwalk)

      if (lrigid(imolty)) then
!     --- calculate new vector from initial bead

         do j = beg,iunit
            rxnew(j) = rxnew(beg)  - (rxu(irem,beg) - rxu(irem,j))
            rynew(j) = rynew(beg)  - (ryu(irem,beg) - ryu(irem,j))
            rznew(j) = rznew(beg)  - (rzu(irem,beg) - rzu(irem,j))
            qqu(iins,j) = qqu(irem,j)
         end do
         call schedule(igrow,imolty,ifrom,iutry,0,4)
      else if (lring(imolty)) then
         lfixnow = .true.
         call safeschedule(igrow,imolty,ifrom,iutry,findex,2)
      else
         call schedule(igrow,imolty,ifrom,iutry,0,2)
      end if
      

!     ------------------------------------------------------------------
      
      waddnew = 1.0d0
      lterm = .false.

      call rosenbluth( .true.,lterm,iins,iins,imolty,ifrom ,boxins,igrow,waddnew,lfixnow,ctorfn,2 )

!     --- termination of cbmc attempt due to walk termination ---
      if ( lterm ) return

      if ( ldual .or. lewald .or. iunit .ne. igrow  .or. ((.not. lchgall) .and. lelect(imolty)) ) then
!     --- Put on hydrogens for explicit AA model for calculation of COM
!     --- and assign all of the grown new and old beads to rxuion
!     --- with rxuion: new = 2
         iii = 2
         do j=1,igrow
            rxuion(j,iii) = rxnew(j)
            ryuion(j,iii) = rynew(j)
            rzuion(j,iii) = rznew(j)
            qquion(j,iii) = qqu(iins,j)
         end do
         moltion(iii) = imolty
         
         ibox=boxins
         nboxi(iins) = ibox
      end if
      
!     --- Begin DC-CBMC Corrections for NEW configuration
!      if (ldual .or. ((.not. lchgall) .and. lelect(imolty))) then 

      if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) .or. (lchgall .and. lewald .and. (.not. ldual))) then

!     calculate the true site-site energy
         istt = 1
         iett = igrow
!         write(io_output,*) igrow
         
         call energy (iins,imolty, v, vintra,vinter,vext ,velect,vewald,iii,ibox, istt, iett, .true.,ovrlap ,.false.,vdum,.false.,lfavor)
         
         if (ovrlap) then
            write(io_output,*) 'iins',iins,'irem',irem
            call cleanup('strange screwup in DC-CBMC swap')
         end if
! v1insewd, vnewewald and vnewintra now accounted for in v from energy
!$$$         delen = v - ( vnewinter + vnewext +vnewelect) 
!$$$     &        - (v1ins - v1insewd)
         delen = v - ( vnewinter + vnewext +vnewelect + vnewintra + vnewewald + v1ins) 
         waddnew = waddnew*exp(-beta*delen)
         vnewt     = vnewt + delen
         vnewinter = vinter - v1insint
         vnewext   = vext - v1insext
         vnewelect = velect - v1inselc
         vnewewald = vewald - v1insewd
         vnewintra = vintra
      end if
!     End DC-CBMC Corrections for NEW configuration

!     Begin Ewald-sum Corrections
      if ( lewald ) then
! --- reciprocal space sum
! --- prepare qquion(jj,1) etc
         moltion(1) = imolty
         if ( lswapinter ) then
            do j = 1,iunit
               qquion(j,1) = 0.0d0
            end do
            call recip(boxins,vrecipn,vrecipo,1)
            delen = vrecipn - vrecipo 
            waddnew = waddnew*exp(-beta*delen)
            vnewelect = vnewelect + delen
            vnewt = vnewt + delen
         end if
      end if
      
!     End Ewald-sum Corrections

      if (lpbcz) then
          if (lsolid(boxins) .and. .not. lrect(boxins)) then
             ibox = boxins
             volins = cell_vol(ibox)
         else
            volins=boxlx(boxins)*boxly(boxins)*boxlz(boxins)
         end if
      else
         volins=boxlx(boxins)*boxly(boxins)
      end if
      
!     Begin Tail corrections for BOXINS with inserted particle

      if (ltailc .and. lswapinter) then
         vinsta = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty


               if ( jmt .eq. imolty ) then
                  rho = dble( ncmt(boxins,jmt) + 1 ) / volins
               else
                  rho = dble( ncmt(boxins,jmt) ) / volins
               end if
               if ( imt .eq. imolty ) then
                  vinsta = vinsta +  dble( ncmt(boxins,imt) + 1 ) * coru(imt,jmt,rho, boxins)
               else
                  vinsta = vinsta +  dble( ncmt(boxins,imt) ) * coru(imt,jmt,rho, boxins)
               end if
            end do
         end do

         vinsta = vinsta - vtailb( boxins )
         waddnew = waddnew*exp(-beta*vinsta)
         vnewt = vnewt + vinsta
         vnewinter = vnewinter + vinsta
      else
         vinsta = 0.0d0
      end if

!     End Tail corrections for BOXINS with inserted particle

      if ( .not. lanes ) then
         if ( lswapinter ) then
            arg = w1ins * waddnew * weight * volins  / dble( ncmt(boxins,imolty)+1 )
            acchem(boxins,imolty) = acchem(boxins,imolty)+arg
         end if
      end if

      bsswap(imolty,ipairb,boxins+nbox) =  bsswap(imolty,ipairb,boxins+nbox) + 1.0d0

!     Compute weights for the molecule to be removed from boxrem

! *** check that there is at least one molecule in BOXREM ***
      if ( lempty ) then
!         write(io_output,*) 'no molecule in BOXREM'
         if (lgrand) then
	    if (boxrem.eq.2) write(io_output,*) ' ERROR ***** array too low !'
         end if
         return
      end if


! *** select a position of the first/starting unit at RANDOM ***
! *** and calculate the boltzmann weight                     ***
! *** for the chain to be REMOVED                           ***

      rxp(1,1) = rxu(irem,beg)
      ryp(1,1) = ryu(irem,beg)
      rzp(1,1) = rzu(irem,beg)

      ichoi = nchoi1(imolty)

      if (lsolid(boxrem) .and. .not. lrect(boxrem)) then
            ibox = boxrem
            do icbu = 2,ichoi
               sx = random()
               sy = random()
               sz = random()
               rxp(1,icbu) = sx*hmat(ibox,1)+sy*hmat(ibox,4) +sz*hmat(ibox,7)
               ryp(1,icbu) = sx*hmat(ibox,2)+sy*hmat(ibox,5) +sz*hmat(ibox,8)
               rzp(1,icbu) = sx*hmat(ibox,3)+sy*hmat(ibox,6) +sz*hmat(ibox,9)
            end do
      else
            do icbu = 2,ichoi
               rxp(1,icbu) = boxlx(boxrem) * random()
               ryp(1,icbu) = boxly(boxrem) * random()
               if (lpbcz .or. lslit) then
                  rzp(1,icbu) = boxlz(boxrem) * random()
               else if ( lsami .or. lmuir .or. ljoe ) then
                  if ( lempty ) then
                     rzp(1,icbu) = 20*random()-10
                  else
                     rzp(1,icbu) = rzu(irem,1)
                  end if
               else
                  rzp(1,icbu) = 0.0d0
               end if
            end do
      end if
         
! *** calculate the boltzmann weight of first bead          ***

      lnew = .false.
      call boltz(lnew,.true.,ovrlap,irem,irem,imolty,boxrem,ichoi,idum ,1,glist,0.0d0)

      if ( ovrlap ) then
         write(io_output,*) 'disaster: overlap for 1st bead in SWAP'
      end if
! *** calculate the correct weight for the  old  walk ***

      w1rem = 0.0d0
      do ip = 1, ichoi
         w1rem = w1rem + bfac(ip)
      end do

! --- check for termination of walk ---
      if ( w1rem .lt. softlog ) then 
         write(io_output,*) ' run problem : soft overlap in old'
      end if

      v1rem = vtry(1)
      v1remint = vtrinter(1)
      v1remext = vtrext(1)
      v1remelc = vtrelect(1)
      v1remewd = vtrewald(1)

      waddold = 1.0d0

!     --- call rosenbluth for old conformation

      call rosenbluth(.false.,lterm,irem,irem,imolty,ifrom ,boxrem,igrow,waddold,lfixnow,ctorfo,2 )

      if ( lterm ) then 
         write(io_output,*) 'SWAP: rosenbluth old rejected'
         return
      end if

      if ( ldual .or. lewald .or. igrow .ne. iunit  .or. ((.not. lchgall) .and. lelect(imolty)) ) then
!     --- store the old grown beads and explict placed beads positions
!     --- 1 = old conformation
         iii = 1
         do j = 1,iunit
            rxuion(j,1) = rxu(irem,j)
            ryuion(j,1) = ryu(irem,j)
            rzuion(j,1) = rzu(irem,j)
            qquion(j,1) = qqu(irem,j)
         end do
      end if

!     Begin Correction for DC-CBMC for OLD configuration
!      if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) ) then 


      if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) .or. (lchgall .and. lewald .and. (.not. ldual))) then


!     --- correct the acceptance rules 
!     --- calculate the Full rcut site-site energy
         istt=1
         iett = igrow

         call energy (irem,imolty, v, vintra,vinter,vext,velect ,vewald,iii, boxrem, istt, iett, .true.,ovrlap ,.false.,vtorold,.false.,lfavor)
            
         if (ovrlap) call cleanup('disaster ovrlap in old conf SWAP')
! v now includes vnewintra,v1remewd and voldewald, take out
!$$$         deleo = v - ( voldinter + voldext +voldelect) 
!$$$     &        - (v1rem - v1remewd)
         deleo = v - ( voldinter + voldext +voldelect + voldintra  + voldewald + v1rem) 
         waddold = waddold*exp(-beta*deleo)
         voldt     = voldt + deleo
         voldinter = vinter - v1remint
         voldext   = vext - v1remext
         voldelect = velect - v1remelc
         voldewald = vewald - v1remewd
         voldintra = vintra
      end if
!     End Correction for DC-CBMC for OLD configuration

!     Begin Ewald-sum Corrections for OLD configuration
      if ( lewald ) then
! --- reciprocal space sum on r*uion
! --- prepare qquion(jj,1) etc
         if ( lswapinter ) then
            do j = 1,iunit
               qquion(j,2) = 0.0d0
            end do
         end if
         call recip(boxrem,vrecipn,vrecipo,1)
         deleo = vrecipo - vrecipn
         voldt = voldt + deleo 
         voldelect = voldelect + deleo
         waddold = waddold * exp(-beta*deleo)
      end if

!     End Ewald-sum Corrections for OLD configuration
      
      if (lpbcz) then
         if (lsolid(boxrem) .and. .not. lrect(boxrem)) then
            ibox = boxrem
            volrem = cell_vol(ibox)
         else
            volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)
         end if
      else
         volrem=boxlx(boxrem)*boxly(boxrem)
      end if
      
!     Start of intermolecular tail correction for boxrem

      if (ltailc .and. lswapinter) then
!     --- BOXREM without removed particle
         vremta = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty
                  
               if ( jmt .eq. imolty ) then
                  rho = dble( ncmt(boxrem,jmt) - 1 ) / volrem
               else
                  rho = dble( ncmt(boxrem,jmt) ) / volrem
               end if
               if ( imt .eq. imolty ) then
                  vremta = vremta +  dble( ncmt(boxrem,imt) - 1 ) * coru(imt,jmt,rho, boxrem)
               else
                  vremta = vremta +  dble( ncmt(boxrem,imt) ) * coru(imt,jmt,rho, boxrem)
               end if
            end do
         end do
   
         vremta = - vremta + vtailb( boxrem )
         waddold=waddold*exp(-beta*vremta) 
         voldt = voldt + vremta
         voldinter = voldinter + vremta
      else
         vremta = 0.0d0
      end if
!     End of intermolecular tail correction for boxrem

!     --- Add contributions of the first bead and additional beads:

      vnewt     = vnewt  + v1ins
      vnewinter = vnewinter + v1insint
      vnewext   = vnewext + v1insext
      vnewelect = vnewelect + v1inselc
      vnewewald = vnewewald + v1insewd
      
      voldt     = voldt  + v1rem
      voldinter = voldinter+(v1remint)
      voldext   = voldext+v1remext
      voldelect = voldelect + v1remelc
      voldewald = voldewald + v1remewd

      weight= w1ins * waddnew * weight
      weiold= w1rem * waddold * weiold

      wnlog = log10 ( weight )
      wolog = log10 ( weiold )
      wdlog = wnlog - wolog
      
      if ( wdlog .lt. -softcut ) then
!         write(io_output,*) '### underflow in wratio calculation ###'
         return
      end if

      if ( lswapinter ) then
         if (lgibbs) then
!     --- Note: acceptance based on only molecules of type imolty
            wratio = ( weight / weiold ) * ( volins * dble( ncmt(boxrem,imolty) ) /  ( volrem * dble( ncmt(boxins,imolty) + 1 ) ) ) * exp(beta*(eta2(boxrem,imolty)- eta2(boxins,imolty)))

         else if (lgrand) then
            if (boxins.eq.1) then
!              --- molecule added to box 1
               wratio = (weight /  weiold ) *  volins * B(imolty) / (ncmt(boxins,imolty)+1) 
            else
!              --- molecule removed from box 1
               wratio = (weight /weiold)* dble(ncmt(boxrem,imolty))/(volrem*B(imolty))

            end if
         end if
      else
         wratio = (weight*wbias_ins)/(weiold*wbias_rem)
      end if

       
!         wratio = 1.0   

      if ( random() .le. wratio ) then
!         write(io_output,*) 'SWAP MOVE ACCEPTED',irem
! *** we can now accept !!!!! ***
         bnswap(imolty,ipairb,boxins+nbox) =  bnswap(imolty,ipairb,boxins+nbox) + 1.0d0
         if ( .not. lswapinter .and. lbias(imolty) ) then
            if ( lrem_out .and. lins_in ) then
               bnswap_in(imolty,2) = bnswap_in(imolty,2) + 1.0d0
            else if ( (.not. lrem_out) .and. (.not. lins_in ) ) then
               bnswap_out(imolty,2) = bnswap_out(imolty,2) + 1.0d0
            end if
         end if

!---update the position, it will be used to get the bonded energy
         do ic = 1,igrow
            rxu(irem,ic) = rxnew(ic)
            ryu(irem,ic) = rynew(ic)
            rzu(irem,ic) = rznew(ic)
         end do
         do ic = igrow+1,iunit
            rxu(irem,ic) = rxuion(ic,2)
            ryu(irem,ic) = ryuion(ic,2)
            rzu(irem,ic) = rzuion(ic,2)
         end do

!         call Intra_energy(irem,imolty, vdum ,vintran, vdum,vdum,
!     &     velectn,vewaldn,flagon, boxins, 1, iunit,.true.,ovrlap,
!     &      .false.
!     &     ,vdum,.false.,.false.,vvibn,vbendn,vtgn) 


!         total_NBE = vintran+velectn+vewaldn+vtgn+vbendn+vvibn 
          total_NBE = 0.0d0
          vtgn      = 0.0d0
          vbendn    = 0.0d0
          vvibn     = 0.0d0
!         write(io_output,*) vintran,velectn,vewaldn       

!         write(io_output,*) 'irem', irem  

! ---    update energies:

         vbox(boxrem)     = vbox(boxrem)     - voldt - total_NBE
	 vinterb(boxrem)  = vinterb(boxrem)  - voldinter
         vtailb(boxrem)   = vtailb(boxrem)   - vremta
	 vintrab(boxrem)  = vintrab(boxrem)  - voldintra
	 vvibb(boxrem)    = vvibb(boxrem)    - voldbvib - vvibn
	 vtgb(boxrem)     = vtgb(boxrem)     - voldtg - vtgn
	 vextb(boxrem)    = vextb(boxrem)    - voldext
	 vbendb(boxrem)   = vbendb(boxrem)   - voldbb - vbendn
         velectb(boxrem)  = velectb(boxrem)  - (voldelect+voldewald)
         vflucqb(boxrem)  = vflucqb(boxrem)  - voldflucq
	
 
         vbox(boxins)     = vbox(boxins)     + vnewt + total_NBE
	 vinterb(boxins)  = vinterb(boxins)  + vnewinter
	 vtailb(boxins)   = vtailb(boxins)   + vinsta
	 vintrab(boxins)  = vintrab(boxins)  + vnewintra
	 vvibb(boxins)    =  vvibb(boxins)   + vnewbvib + vvibn
	 vtgb(boxins)     = vtgb(boxins)     + vnewtg + vtgn
	 vextb(boxins)    = vextb(boxins)    + vnewext
	 vbendb(boxins)   = vbendb(boxins)   + vnewbb + vbendn
         velectb(boxins)  = velectb(boxins)  + (vnewelect+vnewewald)
         vflucqb(boxins)  = vflucqb(boxins)  + vnewflucq

! ---    update book keeping

         if ( lswapinter ) then
            nboxi(irem) = boxins
         
            parbox(ncmt(boxins,imolty)+1,boxins,imolty)= irem
            parbox(pointp,boxrem,imolty)= parbox(ncmt(boxrem,imolty),boxrem,imolty)
            parbox(ncmt(boxrem,imolty),boxrem,imolty)=0
            
            nchbox(boxins) = nchbox(boxins) + 1
            nchbox(boxrem) = nchbox(boxrem) - 1
            ncmt(boxins,imolty) = ncmt(boxins,imolty) + 1
            ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1

            if ( lexpand(imolty) ) then
               itype = eetype(imolty)
               ncmt2(boxins,imolty,itype) =  ncmt2(boxins,imolty,itype) + 1
               ncmt2(boxrem,imolty,itype) =  ncmt2(boxrem,imolty,itype) + 1
            end if
         end if

!         do ic = 1,igrow
!            rxu(irem,ic) = rxnew(ic)
!            ryu(irem,ic) = rynew(ic)
!            rzu(irem,ic) = rznew(ic)
!         end do
!         do ic = igrow+1,iunit
!            rxu(irem,ic) = rxuion(ic,2)
!            ryu(irem,ic) = ryuion(ic,2)
!            rzu(irem,ic) = rzuion(ic,2)
!         end do
         
         if ( lewald ) then
! *** update reciprocal-space sum
            if ( lswapinter ) then
               call recip(boxins,vdum,vdum,2)
               call recip(boxrem,vdum,vdum,2)
            else
               call recip(boxins,vdum,vdum,2)
            end if
         end if

! ---    update center of mass
         call ctrmas(.false.,boxins,irem,3)
! *** update linkcell, if applicable
         if ( licell .and. ((boxins .eq. boxlink) .or. (boxrem .eq.   boxlink))) then
            call linkcell(2,irem,vdum,vdum,vdum,ddum)
         end if
         
!         write(io_output,*) 'lneighbor:',lneighbor

         if ( lneigh ) call updnn( irem )
         if ( lneighbor ) then
            neigh_old = neigh_cnt(irem)
            if ( neigh_old .le. 6 .and. neigh_icnt .le. 6 ) then
               cnt_wf2(neigh_old,neigh_icnt,ip) = cnt_wf2(neigh_old,neigh_icnt,ip)+1
            end if

            do 10 ic = 1, neigh_old
               j = neighbor(ic,irem)
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. irem ) then
                     neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                     neigh_cnt(j) = neigh_cnt(j)-1
                     goto 10
                  end if
               end do
 10         continue
            neigh_cnt(irem) = neigh_icnt
!            write(io_output,*) 'irem:',irem,neigh_icnt
            do ic = 1,neigh_icnt
               j = neighi(ic)
               neighbor(ic,irem)=j
               lneighij = .false.
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. irem ) then
                     lneighij = .true.
                  end if
               end do
               if ( .not. lneighij ) then
                  neigh_cnt(j) = neigh_cnt(j)+1
                  neighbor(neigh_cnt(j),j) = irem
               end if
            end do
         end if

!         write(io_output,*) irem,'end SWAP'

      end if
 
! -----------------------------------------------------------------
      return
      end

