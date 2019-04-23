      subroutine analysis(switch)

!     *******************************************************************
!     *** Modified to perform analysis on the fly based on anal10.f   ***
!     *** [Marcus Martin] by Neeraj Rai 07/14/04                      ***
!     *******************************************************************

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
!$$$      include 'inputdata.inc'
!$$$      include 'connect.inc'
!$$$      include 'system.inc'
!$$$      include 'cell.inc' 
!$$$      include 'gor.inc'
!$$$      include 'peboco.inc'

      character(LEN=default_path_length)::fileout
      character(LEN=default_string_length)::ftemp,fname2,fname3
      integer(KIND=normal_int)::switch,ichain,fname
      logical::lskip
      integer(KIND=normal_int)::bin,k,kk,box,z
      integer(KIND=normal_int)::tempbx,dummy,xx,yy,g,gg
      integer(KIND=normal_int)::chnum,imolty,i,ii,j,jj ,jmolty,binadj,ntij,ntii,ntjj,jstart,istart,ntji
      integer(KIND=normal_int)::iivib,ip1,ip2,ip3,gaudef,uu,dum,units,zzz,zz1

      real(KIND=double_precision)::avolume,rho,binstep,vec_hist ,rxui,ryui,rzui,xxideal ,rxuij,ryuij,rzuij,ruijsq,ruij ,numxx,numyy,count,const,rlower,rupper,nideal,analhist
      real(KIND=double_precision)::comanalhist
      real(KIND=double_precision)::shlsumx,shlsumy,rcutsq

!   NOW DEFINED IN CONTROL.INC

!c     ntmax = number of moltyps, ntdifmx = max number of diff beads in sim
!      parameter (nbinmx=200,nbxmax=1,ntmax=5,nmax=1600,numax=18
!     & ,ntdifmx=8)

      integer(KIND=normal_int)::bend,iv,iuvib,iuv,iutest
      real(KIND=double_precision)::value,total,degree
      integer(KIND=normal_int)::torsion, tor_code,itor,iutor,patt,bthree ,decimal,power
      real(KIND=double_precision)::xcc,ycc,zcc,tcc,fplus,fminus ,ftrans
      real(KIND=double_precision)::xvec,yvec,zvec,distij,xaa1,yaa1,zaa1,xa1a2 ,ya1a2,za1a2,daa1,da1a2,dot,thetac,theta ,ratio ,psum,tempzcm,tempmasst,slab_vol

      dimension patt(tor_max),tempzcm(nbxmax),tempmasst(nbxmax)
      dimension tor_code(numax,numax)
      dimension vec_hist(nbxmax,ntmax,nbinmax_ete)  

!******** Taking out charge part*************************

!Cc     --- variables used in the charge parts
!C     logical::qhere
!C      integer(KIND=normal_int)::qbin,qbinmax,qbins,qqcode,imol,iunit
!C      integer(KIND=normal_int)::qbin,qbins,qqcode,imol,iunit
!C    NOW DEFINED IN CONTROL.INC
!C      parameter (qbinmax=1000)

!C      real(KIND=double_precision)::qdisp
!C     real(KIND=double_precision)::qmin,qmax,qdiff,qstep,qdummy
!C      dimension qanalhist(numax*ntdifmx+numax,nbxmax,qbinmax)
!C      dimension qcount(numax*ntdifmx+numax,nbxmax)

!Cc     --- variables used in the dipole parts
!C     logical::ldipole
!C     integer(KIND=normal_int)::qblock,block,nblock
!C      real(KIND=double_precision)::dipole,dicount,dipx,dipy,dipz,diconv,stddev
!C     &     ,avera,diprev,dcprev,dipblk
!C     dimension dipole(ntdifmx,nbxmax),dicount(ntdifmx,nbxmax)
!C     dimension diprev(ntdifmx,nbxmax),dcprev(ntdifmx,nbxmax)
!C      dimension dipblk(ntdifmx,nbxmax,20)
!********* end charge part ********************************

      dimension analhist(ntdifmx*ntdifmx*ntmax*ntmax,nbinmx)
      dimension comanalhist(ntmax*ntmax,nbinmx)
      dimension avolume(nbxmax)
      dimension rho(ntmax*ntmax*ntdifmx)
      dimension count(nbxmax,ntdifmx*ntmax)
      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax) ,distij(numax,numax)


      if(switch.eq.0) then

      if(nhere.gt.ntdifmx) then
        write(io_output,*) 'nhere greater than ntdifmax', nhere, ntdifmx
        call cleanup('choose a larger ntdifmx in control.inc')
      end if
 
      if(nbin.gt.nbinmx) then
         write(io_output,*) 'number of bins "nbin" .gt. nbinmx', nbin, nbinmx
         call cleanup('choose a larger nbinmx in control.inc')
       end if  
      
      do kk=1,nbxmax
         max_boxlz(kk)=boxlz(kk)
      end do
      

      if (lrdf) then
         if ( lintra ) then
            nskip = 0.5d0
         else
            nskip = 1.5d0
         end if

      end if

      lstretch = .false.

!************************************************************************
!  THIS PART SHOULD GO TO READDAT

!      if (lcharge) then
!         write(io_output,*) 'input minimum charge for dist (qmin)'
!         read(5,*) qmin
!         write(io_output,*) 'input maximum charge for dist (qmax)'
!         read(5,*) qmax
!         write(io_output,*) 'input qbins must be less than',qbinmax
!         read(5,*) qbins
!         if ( qmax .lt. qmin ) call cleanup('qmax cannont be less than qmin')
!         qdiff = qmax - qmin
!         qstep = qdiff/dble(qbins)
!         write(io_output,*) 'input additional disp for charges (0.0 for none)'
!         read(5,*) qdisp
!      end if
!      write(io_output,*) 'do you want the average dipole moments?'
!      read(5,*) ldipole
!      if ( ldipole ) then
!         write(io_output,*) 'number of blocks for dipole?'
!         read(5,*) qblock
!      end if
!************************************************************************



!  end to end vector probability distribution

      if(lete) then
         do imolty = 1,nmolty
            max_length(imolty)=0.0d0
            do i = 1,nunit(imolty)
               do j = 1,invib(imolty,i)
                  if(i<ijvib(imolty,i,j)) then
                     max_length(imolty)=max_length(imolty) +brvib(itvib(imolty,i,j))
                  end if 
               end do
            end do
         end do  
      
   
         do imolty=1,nmolty
            i = (aint((max_length(imolty))/bin_width)+1) 
            if(i.gt.nbinmax_ete) then
	       write(io_output,*) 'Stopped in Analysis end-to-end dist calc'
               write(io_output,*) 'number of bins greater than nbinmax_ete',i, nbinmax_ete
               call cleanup('choose larger nbinmax_ete in control.inc')
            end if 
         end do
         
!      initializing the arrays

         do kk= 1,nbxmax
            do i=1,nmolty
               do j= 1,nbinmax_ete
                  end_to_end(kk,i,j)=0.0d0
               end do
            end do   
         end do
         
      end if 

!      print*, '*************'
!      print*, '*************'
!      print*,  max_length(1)
!      print*,  max_length(2)
!      print*,  max_length(3)
!      print*, '*************'
!      print*, '*************'


      if(lrhoz) then
         
         
         do kk=1,nbox
            i=(aint(boxlz(kk)/bin_width)+1)
            if(i.gt.nbinmax_ete) then
               write(io_output,*) 'Stopped in Analysis Z profile calculation'  
               write(io_output,*) 'number of bins greater than nbinmax_ete',i, nbinmax_ete
               call cleanup('choose larger nbinmax_ete in control.inc')
            end if
         end do 
         
         do kk=1,nbox
            do i = 1,nmolty 
               do j = 1, nbinmax_ete
                  bigrhoz(kk,i,j)=0.0d0
               end do
               do j=-nbinmax_ete,nbinmax_ete
                  bigboxcom_rhoz(kk,i,j)=0.0d0
               end do
            end do
         end do
      end if 
      



      if ( lbend ) then
!     --- compute ang_bin_size
         ang_bin_size = onepi/dble(ang_bin_max)
!     --- figure out how many angles exist and assign them numbers
         do imolty = 1,nmolty
            bend = 0
            do ii = 1,nunit(imolty)
               do iv = 1,invib(imolty,ii)
                  iuvib = ijvib(imolty,ii,iv)
!     --- determine whether the beads connected to this unit
!     --- are of higher index than ii
                  do iuv = 1,invib(imolty,iuvib)
                     iutest = ijvib(imolty,iuvib,iuv)
                     if ( iutest .gt. ii ) then
                        bend = bend + 1
                        angle_1(imolty,bend) = ii
                        angle_2(imolty,bend) = iuvib
                        angle_3(imolty,bend) = iutest
                     end if
                  end do
               end do
            end do
            angle_num(imolty) = bend
            if(bend.gt.angle_max) then
               write(io_output,*) 'number of bends greater than angle_max', 'molecule type', imolty,' bends', bend   
               call cleanup('choose a larger angle_max in control.inc')
            end if
         end do
      end if
      
      
      if ( lgvst ) then
!     --- compute tor_bin_size
         tor_bin_size = 360.0d0/dble(tor_bin_max)
!     --- figure out how many torsions exist and assign them numbers
         do imolty = 1,nmolty
            torsion = 0
            do ii = 1,nunit(imolty)
               do itor = 1,intor(imolty,ii)
                  iutor = ijtor4(imolty,ii,itor)
!     --- determine whether final bead connected to this unit
!     --- is of higher index than ii
                  if ( iutor .gt. ii ) then
                     torsion = torsion + 1
                     tor_1(imolty,torsion) = ii
                     tor_2(imolty,torsion) = ijtor2(imolty,ii,itor)
                     tor_3(imolty,torsion) = ijtor3(imolty,ii,itor)
                     tor_4(imolty,torsion) = iutor
                     tor_code(ii,iutor) = torsion
                  end if
               end do
            end do
            tor_num(imolty) = torsion
            if(torsion.gt.tor_max) then
               write(io_output,*) 'number of torsions greater than tor_max', 'molecule type', imolty,'torsions', torsion
               call cleanup('choose a larger tor_max in control.inc*** requires memory 3**torsion so choose it equal to torsions')
             end if
         end do
      end if


! --- initialize the arrays

      if ( lbend ) then
         do kk = 1,nbox
            do imolty = 1,nmolty
               do bend = 1,angle_num(imolty)
                  do bin = 1,ang_bin_max
                     angle_bin(kk,imolty,bend,bin) = 0.0d0
                     angle_tot(kk,imolty,bend) = 0.0d0
                  end do
               end do
            end do
         end do
      end if

      if ( lgvst ) then
         do kk = 1,nbox
            do imolty = 1,nmolty
               do torsion = 1,tor_num(imolty)
                  do bin = 1,tor_bin_max
                     tor_bin(kk,imolty,torsion,bin) = 0.0d0
                     tor_tot(kk,imolty,torsion) = 0.0d0
                  end do

               end do
            end do
         end do

         do xx = 1,ntmax
            do zzz = 1,nbox
               do uu = 1,tor_max
                  gdefect(zzz,xx,uu) = 0.0d0
                  btrans(zzz,xx,uu) = 0.0d0
                  bg_plus(zzz,xx,uu) = 0.0d0
                  bg_minus(zzz,xx,uu) = 0.0d0
               end do
               gdefect(zzz,xx,tor_max+1) = 0.0d0
            end do
         end do

      end if

      if ( lrdf ) then
         dummy = ntdifmx*ntdifmx*ntmax*ntmax
         do ntij = 1,dummy
            do kk=1,nbox
               nnone(kk,ntij) = 0.0d0
               binadj = (kk-1)*dummy + ntij
               do bin = 1,nbin
                  biganalhist(binadj,bin) = 0.0d0
                  shell(binadj,bin,1) = 0.0d0
                  shell(binadj,bin,2) = 0.0d0
               end do
            end do
         end do

         dummy = ntmax*ntmax
         do ntij = 1,dummy
            do kk=1,nbox
               comnone(kk,ntij) = 0.0d0
               binadj = (kk-1)*dummy + ntij
               do bin = 1,nbin
                  combiganalhist(binadj,bin) = 0.0d0
                  comshell(binadj,bin,1) = 0.0d0
                  comshell(binadj,bin,1) = 0.0d0 
               end do
            end do
         end do
      end if

!      if ( lcharge ) then
!         do kk = 1, nbox
!            do imol = 1,ntdifmx
!               do iunit = 1,numax
!                  qqcode = numax*(imol-1) + iunit
!                  do qbin = 1,qbins
!                     qanalhist(qqcode,kk,qbin) = 0.0d0
!                     qcount(qqcode,kk) = 0.0d0
!                  end do
!               end do
!            end do
!         end do
!      end if
      
!      if ( ldipole ) then
!         do kk = 1,nbox
!            do imol = 1,nmolty
!               dipole(imol,kk) = 0.0d0
!               dicount(imol,kk) = 0.0d0
!               diprev(imol,kk) = 0.0d0
!               dcprev(imol,kk) = 0.0d0
!            end do
!         end do
!      end if
      return
      end if 




!     cycle through all of the frames and calculate the radial dist func

      if(switch.eq.1) then
         
         do kk=1,nbox
             if(boxlz(kk).gt.max_boxlz(kk)) then
                max_boxlz(kk)=boxlz(kk)
             end if
         end do

         nframe = nframe + 1.0d0
         
         do kk = 1,nbox
           if(lsolid(kk).and.(.not.(lrect(kk)))) then
                 avolume(kk) = hmat(kk,1)*hmat(kk,5)*hmat(kk,9)
	   else      
                 avolume(kk) = boxlx(kk)*boxly(kk)*boxlz(kk)
 	  end if 
         end do

!     initialize frame specific arrays
         do ii = 1,nmolty
            do gg = 1,nbox
               cmolec(ii,gg) = 0.0d0
               do g = 1,nhere
                  dummy = nhere*(ii-1) + g
                  count(gg,dummy) = 0.0d0
               end do
            end do
         end do
         
         do ichain = 1,nchain
            imolty= moltyp(ichain)
            tempbx = nboxi(ichain)
            cmolec(imolty,tempbx) = cmolec(imolty,tempbx) + 1.0d0
            do z = 1, nunit(imolty)
               dummy = nhere*(imolty-1)+decode(ntype(imolty,z))
               count(tempbx,dummy) = count(tempbx,dummy) + 1.0d0
            end do
         end do

!         if (lcharge) then
!            do kk = 1,nbox
!               do i = 1,nchain
!                  if ( nboxi(i) .eq. kk ) then
!                     imolty = moltyp(i)
!                     do ii = 1,nunit(imolty)
!                        qdummy = qqu(i,ii)
!                        if ( qdummy .lt. qmin .or. qdummy .gt. qmax) then
!                           write(io_output,*) 'q value out of range'
!                           write(io_output,*) 'frame,i,ii,q',k,i,ii,qdummy
!                           call cleanup('q value out of range')
!                        end if
!                        qdummy = qdummy - qmin
!                        qbin = aint(qdummy/qstep) + 1
!                        qqcode = numax*(imolty-1) + ii
!                        qanalhist(qqcode,kk,qbin) = 
!     &                       qanalhist(qqcode,kk,qbin) + 1.0d0
!                        qcount(qqcode,kk) = qcount(qqcode,kk) + 1.0d0
!                     end do
!                  end if
!               end do
!            end do
!         end if

!         if ( ldipole ) then
!            do i = 1,nchain
!               kk = nboxi(i)
!               imolty = moltyp(i)
!               dipx = 0.0d0
!               dipy = 0.0d0
!               dipz = 0.0d0
!               do iunit = 1,nunit(imolty)
!                  qdummy = qqu(i,iunit)
!                  dipx = dipx + qdummy*rxu(i,iunit)
!                  dipy = dipy + qdummy*ryu(i,iunit)
!                  dipz = dipz + qdummy*rzu(i,iunit)
!               end do
!               dipole(imolty,kk) = dipole(imolty,kk) 
!     &              + sqrt(dipx*dipx + dipy*dipy + dipz*dipz)
!               dicount(imolty,kk) = dicount(imolty,kk) + 1.0d0
!            end do
!
!            if ( mod(k,block) .eq. 0 ) then
!               nblock = nblock + 1
!               do kk = 1,nbox
!                  do imol = 1,nmolty
!                     if ( dicount(imolty,kk) .gt. 0.5d0 ) then
!                        dipblk(imol,kk,nblock) = 
!     &                       (dipole(imolty,kk) - diprev(imolty,kk))/
!     &                       (dicount(imolty,kk) - dcprev(imolty,kk))
!                     end if
!                     diprev(imolty,kk) = dipole(imolty,kk)
!                     dcprev(imolty,kk) = dicount(imolty,kk)
!                  end do
!               end do
!            end if
!
!         end if
        
!   Compute end to end vector distribution


        if(lete) then
      
           do kk=1,nbox
              do imolty=1,nmolty
                  do i=1,nbinmax_ete
                      vec_hist(kk,imolty,i) =0.0d0
                  end do
               end do
            end do      

           do i=1,nchain
              kk = nboxi(i)
              imolty = moltyp(i)
 
              rxuij = rxu(i,1)- rxu(i,nunit(imolty))
              ryuij = ryu(i,1)- ryu(i,nunit(imolty))
              rzuij = rzu(i,1)- rzu(i,nunit(imolty))

              ruijsq = rxuij*rxuij + ryuij*ryuij+rzuij*rzuij
              
              ruij   = sqrt(ruijsq)
            
              vec_hist(kk,imolty,(int(ruij/bin_width)+1))= vec_hist(kk,imolty,(int(ruij/bin_width)+1))+1

           end do
 
           do kk=1,nbox
             do imolty=1,nmolty
                 xx =  (int(max_length(imolty)/bin_width)+1)
                do bin=1,xx 
                  vec_hist(kk,imolty,bin)=vec_hist(kk,imolty,bin)/ ncmt(kk,imolty)
                
                  end_to_end(kk,imolty,bin)=end_to_end(kk,imolty,bin)+ vec_hist(kk,imolty,bin)
                end do    
             end do
           end do 

        end if 

        if(lrhoz) then

         do kk=1,nbox
           do i = 1,nmolty
             do j = 1, nbinmax_ete
                 rhoz(kk,i,j)=0.0d0
             end do
             do j=-nbinmax_ete,nbinmax_ete
                 boxcom_rhoz(kk,i,j)=0.0d0
             end do
           end do
         end do
           
!! CALCULATE CENTER OF MASS OF EACH BOX

         do kk=1,nbox
            tempzcm(kk)=0.0d0 
            tempmasst(kk) = 0.0d0
            do i =1,nchain
               if(nboxi(i).eq.kk) then
                   tempzcm(kk)   = tempzcm(kk)+mass(moltyp(i))*zcm(i) 
                   tempmasst(kk) = tempmasst(kk)+mass(moltyp(i))
               end if
            end do
            tempzcm(kk)=tempzcm(kk)/tempmasst(kk)  
         end do
         
         do i=1,nchain
            kk=nboxi(i)
            imolty=moltyp(i) 
            rhoz(kk,imolty,int(zcm(i)/bin_width)+1)=rhoz(kk,imolty, int(zcm(i)/bin_width)+1) +1

            rzuij=zcm(i)-tempzcm(kk)
            boxcom_rhoz(kk,imolty,int(rzuij/bin_width))=boxcom_rhoz(  kk,imolty,int(rzuij/bin_width))+1          
         end do  

         do kk=1,nbox
            do imolty=1,nmolty
               xx = (int(max_boxlz(kk)/bin_width))+1
               do bin=1,xx
                 slab_vol = hmat(kk,1)*hmat(kk,5)*bin_width
                 rhoz(kk,imolty,bin)=rhoz(kk,imolty,bin)/slab_vol
                 bigrhoz(kk,imolty,bin)=bigrhoz(kk,imolty,bin)+ rhoz(kk,imolty,bin)
               end do
              
               do bin=-xx,xx
                 boxcom_rhoz(kk,imolty,bin)=boxcom_rhoz(kk,imolty,bin)/ slab_vol
                 bigboxcom_rhoz(kk,imolty,bin)=bigboxcom_rhoz(kk,imolty, bin)+boxcom_rhoz(kk,imolty,bin)
               end do 
            end do
         end do 
        end if

 
       if (lrdf) then
!        --- compute the radial distribution analhistograms for this frame
         do kk = 1,nbox
            rcutsq = rcut(kk)*rcut(kk)
            binstep = rcut(kk)/dble(nbin)
!           initiallize analhist
            do g = 1,nmolty*nmolty*nhere*nhere
               do gg = 1,nbin
                  analhist(g,gg) = 0.0d0
               end do
            end do
!           initiallize comanalhist
            do g = 1,nmolty*nmolty
               do gg = 1,nbin
                  comanalhist(g,gg) = 0.0d0
               end do
            end do

            do i = 1, nchain
 
!  check if i is in relevant box 
               if ( nboxi(i) .eq. kk ) then
                  imolty = moltyp(i)
!     --- loop over all beads ii of chain i 
                  do ii = 1, nunit(imolty)
                     ntii = nhere*(imolty-1) + decode(ntype(imolty,ii))
                     rxui = rxu(i,ii)
                     ryui = ryu(i,ii)
                     rzui = rzu(i,ii)

!     --- loop over all chains j with j>=i 
                     if ( lintra ) then
                        istart = i
                     else
                        istart = i+1
                     end if
                     do j = istart, nchain
!                       ### check for simulation box ###
                        if ( nboxi(j) .eq. kk ) then
                           jmolty = moltyp(j)
!                       --- loop over all beads jj of chain j 
                           if (i .eq. j) then 
                              jstart = ii+1
                           else
                              jstart = 1
                           end if
                           do jj = jstart, nunit(jmolty)
                        
                              ntjj = nhere*(jmolty-1) +  decode( ntype(jmolty,jj) )
                              if ( ntii .gt. ntjj ) then
                                 ntij = (ntjj-1)*nhere*nmolty + ntii
                              else
                                 ntij = (ntii-1)*nhere*nmolty + ntjj
                              end if

                              rxuij = rxui - rxu(j,jj)
                              ryuij = ryui - ryu(j,jj)
                              rzuij = rzui - rzu(j,jj)

                              if ( lpbc )  call mimage (rxuij,ryuij,rzuij,kk)

! *** minimum image the pair separations ***
!                              if ( rxuij .gt. hbx ) then
!                                 rxuij=rxuij-bx
!                              else
!                                 if (rxuij.lt.-hbx) rxuij=rxuij+bx
!                              end if
!
!                              if ( ryuij .gt. hby ) then
!                                 ryuij=ryuij-by
!                              else
!                                 if (ryuij.lt.-hby) ryuij=ryuij+by
!                              end if

!                              if (rzuij.gt.hbz) then
!                                 rzuij=rzuij-bz
!                              else
!                                 if (rzuij.lt.-hbz) rzuij=rzuij+bz
!                              end if
!
                              ruijsq = rxuij*rxuij + ryuij*ryuij  + rzuij*rzuij

                              if (ruijsq .lt. rcutsq) then
                                 ruij = sqrt(ruijsq)

                                 bin = aint(ruij/binstep) + 1 
                                 analhist(ntij,bin) = analhist(ntij, bin) + 1.0d0
                              end if

                           end do
                        end if
                     end do
                  end do
               end if
            end do

!     normalize the analhistogram and add it to the big analhistogram
            
            do ii = 1,nmolty
               do xx  = 1,nhere
                  ntii = nhere*(ii-1) + xx
                  rho(ntii) = (4.0d0 * 3.1415d0 * count(kk,ntii))  / (3.0d0 * avolume(kk) )
               end do
            end do
            
            do ii = 1,nmolty
               do xx = 1,nhere
                  ntii = nhere*(ii-1)+xx
                  numxx = count(kk,ntii)
                  
                  do jj = 1,nmolty
                     do yy = 1,nhere
                        ntjj = nhere*(jj-1)+yy
                        numyy = count(kk,ntjj)
                        ntij = (ntii-1)*nhere*nmolty + ntjj

                        ntji = (ntjj-1)*nhere*nmolty + ntii

                        binadj = (kk-1)*nhere*nhere*nmolty*nmolty  + ntij
                        if ( ntii .eq. ntjj ) then
                           const = rho(ntjj)/2.0d0
                        else
                           const = rho(ntjj)
                        end if
                        
                        shlsumx = 0.0d0
                        shlsumy = 0.0d0  
!     check to see if there are enough molecules to have any interactions
!     if not then lskip is set to true

                        lskip = .false.
                        if ( ntii .eq. ntjj ) then
                           if ( cmolec(ii,kk) .lt. nskip ) then
                              lskip = .true.
                           else if ( numxx .lt. 0.5d0 ) then
                              lskip =.true.
                           end if
                        else
                           if ( cmolec(ii,kk)*cmolec(jj,kk)  .lt. 0.5d0 .or. numxx*numyy  .lt. 0.5d0) then
                              lskip = .true.
                           end if
                        end if
                        
                        if ( .not. lskip  ) then

                           do bin = 1,nbin

                              if ( ntii .eq. ntjj ) then
                                    shlsumx = shlsumx + 2.0d0*analhist(ntij,bin) /numxx
                                    shell(binadj,bin,1) = shell(binadj,bin,1) + shlsumx
                                 else
                                    shlsumx = shlsumx + analhist(ntij,bin)/numxx
                                    shlsumy = shlsumy + analhist(ntij,bin)/numyy
                                    shell(binadj,bin,1) = shell(binadj,bin,1) + shlsumx
                                    shell(binadj,bin,2) = shell(binadj,bin,2) + shlsumy
                                 end if

                              
                              rlower = dble(bin-1)*binstep
                              rupper = rlower + binstep
                              nideal = const*(rupper**3 -rlower**3)
                              xxideal = numxx*nideal

                              biganalhist(binadj,bin) =  biganalhist(binadj,bin)  + analhist(ntij,bin)/xxideal


                           end do
                           
                        else
                           nnone(kk,ntij) = nnone(kk,ntij) + 1.0d0
                        end if
                     end do
                  end do
               end do
            end do

!           --- Center of Mass Radial Distribution Functions
            do i = 1, nchain
 
!  check if i is in relevant box 
               if ( nboxi(i) .eq. kk ) then
                  imolty = moltyp(i)

                  rxui = xcm(i)
                  ryui = ycm(i)
                  rzui = zcm(i)

!     --- loop over all chains j with j>=i 
                  istart = i+1
                  do j = istart, nchain
!                       ### check for simulation box ###
                     if ( nboxi(j) .eq. kk ) then
                        jmolty = moltyp(j)

                        if ( imolty .gt. jmolty ) then
                           ntij = (jmolty-1)*nmolty + imolty
                        else
                           ntij = (imolty-1)*nmolty + jmolty
                        end if

                        rxuij = rxui - xcm(j)
                        ryuij = ryui - ycm(j)
                        rzuij = rzui - zcm(j)

! *** minimum image the pair separations ***
                        if ( lpbc ) call mimage (rxuij,ryuij,rzuij,kk)         

                        ruijsq = rxuij*rxuij + ryuij*ryuij  + rzuij*rzuij
                        
                        if (ruijsq .lt. rcutsq) then
                           ruij = sqrt(ruijsq)
                           
                           bin = aint(ruij/binstep) + 1 
                           comanalhist(ntij,bin) = comanalhist(ntij,bin)  + 1.0d0
                        end if

                     end if
                  end do
               end if
            end do

!     normalize the COM analhistogram and add it to the big analhistogram
            
            do ii = 1,nmolty
               rho(ii) = (4.0d0 * 3.1415d0 * cmolec(ii,kk))  / (3.0d0 * avolume(kk) )
            end do
            
            do ii = 1,nmolty
               ntii = ii
               numxx = cmolec(ntii,kk)
               do jj = 1,nmolty

                  ntjj = jj
                  numyy = cmolec(ntjj,kk)
                  
                  ntij = (ntii-1)*nmolty + ntjj
                  ntji = (ntjj-1)*nmolty + ntii

                  binadj = (kk-1)*nmolty*nmolty + ntij
                  
                  if ( ntii .eq. ntjj ) then
                     const = rho(ntjj)/2.0d0
                  else
                     const = rho(ntjj)
                  end if
                  
                  shlsumx = 0.0d0
                  shlsumy = 0.0d0

!     check to see if there are enough molecules to have any interactions
!     if not then lskip is set to true

                  lskip = .false.
                  
                  if ( ntii .eq. ntjj ) then
                     if ( cmolec(ii,kk) .lt. nskip ) then
                        lskip = .true.
                     end if
                  else if ( cmolec(ii,kk)*cmolec(jj,kk)  .lt. 0.5d0 .or. numxx*numyy  .lt. 0.5d0) then
                     lskip = .true.
                  end if
                           
                  if ( .not. lskip  ) then

                     do bin = 1,nbin

	               if ( ntii .eq. ntjj ) then
                           shlsumx = shlsumx + 2.0d0*comanalhist(ntij,bin) /numxx
                           comshell(binadj,bin,1) = comshell(binadj,bin,1)+ shlsumx
                        else
                           shlsumx = shlsumx + comanalhist(ntij,bin)/numxx
                           shlsumy = shlsumy + comanalhist(ntij,bin)/numyy
                           comshell(binadj,bin,1) = comshell(binadj,bin,1)+ shlsumx
                           comshell(binadj,bin,2) = comshell(binadj,bin,2)+ shlsumy
                        end if


                        rlower = dble(bin-1)*binstep
                        rupper = rlower + binstep
                        nideal = const*(rupper**3 -rlower**3)
                        xxideal = numxx*nideal

                        if (jj .ge. ii) then

                           combiganalhist(binadj,bin) =   combiganalhist(binadj,bin) +  comanalhist(ntij,bin)/xxideal

                        end if

                     end do
                              
                  else
                     comnone(kk,ntij) = comnone(kk,ntij) + 1.0d0
                  end if
               end do
            end do

         end do
        end if

!ccccccccccccccccccccccccc

!     analyse the torsional angles

        if (lstretch .or. lbend .or. lgvst ) then
           do kk=1,nbox
              do i = 1, nchain
          
               imolty = moltyp(i)

!  check if i is in relevant box 
               if ( nboxi(i) .eq. kk ) then

!                 --- calculate all bonds vectors and lengths

                  do ii = 1, nunit(imolty)
                     rxui=rxu(i,ii)
                     ryui=ryu(i,ii)
                     rzui=rzu(i,ii)
                     do iivib = 1, invib(imolty,ii)
                        jj = ijvib(imolty,ii,iivib)
                        xvec(ii,jj) = rxu(i,jj) - rxui
                        yvec(ii,jj) = ryu(i,jj) - ryui
                        zvec(ii,jj) = rzu(i,jj) - rzui
                        distij(ii,jj) = sqrt( xvec(ii,jj)**2 + yvec(ii,jj)**2 + zvec(ii,jj)**2 )
                     end do
                  end do
                  
! - stretching -
                  if ( lstretch ) then
                  end if

! - bending -
                  if ( lbend ) then
                     do bend = 1,angle_num(imolty)
                        j = angle_1(imolty,bend)
                        ip1 = angle_2(imolty,bend)
                        ip2 = angle_3(imolty,bend)

                        thetac = ( xvec(ip1,j)*xvec(ip1,ip2) + yvec(ip1,j)*yvec(ip1,ip2) + zvec(ip1,j)*zvec(ip1,ip2) ) / ( distij(ip1,j)*distij(ip1,ip2) )
                        theta = acos(thetac)
                        if(abs(theta).lt.0.000001) then
                            theta =0.0d0
                        end if
!                        write(io_output,*) 'value of theta',theta
                        angle_avg(kk,imolty,bend) =  angle_avg(kk,imolty,bend) + theta
                        bin = aint(dble(theta)/dble(ang_bin_size))+1
	                if (bin.lt.0) then
                            write(io_output,*) 'theta, ang_bin_size, bin',theta, ang_bin_size, bin
                            bin=1
		        end if
                        angle_bin(kk,imolty,bend,bin) =  angle_bin(kk,imolty,bend,bin) + 1.0d0
                        angle_tot(kk,imolty,bend) =  angle_tot(kk,imolty,bend)+ 1.0d0
                     end do
                  end if

! - torsions -
                  if ( lgvst ) then
!  molecule with dihedral potenials 
                     gaudef = 1
                     dum = 0
                     do torsion = 1, tor_num(imolty)
                        j = tor_1(imolty,torsion)
                        ip1 = tor_2(imolty,torsion)
                        ip2 = tor_3(imolty,torsion)
                        ip3 = tor_4(imolty,torsion)

!              ***  calculate cross product d_a x d_a-1  ***
                        xaa1 = yvec(ip1,j) * zvec(ip2,ip1) + zvec(ip1,j) * yvec(ip1,ip2)
                        yaa1 = zvec(ip1,j) * xvec(ip2,ip1) + xvec(ip1,j) * zvec(ip1,ip2)
                        zaa1 = xvec(ip1,j) * yvec(ip2,ip1) + yvec(ip1,j) * xvec(ip1,ip2)
!              ***  calculate cross product d_a-1 x d_a-2 ***
                        xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) + zvec(ip1,ip2) * yvec(ip3,ip2)
                        ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) + xvec(ip1,ip2) * zvec(ip3,ip2)
                        za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) + yvec(ip1,ip2) * xvec(ip3,ip2)
                              
!     *** calculate lengths of cross products ***
                        daa1 = sqrt(xaa1**2+yaa1**2+zaa1**2)
                        da1a2 = sqrt(xa1a2**2+ya1a2**2+za1a2**2)
!              *** calculate dot product of cross products ***
                        dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                        thetac = - dot / ( daa1 * da1a2 )
                        theta = acos(thetac)

!              *** calculate cross product of cross products ***
                        xcc = yaa1*za1a2 - zaa1*ya1a2
                        ycc = zaa1*xa1a2 - xaa1*za1a2
                        zcc = xaa1*ya1a2 - yaa1*xa1a2
!              *** calculate scalar triple product ***
                        tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2)  + zcc*zvec(ip1,ip2) 
                        if ( tcc .lt. 0.0d0 ) theta = - theta

!              *** convert theta to degrees ***
                        theta = theta*(180.0d0/onepi)

!                       --- bin the torsion --
                        bin = aint( (theta+180.0d0)/tor_bin_size)+1
                        tor_bin(kk,imolty,torsion,bin) =  tor_bin(kk,imolty,torsion,bin)+1.0d0
                        tor_tot(kk,imolty,torsion) =  tor_tot(kk,imolty,torsion) + 1.0d0

                        if ( theta .lt. -60.0d0 ) then
!                          --- gauch minus
                           g_minus(kk,imolty)=g_minus(kk,imolty) + 1.0d0
                           gaudef = gaudef + 1
                           bg_minus(kk,imolty,torsion) =  bg_minus(kk,imolty,torsion) + 1.0d0
                           dum = dum + 0 * 3**(torsion-1)
                        else if ( theta .lt. 60.0d0 ) then
!                          --- trans
                           trans(kk,imolty) = trans(kk,imolty)+1.0d0
                           btrans(kk,imolty,torsion) =  btrans(kk,imolty,torsion) + 1.0d0
                           dum = dum + 1 * 3**(torsion-1)
                        else
!                          --- gauch plus
                           g_plus(kk,imolty)=g_plus(kk,imolty) + 1.0d0
                           gaudef = gaudef + 1
                           bg_plus(kk,imolty,torsion) =  bg_plus(kk,imolty,torsion) + 1.0d0
                           dum = dum + 2 * 3**(torsion-1)
                        end if

                     end do

                     gdefect(kk,imolty,gaudef) =  gdefect(kk,imolty,gaudef) + 1.0d0 
                     gdefect(kk,imolty,tor_max+1) =  gdefect(kk,imolty,tor_max+1) + 1.0d0

!  Being removed because of high memory requirement
!                     pattern(kk,imolty,dum) = pattern(kk,imolty,dum) 
!     &                    + 1.0d0

                  end if
               end if
              end do
           end do
        end if

      return
      end if

      if(switch.eq.2) then
 
        fname = run_num
        write(ftemp,*) fname
        read(ftemp,*) fname3


      if(lete) then
         do i=1,nbox
!             fname = i
             write(ftemp,*) i 
             read(ftemp,*) fname2 
             fileout = 'end2end_box'//fname2(1:len_trim(fname2))
             open (140+i,FILE=fileout,STATUS="unknown")
          end do



!         open (142,FILE="end2end_box2",STATUS="unknown")
!         open (143,FILE="end2end_box3",STATUS="unknown")

         do kk=1,nbox
           do imolty=1,nmolty
              xx =  (aint(max_length(imolty)/bin_width)+1)
              write(140+kk,*) 'molecule type', imolty
              write(140+kk,*) 
              do bin=1,xx
                  rxuij=bin_width*(dble(bin)-0.5d0)
                  write(140+kk,*) rxuij,end_to_end(kk,imolty,bin)/ nframe
              end do
           end do
         end do              
         
         do i = 1,nbox          
            close(140+i)
         end do 


!         close(142)
!         close(143)

      end if


      if(lrhoz) then
         
         do i=1,nbox
             write(ftemp,*) i
             read(ftemp,*) fname2
             fileout = 'rhoz_box'//fname2(1:len_trim(fname2))
             open (150+i,FILE=fileout,STATUS="unknown")
             fileout = 'comrhoz_box'//fname2(1:len_trim(fname2))
             open (153+i,FILE=fileout,STATUS="unknown") 
         end do
                     
         do kk=1,nbox
           do imolty=1,nmolty
              xx =  (aint(max_boxlz(kk)/bin_width)+1)
              write(150+kk,*) 'molecule type',imolty
              write(153+kk,*) 'molecule type',imolty
              write(153+kk,*) 
              write(150+kk,*)
              do bin=1,xx
                  rxuij=bin_width*(dble(bin)-0.5d0)
                  write(150+kk,*) rxuij,bigrhoz(kk,imolty,bin)/ nframe
              end do
              do bin=-xx,xx
                 if (bin.eq.0) then
                     rxuij=0.0
                 else
                     rxuij=bin_width*(dble(bin)+0.5d0) 
                 end if
                 write(153+kk,*) rxuij,bigboxcom_rhoz(kk,imolty,bin)/ nframe
              end do 
           end do
         end do

         close(151)
         close(152)
         close(153)
         close(154)
         close(155)
         close(156)

      end if


!     whew, now that is done we just have to divide out the number of
!     frames in biganalhist


      if (lrdf) then

         do i=1,nbox
             write(ftemp,*) i
             read(ftemp,*) fname2
             fileout = 'beadrdf_box'//fname2(1:len_trim(fname2))// '_'//fname3(1:len_trim(fname3))//suffix//'.dat'
             open (100+i,FILE=fileout,STATUS="unknown")
             fileout = 'comrdf_box'//fname2(1:len_trim(fname2))// '_'//fname3(1:len_trim(fname3))//suffix//'.dat'
             open (103+i,FILE=fileout,STATUS="unknown")
             fileout = 'beadnum_box'//fname2(1:len_trim(fname2))// '_'//fname3(1:len_trim(fname3))//suffix//'.dat'
             open (110+i,FILE=fileout,STATUS="unknown")
             fileout = 'comnum_box'//fname2(1:len_trim(fname2))// '_'//fname3(1:len_trim(fname3))//suffix//'.dat'
             open (113+i,FILE=fileout,STATUS="unknown")
         end do
         

!         write(io_output,*)
!         write(io_output,*) 'bead bead radial distribution functions in *gor*'
!         write(io_output,*)
!         write(io_output,*) 'avg. number of beads vs. shell size in *num*'
         
!         print*, 'nhere', nhere
         
      do kk = 1,nbox
        binstep = rcut(kk)/dble(nbin)
         do ii = 1,nmolty
            do xx = 1,nhere
               ntii = nhere*(ii-1)+xx
               do jj = 1,nmolty
                  do yy = 1,nhere
                   ntjj = nhere*(jj-1)+yy
                   if ( ntii .le. ntjj ) then 
                      ntij = (ntii-1)*nhere*nmolty + ntjj
                   
                      binadj = (kk-1)*nhere*nhere*nmolty*nmolty + ntij

                      aframe = nframe-nnone(kk,ntij)

                      if ( aframe .gt. 0.5d0 ) then

                         write(100+kk,'(2f7.2,4i5)')0.0d0, 0.0d0,ii ,beadtyp(xx),jj,beadtyp(yy) 

                         write(110+kk,'(2f7.2,4i5)') 0.0d0,0.0d0,ii ,beadtyp(xx),jj,beadtyp(yy) 
                      do bin = 1,nbin

                            rxuij =  binstep*(dble(bin)-0.5d0)
 
                            biganalhist(binadj,bin) = biganalhist(binadj ,bin)/aframe

                            write(100+kk,*) rxuij,biganalhist(binadj, bin)

                             shell(binadj,bin,1) = shell(binadj,bin,1) /aframe
                         write(110+kk,*) rxuij,shell(binadj,bin,1)
                         
                      end do
                      write(100+kk,*)
	              write(110+kk,*)
 
                      if (ntii .ne. ntjj) then
                         write(110+kk,'(2f7.2,4i5)') 0.0d0,0.0d0,jj ,beadtyp(yy),ii,beadtyp(xx)
                          do bin = 1,nbin
                              rxuij =  binstep*(dble(bin)-0.5d0)
                              shell(binadj,bin,2) = shell(binadj,bin,2) /aframe
                              write(110+kk,*) rxuij,shell(binadj,bin,2)
                            end do
                            write(110+kk,*)
                      end if
                      
                   else if( aframe .lt. -0.5d0 ) then
                      write(io_output,*) 'aframe',aframe
                      write(io_output,*) 'nframe',nframe
                      write(io_output,*) 'nnone(kk,ntij),kk,ntij', nnone(kk,ntij),kk,ntij
                      call cleanup('srewup aframe')
                     end if
                   end if
                 end do
               end do
            end do
         end do
       end do

!     --- same thing for the COM rdf
       do kk = 1,nbox
         binstep = rcut(kk)/dble(nbin)
         do ii = 1,nmolty
            ntii = ii
            do jj = 1,nmolty

               ntjj = jj
               ntij = (ntii-1)*nmolty + ntjj

               binadj = (kk-1)*nmolty*nmolty + ntij
               aframe = nframe-comnone(kk,ntij)

               if ( aframe .gt. 0.5d0 ) then
                  write(103+kk,'(2f7.2,2i5)') 0.0d0,0.0d0,ii,jj
                  write(113+kk,'(2f7.2,2i5)') 0.0d0,0.0d0,ii,jj
                  do bin = 1,nbin
                     rxuij =  binstep*(dble(bin)-0.5d0)
                     
                     combiganalhist(binadj,bin) = combiganalhist(bina dj,bin)/aframe
                     write(103+kk,*) rxuij,combiganalhist(binadj,bin)

                     comshell(binadj,bin,1) = comshell(binadj,bin,1) /aframe
                     write(113+kk,*) rxuij,comshell(binadj,bin,1)

                  end do

                  write(103+kk,*)
                  write(113+kk,*)
                  
                  if(ntii.ne.ntjj) then  
                    write(113+kk,'(2f7.2,2i5)') 0.0d0,0.0d0,jj,ii
                    do bin = 1,nbin
                       rxuij =  binstep*(dble(bin)-0.5d0)
                       comshell(binadj,bin,2) = comshell(binadj,bin,2) /aframe
                       write(113+kk,*) rxuij,comshell(binadj,bin,2)
	            end do
                    write(113+kk,*)
                  end if    
               else if( aframe .lt. -0.5d0 ) then
                  write(io_output,*) 'aframe',aframe
                  write(io_output,*) 'nframe',nframe
                  write(io_output,*) 'comnone(kk,ntij),kk,ntij', comnone(kk,ntij),kk,ntij
                  call cleanup('srewup aframe')
               end if
            end do
         end do
       end do

       close(101)
       close(102)
       close(103)
       close(104)
       close(105)
       close(106)
       close(111)
       close(112)
       close(113)
       close(114)
       close(115)
       close(116)
      end if 

!      if ( lcharge ) then
!         write(io_output,*) 'charge distributions in fort.35+box'
!         do kk = 1, nbox
!            do imol = 1,nmolty
!               do iunit = 1,nunit(imolty)
!                  qqcode = numax*(imol-1) + iunit
!                  qhere = .false.
!                  do qbin = 1,qbins
!                     if ( qcount(qqcode,kk) .gt. 0.5d0 ) then
!                        qhere = .true.
!                        qanalhist(qqcode,kk,qbin) = 
!     &                       qanalhist(qqcode,kk,qbin)/qcount(qqcode,kk)
!                     end if
!                  end do
!                  if (qhere) then
!                     qdummy = ( (dble(1)-0.5d0)*qstep)+qmin
!                     write(35+kk,*) qdummy,qanalhist(qqcode,kk,1)+qdisp
!     &                    ,imol,iunit
!                     
!                     do qbin = 2,qbins
!                        qdummy = ( (dble(qbin)-0.5d0)*qstep)+qmin
!                        write(35+kk,*) qdummy
!     &                       ,qanalhist(qqcode,kk,qbin)+qdisp
!                     end do
!                     write(35+kk,*)
!                  end if
!               end do
!            end do
!         end do
!      end if

!      if ( ldipole ) then
!         diconv = (1.602d-19)/((1d10)*(3.336d-30))
!         do kk = 1, nbox
!            write(io_output,*) 'Dipoles [e A], [D] in Box',kk
!            do imol = 1,nmolty
!               if ( dicount(imol,kk) .gt. 0.5d0 ) then
!                  dipole(imol,kk) = dipole(imol,kk)/dicount(imol,kk)
!               end if
!               write(io_output,*) 'Moltyp,dipole',imolty
!     &              ,dipole(imol,kk),dipole(imol,kk)*diconv
!            end do
!         end do
!
!         if ( block .gt. 1 ) then
!            write(io_output,*) 'Number of blocks input and used',block,nblock
!            do kk = 1,nbox
!               write(io_output,*) 'Box     Moltyp   Dipole [D]   Std dev [D]'
!               do imol = 1,nmolty
!                  avera = 0.0d0
!                  do k=1,nblock
!                     avera = avera + dipblk(imol,kk,k)
!                  end do
!                  avera = avera/nblock
!                  stddev = 0.0d0
!                  do k = 1,nblock
!                     stddev = stddev + (dipblk(imol,kk,k)-avera)**2
!                  end do
!                  stddev = sqrt( stddev/dble(nblock-1) )
!                  stddev = stddev*diconv
!                  avera = avera*diconv
!                  write(io_output,'(i3,4x,i3,4x,2f10.4)') kk,imol,avera,stddev
!               end do
!            end do
!         end if

!      end if
               
      if ( lbend ) then
         do i=1,nbox
             write(ftemp,*) i
             read(ftemp,*) fname2
             fileout = 'bendang_dist_box'//fname2
             open (120+i,FILE=fileout,STATUS="unknown")
         end do

!     --- output the bending angle distributions for each angle in the molecule
!     --- write to 37+box

         do kk = 1,nbox
            do imolty = 1,nmolty
               do bend = 1,angle_num(imolty)
                  total = angle_tot(kk,imolty,bend)
                  if ( total .gt. 0.5d0 ) then
                     value = ang_bin_size/2.0d0
                     degree = value*180.0d0/onepi
                     write(120+kk,*) degree ,angle_bin(kk,imolty,bend,1)/total ,'Moltyp ',imolty,' angle ' ,angle_1(imolty,bend),angle_2(imolty,bend) ,angle_3(imolty,bend)
                     do bin = 2,ang_bin_max
                        value = value + ang_bin_size
                        degree = value*180.0d0/onepi
                        write(120+kk,*) degree ,angle_bin(kk,imolty,bend,bin)/total
                     end do
                     write(120+kk,*)
!                    --- output the average angle to the screen
                     value = angle_avg(kk,imolty,bend)/total
                     degree = value*180.0d0/onepi
                     write(io_output,*)  'Moltyp ',imolty,' angle ' ,angle_1(imolty,bend),angle_2(imolty,bend) ,angle_3(imolty,bend),' average ',degree
                     
                  end if
               end do
            end do
         end do                
         close(121)
         close(122)
         close(123)
      end if

      if ( lgvst ) then

         do i=1,nbox
             write(ftemp,*) i
             read(ftemp,*) fname2
             fileout = 'tor_frac_box'//fname2
             open (130+i,FILE=fileout,STATUS="unknown")
             fileout = 'torsprob_box'//fname2
             open (133+i,FILE=fileout,STATUS="unknown")
         end do
 


!     output the gauch versus trans data for each moltype
         write(io_output,*) 'moltyp  box  trans      g+        g-       g frac'
         do kk  = 1,nbox
            do imolty = 1,nmolty
               write(130+kk,*) 'Moltype ',imolty
               write(130+kk,"('Units',8x,'Frac g+',1x,'Frac g-',1x ,'Frac t')")
               gaudef = g_plus(kk,imolty)+g_minus(kk,imolty)
       
               if((aint(trans(kk,imolty))+gaudef).ne.0) then
                  ratio = gaudef/(trans(kk,imolty)+gaudef)
               else
!                  ratio = 0.0d0
                  write(io_output,*) 'ratio not defined'
               end if
               
               write(io_output,'(i3,5x,i3,3(2x,e8.2),f8.4)')imolty,kk ,trans(kk ,imolty),g_plus(kk,imolty),g_minus(kk,imolty),ratio

               do torsion = 1,tor_num(imolty)
                  total = bg_plus(kk,imolty,torsion) +  bg_minus(kk,imolty,torsion) + btrans(kk,imolty,torsion)

                  if ( total .gt. 0.5d0) then
!                    --- output torsion type fractions 
                     fplus = bg_plus(kk,imolty,torsion)/total
                     fminus = bg_minus(kk,imolty,torsion)/total
                     ftrans = btrans(kk,imolty,torsion)/total
                     write(130+kk,'(4(i2,1x),1x,3(f5.3,3x))') tor_1(imolty,torsion),tor_2(imolty,torsion) ,tor_3(imolty,torsion),tor_4(imolty,torsion),fplus ,fminus,ftrans

!                    --- output the torsion probablility distributions
                     value = -180.0d0 + (0.5d0*tor_bin_size)
                     write(133+kk,*) value,tor_bin(kk,imolty,torsion,1) ,imolty,tor_1(imolty,torsion) ,tor_2(imolty,torsion),tor_3(imolty,torsion) ,tor_4(imolty,torsion)

                     do bin = 2,tor_bin_max
                        value = value + tor_bin_size
                        write(133+kk,*) value ,tor_bin(kk,imolty,torsion,bin)
                     end do
                     write(133+kk,*)
                  end if
               end do
            end do
         end do

         close(131)
         close(132)
         close(133)
         close(134)
         close(135)
         close(136)

!         write(io_output,*)
!         write(io_output,*) 'gauch fraction vs. torsion in fort.90+box'
!         write(io_output,*)
!         write(io_output,*) 'torsion angle distribution in fort.95+box'
!         write(io_output,*)
!         write(io_output,*) 'probabiltiy vs. number of defects per chain'
!         write(io_output,*) 'shown in fort.5*'
!         write(io_output,*)
         
         do i=1,nbox
             write(ftemp,*) i
             read(ftemp,*) fname2
             fileout = 'Gauchedefects_box'//fname2
             open (136+i,FILE=fileout,STATUS="unknown")
!             fileout = 'pattern_box'//fname2
!             open (124+i,FILE=fileout,STATUS="unknown")
         end do

         open(124,FILE="decoder",status="unknown")



!     analyse the number of gauch defects per chain
         do kk = 1,nbox
            do imolty = 1,nmolty
               total = gdefect(kk,imolty,tor_max+1)
               if ( total .gt. 0.01d0 ) then
                  write(136+kk,*) 'molecule type ',imolty
                  do torsion = 1,tor_num(imolty)
                     gdefect(kk,imolty,torsion) =  gdefect(kk,imolty,torsion)/total
                     write(136+kk,*) torsion,gdefect(kk,imolty,torsion)
                  end do
               end if

               psum = 0.0d0
!               if ( tor_num(imolty) .gt. 0 ) then
!c                 --- write out a decoder for the torsions
!                  if ( kk .eq. 1 ) then
!                     write(124,*) 'Moltyp ',imolty,' torsions'
!                     do torsion = 1,tor_num(imolty)
!                        write(124,*) torsion,'   ',tor_1(imolty,torsion)
!     &                       ,tor_2(imolty,torsion)
!     &                       ,tor_3(imolty,torsion)
!     &                       ,tor_4(imolty,torsion)
!                     end do
!                  end if
!                  units = 3**( tor_num(imolty) )
!                  do uu = 0,units
!                     psum = psum + pattern(kk,imolty,uu)
!                  end do
!                  write(124+kk,*) 'Code     Prob.     Torsion pattern'
!                  do uu = 0,units-1
!                     if ( pattern(kk,imolty,uu) .gt. 0.5d0 ) then
!                        decimal = uu
!                        power = units
!                        torsion = tor_num(imolty)
!                        do dummy = torsion,1,-1
!                           power = power/3
!                           bthree = decimal/power 
!                           decimal = decimal - bthree*power
!                           patt(dummy) = bthree-1
!                        end do

!                        write(124+kk,'(i8,2x,f5.3,2x,20(i2,1x))') uu,
!     &                       pattern(kk,imolty,uu)/psum
!     &                       ,(patt(dummy),dummy=1,torsion)
!                     end if
!                  end do 
!               end if

            end do
         end do
         close(137)
         close(138)
         close(139)
!         close(124)
!         close(125)
!         close(126)
!         close(127)

!         write(io_output,*) 'patterns of gauch defects in fort.2*'
!         write(io_output,*) 'transform first number into base 3 to get the'
!         write(io_output,*) 'pattern of gauch defects where -1=g-, 0=trans 1=g+'
!         write(io_output,*)
!         write(io_output,*) 'the maximum rcut value that could be used is'
!         write(io_output,*) boxmin
!         write(io_output,*)
         

      end if
      return 
      end if
       
      end subroutine analysis
