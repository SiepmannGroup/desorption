MODULE energy_kspace
  use var_type,only:dp
  use const_math,only:onepi,twopi
  use const_phys,only:qqfact
  use util_runtime,only:err_exit
  use util_mp,only:mp_sum,mp_allgather,mp_set_displs,mp_bcast
  use sim_system,only:lsolid,lrect,lrecip,boxlx,boxly,boxlz,nchain,moltyp,lelect,nboxi,ntype,lqchg,rxu,ryu,rzu,qqu,myid,numprocs,moltion&
   ,nunit,rxuion,ryuion,rzuion,qquion,xcm,ycm,zcm,groupid,rcut,rootid
  use sim_cell
  implicit none
  private
  save
  public::recipsum,recip,recip_atom,ee_recip,recippress,calp,sself,correct,save_kvector,restore_kvector,allocate_kspace
  public::k_max_l, k_max_m, k_max_n, compute_kmax

  integer,parameter::vectormax=100000 !< the maximum number of reciprocal vectors for Ewald sum
  integer,parameter::vectormax1=100 !< the maximum number of reciprocal vectors in each dimension
  integer,allocatable::numvect(:)& !< the total number of reciprocal vectors
   ,numvecto(:),k_max_l(:),k_max_m(:),k_max_n(:),k_max_lo(:),k_max_mo(:),k_max_no(:)
  real,allocatable::kx(:,:),ky(:,:),kz(:,:),prefact(:,:),ssumr(:,:),ssumi(:,:),ssumrn(:,:),ssumin(:,:),ssumro(:,:),ssumio(:,:)&
   ,kxo(:,:),kyo(:,:),kzo(:,:),prefacto(:,:),calpo(:)
  real,allocatable,target::calp(:) !< calp = kalp / boxlen; kalp is a parameter to control the real space sum
  real::sself,correct
  integer,allocatable::nN0(:),nNl(:,:,:),nNlm(:,:,:,:),nN0o(:),nNlo(:,:,:),nNlmo(:,:,:,:),nnc(:,:),nnco(:,:)

contains
!> \brief calculates the total reciprocal space ewald-sum term for volume
!> moves
!> \par History
!> written in 1998 by Bin Chen \n
!> rewritten in 2001 by Bin Chen \n
!> rewritten again, probably by Bin
!> Bai Xue 2016: Enable efficient calculation of cos and sin if lrecip=.true. 
  subroutine recipsum(ibox,vrecip)
    real(kind=dp)::vrecip
    integer::ibox,i,ii,imolty,ncount,ncount1,m_pre,n_pre

    ! from h-matrix formulation
    integer::l,m,n,m_min,n_min,kmaxl,kmaxm,kmaxn,N0,Nl(100,2),Nlm(100,100,4),nm,nn

    real::alpsqr4,vol,ksqr,sumr,sumi,arg,bx1,by1,bz1,hmatik(9),kx1,ky1,kz1,hmaxsq,calpi
    ! real::sum_sumr,sum_sumi

    ! RP added for calculating time for communication step
    integer::mystart,myend,blocksize
    integer::rcounts(numprocs),displs(numprocs)
    real::my_kx(vectormax),my_ky(vectormax),my_kz(vectormax),my_ssumr(vectormax),my_ssumi(vectormax),my_prefact(vectormax)
    real::A(vectormax1,vectormax1,vectormax1),B(vectormax1,vectormax1,vectormax1),sumr1(vectormax),sumi1(vectormax)

    ! Set up the reciprocal space vectors ***
    ncount = 0
    vrecip = 0.0E0_dp

    calpi = calp(ibox)

    if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
       bx1 = boxlx(ibox)
       by1 = boxly(ibox)
       bz1 = boxlz(ibox)
       hmat(ibox,1) = bx1
       hmat(ibox,5) = by1
       hmat(ibox,9) = bz1
       do i = 1,9
          hmatik(i) = 0.0E0_dp
       end do
       hmatik(1) = twopi/bx1
       hmatik(5) = twopi/by1
       hmatik(9) = twopi/bz1
    else
       do i = 1,9
          hmatik(i) = twopi*hmati(ibox,i)
       end do
    end if
    call compute_kmax(ibox)
    kmaxl = k_max_l(ibox)
    kmaxm = k_max_m(ibox)
    kmaxn = k_max_n(ibox)

    alpsqr4 = 4.0E0_dp*calpi*calpi

    vol = hmat(ibox,1)* (hmat(ibox,5)*hmat(ibox,9) - hmat(ibox,8)*hmat(ibox,6))&
     + hmat(ibox,4)* (hmat(ibox,8)*hmat(ibox,3) - hmat(ibox,2)*hmat(ibox,9))&
     + hmat(ibox,7)* (hmat(ibox,2)*hmat(ibox,6) - hmat(ibox,5)*hmat(ibox,3))

    vol = vol/(4.0E0_dp*onepi)

    hmaxsq = alpsqr4*calpi*rcut(ibox)*calpi*rcut(ibox)

    if(lrecip .and. (.not. lsolid(ibox)) .or. lrect(ibox)) then

       if(kmaxl.gt.vectormax1.or.kmaxm.gt.vectormax1.or.kmaxn.gt.vectormax1) call err_exit(__FILE__,__LINE__,'choose a larger vectormax1',myid+1)

       ! First find out which k vectors will be calculated and construct their structure
       ncount = 0
       N0 = 0
       Nl = 0
       Nlm = 0
       do l = 0,kmaxl
      
          if( l .gt. 0) then
             nnc(l,ibox) = ncount
             kx1 = real(l,dp)*hmatik(1)
             ky1 = real(l,dp)*hmatik(4)
             kz1 = real(l,dp)*hmatik(7)
             ksqr = kx1*kx1+ky1*ky1+kz1*kz1
             if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq) .gt. 1E-9_dp ) then
                ncount = ncount +1
                N0 = N0 + 1
                kx(ncount,ibox) = kx1
                ky(ncount,ibox) = ky1
                kz(ncount,ibox) = kz1
                prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
             end if
          end if

          do m = 0, kmaxm

             if( m .gt. 0) then
                kx1 = real(l,dp)*hmatik(1)+real(m,dp)*hmatik(2)
                ky1 = real(l,dp)*hmatik(4)+real(m,dp)*hmatik(5)
                kz1 = real(l,dp)*hmatik(7)+real(m,dp)*hmatik(8)
                ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq) .gt. 1E-9_dp ) then
                   ncount = ncount + 1
                   Nl(l+1,1) = Nl(l+1,1) + 1
                   kx(ncount,ibox) = kx1
                   ky(ncount,ibox) = ky1
                   kz(ncount,ibox) = kz1
                   prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                end if
             end if

             do n = 1, kmaxn

                kx1 = real(l,dp)*hmatik(1)+real(m,dp)*hmatik(2)+real(n,dp)*hmatik(3)
                ky1 = real(l,dp)*hmatik(4)+real(m,dp)*hmatik(5)+real(n,dp)*hmatik(6)
                kz1 = real(l,dp)*hmatik(7)+real(m,dp)*hmatik(8)+real(n,dp)*hmatik(9)
                ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq) .gt. 1E-9_dp ) then
                   ncount = ncount + 1
                   Nlm(l+1,m+1,1) = Nlm(l+1,m+1,1) + 1
                   kx(ncount,ibox) = kx1
                   ky(ncount,ibox) = ky1
                   kz(ncount,ibox) = kz1
                   prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                end if

             end do

             if (l.ne.0.or.m.ne.0) then

                do n = 1,kmaxn
                   kx1 = real(l,dp)*hmatik(1)+real(m,dp)*hmatik(2)-real(n,dp)*hmatik(3)
                   ky1 = real(l,dp)*hmatik(4)+real(m,dp)*hmatik(5)-real(n,dp)*hmatik(6)
                   kz1 = real(l,dp)*hmatik(7)+real(m,dp)*hmatik(8)-real(n,dp)*hmatik(9)
                   ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                   if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq) .gt. 1E-9_dp ) then
                      ncount = ncount + 1
                      Nlm(l+1,m+1,2) = Nlm(l+1,m+1,2) + 1
                      kx(ncount,ibox) = kx1
                      ky(ncount,ibox) = ky1
                      kz(ncount,ibox) = kz1
                      prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                   end if
                end do

             end if
          end do

          if (l.ne.0) then
             do m = 1,kmaxm

                 kx1 = real(l,dp)*hmatik(1)-real(m,dp)*hmatik(2)
                 ky1 = real(l,dp)*hmatik(4)-real(m,dp)*hmatik(5)
                 kz1 = real(l,dp)*hmatik(7)-real(m,dp)*hmatik(8)
                 ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                 if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq).gt.1E-9_dp) then
                    ncount = ncount + 1
                    Nl(l+1,2) = Nl(l+1,2) + 1
                    kx(ncount,ibox) = kx1
                    ky(ncount,ibox) = ky1
                    kz(ncount,ibox) = kz1
                    prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                 end if

                 do n=1,kmaxn

                     kx1 = real(l,dp)*hmatik(1)-real(m,dp)*hmatik(2)+real(n,dp)*hmatik(3)
                     ky1 = real(l,dp)*hmatik(4)-real(m,dp)*hmatik(5)+real(n,dp)*hmatik(6)
                     kz1 = real(l,dp)*hmatik(7)-real(m,dp)*hmatik(8)+real(n,dp)*hmatik(9)
                     ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                     if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq).gt.1E-9_dp) then
                        ncount = ncount + 1
                        Nlm(l+1,m+1,3) = Nlm(l+1,m+1,3) + 1
                        kx(ncount,ibox) = kx1
                        ky(ncount,ibox) = ky1
                        kz(ncount,ibox) = kz1
                        prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                     end if
                   
                 end do
          
                 do n = 1,kmaxn

                    kx1 = real(l,dp)*hmatik(1)-real(m,dp)*hmatik(2)-real(n,dp)*hmatik(3)
                    ky1 = real(l,dp)*hmatik(4)-real(m,dp)*hmatik(5)-real(n,dp)*hmatik(6)
                    kz1 = real(l,dp)*hmatik(7)-real(m,dp)*hmatik(8)-real(n,dp)*hmatik(9)
                    ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                    if ( ksqr .lt. hmaxsq .and.abs(ksqr-hmaxsq).gt.1E-9_dp) then
                       ncount = ncount + 1
                       Nlm(l+1,m+1,4) = Nlm(l+1,m+1,4) + 1
                       kx(ncount,ibox) = kx1
                       ky(ncount,ibox) = ky1
                       kz(ncount,ibox) = kz1
                       prefact(ncount,ibox) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                    end if
                 end do

             end do
          end if
       end do

       nN0(ibox) = N0
       nNl(:,:,ibox) = Nl
       nNlm(:,:,:,ibox) = Nlm

       ! Using cos(A+B)=cosAcosB-sinAsinB and sin(A+B)=sinAcosB+cosAsinB to calculate all cos and sin 
       blocksize = nchain/numprocs
       mystart = myid * blocksize + 1
       if (myid .eq. (numprocs-1)) then
          myend = nchain
       else
          myend = (myid + 1) * blocksize 
       end if
       vrecip = 0.0E0_dp
       sumr1 = 0.0E0_dp
       sumi1 = 0.0E0_dp
       do i = mystart,myend
          imolty = moltyp(i)
          if(.not.lelect(imolty).or.nboxi(i).ne.ibox) cycle
          do ii = 1,nunit(imolty)
             if ( lqchg(ntype(imolty,ii)) .and. abs(qqu(i,ii)).gt.1e-6) then

                ncount1 = 0
                kx1 = hmatik(1)*rxu(i,ii)+hmatik(4)*ryu(i,ii)+hmatik(7)*rzu(i,ii)
                ky1 = hmatik(2)*rxu(i,ii)+hmatik(5)*ryu(i,ii)+hmatik(8)*rzu(i,ii)
                kz1 = hmatik(3)*rxu(i,ii)+hmatik(6)*ryu(i,ii)+hmatik(9)*rzu(i,ii)
                A(1,1,1) = 1.0E0_dp
                B(1,1,1) = 0.0E0_dp
                A(2,1,1) = cos(kx1)
                B(2,1,1) = sin(kx1)
                A(1,2,1) = cos(ky1)
                B(1,2,1) = sin(ky1)
                A(1,1,2) = cos(kz1)
                B(1,1,2) = sin(kz1)
                do l = 0,N0

                   if( l .gt. 0) then
                      ncount1 = ncount1 + 1
                      A(l+1,1,1) = A(l,1,1)*A(2,1,1)-B(l,1,1)*B(2,1,1)
                      B(l+1,1,1) = B(l,1,1)*A(2,1,1)+A(l,1,1)*B(2,1,1)
                      sumr1(ncount1) = sumr1(ncount1) + A(l+1,1,1)*qqu(i,ii)
                      sumi1(ncount1) = sumi1(ncount1) + B(l+1,1,1)*qqu(i,ii)
                   end if

                   do m = 0, Nl(l+1,1)

                      if( m .gt. 0) then
                         ncount1 = ncount1 + 1
                         A(l+1,m+1,1) = A(l+1,m,1)*A(1,2,1)-B(l+1,m,1)*B(1,2,1)
                         B(l+1,m+1,1) = B(l+1,m,1)*A(1,2,1)+A(l+1,m,1)*B(1,2,1)                  
                         sumr1(ncount1) = sumr1(ncount1) + A(l+1,m+1,1)*qqu(i,ii)
                         sumi1(ncount1) = sumi1(ncount1) + B(l+1,m+1,1)*qqu(i,ii)
                      end if

                      do n = 1, Nlm(l+1,m+1,1)

                         ncount1 = ncount1 + 1
                         A(l+1,m+1,n+1) = A(l+1,m+1,n)*A(1,1,2)-B(l+1,m+1,n)*B(1,1,2)
                         B(l+1,m+1,n+1) = B(l+1,m+1,n)*A(1,1,2)+A(l+1,m+1,n)*B(1,1,2)
                         sumr1(ncount1) = sumr1(ncount1) + A(l+1,m+1,n+1)*qqu(i,ii)
                         sumi1(ncount1) = sumi1(ncount1) + B(l+1,m+1,n+1)*qqu(i,ii)

                      end do

                      if (l.ne.0.or.m.ne.0) then

                         do n = 1,Nlm(l+1,m+1,2)
                            ncount1 = ncount1 + 1
                            nn = Nlm(l+1,m+1,2)
                            if(n .eq. 1) then
                               n_pre = 0
                            else
                               n_pre = nn
                            end if
                            A(l+1,m+1,nn+n+1) = A(l+1,m+1,n_pre+n)*A(1,1,2)+B(l+1,m+1,n_pre+n)*B(1,1,2)
                            B(l+1,m+1,nn+n+1) = B(l+1,m+1,n_pre+n)*A(1,1,2)-A(l+1,m+1,n_pre+n)*B(1,1,2)
                            sumr1(ncount1) = sumr1(ncount1) + A(l+1,m+1,nn+n+1)*qqu(i,ii)
                            sumi1(ncount1) = sumi1(ncount1) + B(l+1,m+1,nn+n+1)*qqu(i,ii)
                         end do
                      end if
                   end do

                   if (l.ne.0) then
                      do m = 1,Nl(l+1,2)
                         ncount1 = ncount1 + 1
                         nm = Nl(l+1,2)
                         if(m .eq. 1) then
                            m_pre = 0
                         else
                            m_pre = nm
                         end if
                         A(l+1,nm+m+1,1) = A(l+1,m_pre+m,1)*A(1,2,1)+B(l+1,m_pre+m,1)*B(1,2,1)
                         B(l+1,nm+m+1,1) = B(l+1,m_pre+m,1)*A(1,2,1)-A(l+1,m_pre+m,1)*B(1,2,1)
                         sumr1(ncount1) = sumr1(ncount1) + A(l+1,nm+m+1,1)*qqu(i,ii)
                         sumi1(ncount1) = sumi1(ncount1) + B(l+1,nm+m+1,1)*qqu(i,ii)

                         do n = 1, Nlm(l+1,m+1,3)
                            ncount1 = ncount1 + 1
                            A(l+1,nm+m+1,n+1) = A(l+1,nm+m+1,n)*A(1,1,2)-B(l+1,nm+m+1,n)*B(1,1,2)
                            B(l+1,nm+m+1,n+1) = B(l+1,nm+m+1,n)*A(1,1,2)+A(l+1,nm+m+1,n)*B(1,1,2)
                            sumr1(ncount1) = sumr1(ncount1) + A(l+1,nm+m+1,n+1)*qqu(i,ii)
                            sumi1(ncount1) = sumi1(ncount1) + B(l+1,nm+m+1,n+1)*qqu(i,ii)
                         end do

                         do n = 1,Nlm(l+1,m+1,4)
                            ncount1 = ncount1 + 1
                            nn = Nlm(l+1,m+1,4)
                            if(n .eq. 1) then
                               n_pre = 0
                            else
                               n_pre = nn
                            end if
                            A(l+1,nm+m+1,nn+n+1) = A(l+1,nm+m+1,n_pre+n)*A(1,1,2)+B(l+1,nm+m+1,n_pre+n)*B(1,1,2)
                            B(l+1,nm+m+1,nn+n+1) = B(l+1,nm+m+1,n_pre+n)*A(1,1,2)-A(l+1,nm+m+1,n_pre+n)*B(1,1,2)
                            sumr1(ncount1) = sumr1(ncount1) + A(l+1,nm+m+1,nn+n+1)*qqu(i,ii)
                            sumi1(ncount1) = sumi1(ncount1) + B(l+1,nm+m+1,nn+n+1)*qqu(i,ii)
                         end do
                      end do
                   end if
                end do
             end if              
          end do
       end do

       call mp_sum(sumr1,ncount,groupid)
       call mp_sum(sumi1,ncount,groupid)     

       do i=1,ncount
          vrecip = vrecip + (sumr1(i)*sumr1(i)+sumi1(i)*sumi1(i)) * prefact(i,ibox)
          ssumr(i,ibox) = sumr1(i)
          ssumi(i,ibox) = sumi1(i)
       end do
       numvect(ibox) = ncount

       return

    else

       ! RP added for MPI
       blocksize = kmaxl/numprocs
       mystart = myid * blocksize
       if (myid .eq. (numprocs-1)) then
          myend = kmaxl
       else
          myend = (myid + 1) * blocksize - 1
       end if
       ! generate the reciprocal-space
       ! here -kmaxl,-kmaxl+1,...,-1 are skipped, so no need to divide by 2 for the prefactor
       do l = mystart,myend
       ! do l = 0,kmaxl
          if ( l .eq. 0 ) then
             m_min = 0
          else
             m_min = -kmaxm
          end if
          do m = m_min, kmaxm
             if (l .eq. 0 .and. m .eq. 0) then
                n_min = 1
             else
                n_min = -kmaxn
             end if
             do n = n_min, kmaxn
                kx1 = real(l,dp)*hmatik(1)+real(m,dp)*hmatik(2)+real(n,dp)*hmatik(3)
                ky1 = real(l,dp)*hmatik(4)+real(m,dp)*hmatik(5)+real(n,dp)*hmatik(6)
                kz1 = real(l,dp)*hmatik(7)+real(m,dp)*hmatik(8)+real(n,dp)*hmatik(9)
                ksqr = kx1*kx1+ky1*ky1+kz1*kz1
                ! if ( ksqr .lt. hmaxsq ) then 
                ! sometimes these are about equal, which can cause different
                ! behavior on 32 and 64 bit machines without this .and. statement
                if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq) .gt. 1E-9_dp ) then
                   ncount = ncount + 1
                   my_kx(ncount) = kx1
                   my_ky(ncount) = ky1
                   my_kz(ncount) = kz1
                   my_prefact(ncount) = exp(-ksqr/alpsqr4)/(ksqr*vol)
                   ! sum up q*cos and q*sin ***
                   sumr = 0.0E0_dp
                   sumi = 0.0E0_dp
                   ! do i = myid+1,nchain,numprocs
                   do i = 1,nchain
                      imolty = moltyp(i)
                      if (.not.lelect(imolty).or.nboxi(i).ne.ibox ) cycle
                      do ii = 1,nunit(imolty)
                         if ( lqchg(ntype(imolty,ii)) ) then
                            arg=kx1*rxu(i,ii)+ky1*ryu(i,ii)+kz1*rzu(i ,ii)
                            sumr = sumr + cos(arg)*qqu(i,ii)
                            sumi = sumi + sin(arg)*qqu(i,ii)
                         end if
                      end do
                   end do

                   my_ssumr(ncount) = sumr
                   my_ssumi(ncount) = sumi
                   ! Potential energy ***
                   vrecip = vrecip + (sumr*sumr + sumi*sumi) * my_prefact(ncount)
                end if
             end do
          end do
       end do

       call mp_sum(vrecip,1,groupid)
       call mp_allgather(ncount,rcounts,groupid)
       call mp_set_displs(rcounts,displs,numvect(ibox),numprocs)
       if ( numvect(ibox) .gt. vectormax ) call err_exit(__FILE__,__LINE__,'choose a larger vectormax',myid+1)
       call mp_allgather(my_kx,kx(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_ky,ky(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_kz,kz(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_ssumr,ssumr(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_ssumi,ssumi(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_prefact,prefact(:,ibox),rcounts,displs,groupid)

       return
   end if

  end subroutine recipsum

!> \brief calculates the reciprocal ewald-sum term for trans, rot, flucq,
!> swatch and swap moves, and update the reciprocal ewald-sum.
!> \par History
!> rewritten on June 25/99 by Bin Chen
!> Bai Xue 2016: Enable efficient calculation of cos and sin if lrecip=.true. 
  subroutine recip(ibox,vrecipnew,vrecipold,type)
    integer::ic,izz,ii,imolty,ibox,ncount,ncount1,type,count1,count2
    integer::l,m,n,i,nm,nn
    real::vrecipnew,vrecipold,sumr(2),sumi(2),arg,kx1,ky1,kz1,m_pre,n_pre,bx1,by1,bz1,hmatik(9)
    real::A(vectormax1,vectormax1,vectormax1),B(vectormax1,vectormax1,vectormax1),sumr1(vectormax,2),sumi1(vectormax,2)

    ! RP added for MPI
    integer::rcounts(numprocs),displs(numprocs),my_start,my_end,blocksize
    real::my_ssumrn(vectormax),my_ssumin(vectormax)

! if (LSOLPAR.and.(ibox.eq.2))then
! return
! end if
    ncount = numvect(ibox)
    if ( type .eq. 1 ) then
! recalculate the reciprocal space part for one-particle move, translation,
! rotation, swap, flucq, and swatch.
! old conformation izz = 1 (which is 0 for swap inserted molecule)
! new conformation izz = 2 (which is 0 for swap removed molecule)

#ifdef __DEBUG_KSPACE__
       write(io_output,*) myid,' in recip:',moltion(1),moltion(2)
       do izz = 1,2
          imolty = moltion(izz)
          do ii = 1, nunit(imolty)
             write(io_output,*) rxuion(ii,izz),ryuion(ii,izz),rzuion(ii,izz), qquion(ii,izz)
          end do
       end do
#endif

    if(lrecip .and. (.not. lsolid(ibox)) .or. lrect(ibox)) then

       if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
          bx1 = boxlx(ibox)
          by1 = boxly(ibox)
          bz1 = boxlz(ibox)
          do i = 1,9
             hmatik(i) = 0.0E0_dp
          end do
          hmatik(1) = twopi/bx1
          hmatik(5) = twopi/by1
          hmatik(9) = twopi/bz1
       else
          do i = 1,9
             hmatik(i) = twopi*hmati(ibox,i)
          end do
       end if 
                                    
       do ic=1,ncount
          do izz=1,2
             sumr1(ic,izz) = 0.0E0_dp
             sumi1(ic,izz) = 0.0E0_dp
          end do
       end do
       do izz = 1,2
          imolty = moltion(izz)   
          do ii = 1, nunit(imolty)
             if ( lqchg(ntype(imolty,ii)) .and. abs(qquion(ii,izz)).gt.1e-6 ) then
      
                kx1 = hmatik(1)*rxuion(ii,izz)+hmatik(4)*ryuion(ii,izz)+hmatik(7)*rzuion(ii,izz)
                ky1 = hmatik(2)*rxuion(ii,izz)+hmatik(5)*ryuion(ii,izz)+hmatik(8)*rzuion(ii,izz)
                kz1 = hmatik(3)*rxuion(ii,izz)+hmatik(6)*ryuion(ii,izz)+hmatik(9)*rzuion(ii,izz)
                A(1,1,1) = 1.0E0_dp
                B(1,1,1) = 0.0E0_dp
                A(2,1,1) = cos(kx1)
                B(2,1,1) = sin(kx1)
                A(1,2,1) = cos(ky1)
                B(1,2,1) = sin(ky1)
                A(1,1,2) = cos(kz1)
                B(1,1,2) = sin(kz1)         
                ncount1 = 0
                count2 = 0
                do l = 0,nN0(ibox)
                   ! For MPI purpose. Since different l has different workload, try to give each processor equal workload to some extent.
                   if(mod(l,numprocs).ne.mod(numprocs-mod(count2+myid,numprocs),numprocs)) cycle
                   if(int(l/numprocs).ne.count2) cycle
                   count2 = count2 + 1

                   if(l .gt. 0) then
                      ncount1 = nnc(l,ibox)
                      ncount1 = ncount1 + 1
                      if(numprocs .eq. 1) then
                         A(l+1,1,1) = A(l,1,1)*A(2,1,1)-B(l,1,1)*B(2,1,1)
                         B(l+1,1,1) = B(l,1,1)*A(2,1,1)+A(l,1,1)*B(2,1,1)
                      else
                         A(l+1,1,1) = cos(real(l,dp)*kx1)
                         B(l+1,1,1) = sin(real(l,dp)*kx1)
                      end if
                      sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,1,1)*qquion(ii,izz)
                      sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,1,1)*qquion(ii,izz)
                   end if

                   do m = 0, nNl(l+1,1,ibox)

                      if( m .gt. 0) then
                         ncount1 = ncount1 + 1
                         A(l+1,m+1,1) = A(l+1,m,1)*A(1,2,1)-B(l+1,m,1)*B(1,2,1)
                         B(l+1,m+1,1) = B(l+1,m,1)*A(1,2,1)+A(l+1,m,1)*B(1,2,1)
                         sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,m+1,1)*qquion(ii,izz)
                         sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,m+1,1)*qquion(ii,izz)
                      end if

                      do n = 1, nNlm(l+1,m+1,1,ibox)

                         ncount1 = ncount1 + 1
                         A(l+1,m+1,n+1) = A(l+1,m+1,n)*A(1,1,2)-B(l+1,m+1,n)*B(1,1,2)
                         B(l+1,m+1,n+1) = B(l+1,m+1,n)*A(1,1,2)+A(l+1,m+1,n)*B(1,1,2)
                         sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,m+1,n+1)*qquion(ii,izz)
                         sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,m+1,n+1)*qquion(ii,izz)

                      end do

                      if (l.ne.0.or.m.ne.0) then

                         do n = 1,nNlm(l+1,m+1,2,ibox)
                            ncount1 = ncount1 + 1
                            nn = nNlm(l+1,m+1,2,ibox)
                            if(n .eq. 1) then
                               n_pre = 0
                            else
                               n_pre = nn
                            end if
                            A(l+1,m+1,nn+n+1) = A(l+1,m+1,n_pre+n)*A(1,1,2)+B(l+1,m+1,n_pre+n)*B(1,1,2)
                            B(l+1,m+1,nn+n+1) = B(l+1,m+1,n_pre+n)*A(1,1,2)-A(l+1,m+1,n_pre+n)*B(1,1,2)
                            sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,m+1,nn+n+1)*qquion(ii,izz)
                            sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,m+1,nn+n+1)*qquion(ii,izz)
                         end do

                       end if
                    end do

                    if (l.ne.0) then
                       do m = 1,nNl(l+1,2,ibox)
                          ncount1 = ncount1 + 1
                          nm = nNl(l+1,2,ibox)
                          if(m .eq. 1) then
                             m_pre = 0
                          else
                             m_pre = nm
                          end if
                          A(l+1,nm+m+1,1) = A(l+1,m_pre+m,1)*A(1,2,1)+B(l+1,m_pre+m,1)*B(1,2,1)
                          B(l+1,nm+m+1,1) = B(l+1,m_pre+m,1)*A(1,2,1)-A(l+1,m_pre+m,1)*B(1,2,1)
                          sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,nm+m+1,1)*qquion(ii,izz)
                          sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,nm+m+1,1)*qquion(ii,izz)

                          do n = 1,nNlm(l+1,m+1,3,ibox)
                             ncount1 = ncount1 + 1
                             A(l+1,nm+m+1,n+1) = A(l+1,nm+m+1,n)*A(1,1,2)-B(l+1,nm+m+1,n)*B(1,1,2)
                             B(l+1,nm+m+1,n+1) = B(l+1,nm+m+1,n)*A(1,1,2)+A(l+1,nm+m+1,n)*B(1,1,2)
                             sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,nm+m+1,n+1)*qquion(ii,izz)
                             sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,nm+m+1,n+1)*qquion(ii,izz)
                          end do

                          do n = 1,nNlm(l+1,m+1,4,ibox)
                             ncount1 = ncount1 + 1
                             nn = nNlm(l+1,m+1,4,ibox)
                             if(n .eq. 1) then
                                n_pre = 0
                             else
                                n_pre = nn
                             end if
                             A(l+1,nm+m+1,nn+n+1) = A(l+1,nm+m+1,n_pre+n)*A(1,1,2)+B(l+1,nm+m+1,n_pre+n)*B(1,1,2)
                             B(l+1,nm+m+1,nn+n+1) = B(l+1,nm+m+1,n_pre+n)*A(1,1,2)-A(l+1,nm+m+1,n_pre+n)*B(1,1,2)
                             sumr1(ncount1,izz) = sumr1(ncount1,izz) + A(l+1,nm+m+1,nn+n+1)*qquion(ii,izz)
                             sumi1(ncount1,izz) = sumi1(ncount1,izz) + B(l+1,nm+m+1,nn+n+1)*qquion(ii,izz)
                         end do
                      end do
                   end if
                end do
             end if
          end do
       end do

       call mp_sum(sumr1(:,1),ncount,groupid)
       call mp_sum(sumr1(:,2),ncount,groupid)
       call mp_sum(sumi1(:,1),ncount,groupid)
       call mp_sum(sumi1(:,2),ncount,groupid)
!----------------------------------------------------------------------
       vrecipnew = 0.0E0_dp
       vrecipold = 0.0E0_dp
       blocksize = ncount/numprocs
       rcounts = blocksize
       blocksize = ncount - blocksize * numprocs
       if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
       call mp_set_displs(rcounts,displs,blocksize,numprocs)
       my_start = displs(myid+1) + 1
       my_end = my_start + rcounts(myid+1) - 1
       do ic = my_start,my_end
          count1 = ic - my_start + 1
          my_ssumrn(count1) = ssumr(ic,ibox) - sumr1(ic,1) + sumr1(ic,2)
          my_ssumin(count1) = ssumi(ic,ibox) - sumi1(ic,1) + sumi1(ic,2)
          vrecipnew = vrecipnew + (my_ssumrn(count1)*my_ssumrn(count1) + my_ssumin(count1)*my_ssumin(count1))*prefact(ic,ibox)
          vrecipold = vrecipold + (ssumr(ic,ibox)*ssumr(ic,ibox) + ssumi(ic,ibox)*ssumi(ic,ibox))*prefact(ic,ibox)
       end do
       call mp_allgather(my_ssumrn,ssumrn(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_ssumin,ssumin(:,ibox),rcounts,displs,groupid)

    else

       ! RP added for MPI
       blocksize = ncount/numprocs
       rcounts = blocksize
       blocksize = ncount - blocksize * numprocs
       if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
       call mp_set_displs(rcounts,displs,blocksize,numprocs)
       my_start = displs(myid+1) + 1
       my_end = my_start + rcounts(myid+1) - 1

       ! do 30 ic = 1,ncount
       do ic = my_start,my_end
          do izz = 1,2
             ! izz = 1: old configuration
             ! izz = 2: new configuration
             sumr(izz) = 0.0E0_dp
             sumi(izz) = 0.0E0_dp
             imolty = moltion(izz)
             do ii = 1, nunit(imolty)
                if ( lqchg(ntype(imolty,ii)) ) then
                   arg = kx(ic,ibox)*rxuion(ii,izz) + ky(ic,ibox)*ryuion(ii,izz) + kz(ic,ibox)*rzuion(ii,izz)
                   sumr(izz) = sumr(izz) +  qquion(ii,izz)*cos(arg)
                   sumi(izz) = sumi(izz) +  qquion(ii,izz)*sin(arg)
                end if
             end do
          end do

          ! ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1) + sumr(2)
          ! ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1) + sumi(2)
          ! RP added for MPI
          my_ssumrn(ic-my_start + 1) = ssumr(ic,ibox) - sumr(1) + sumr(2)
          my_ssumin(ic-my_start + 1) = ssumi(ic,ibox) - sumi(1) + sumi(2)
       end do

       call mp_allgather(my_ssumrn,ssumrn(:,ibox),rcounts,displs,groupid)
       call mp_allgather(my_ssumin,ssumin(:,ibox),rcounts,displs,groupid)
!----------------------------------------------------------------------
       vrecipnew = 0.0E0_dp
       vrecipold = 0.0E0_dp
       do ic = my_start,my_end
          vrecipnew = vrecipnew + (ssumrn(ic,ibox)*ssumrn(ic,ibox) + ssumin(ic,ibox)*ssumin(ic,ibox))*prefact(ic,ibox)
          vrecipold = vrecipold + (ssumr(ic,ibox)*ssumr(ic,ibox) + ssumi(ic,ibox)*ssumi(ic,ibox))*prefact(ic,ibox)
       end do

    end if

    call mp_sum(vrecipnew,1,groupid)
    call mp_sum(vrecipold,1,groupid)
    vrecipnew = vrecipnew*qqfact
    vrecipold = vrecipold*qqfact

    else if (type .eq. 2) then
       ! update the reciprocal space k vectors
       do ic = 1, ncount
          ssumr(ic,ibox) = ssumrn(ic,ibox)
          ssumi(ic,ibox) = ssumin(ic,ibox)
       end do
    else if (type .eq. 3) then
       ! store the reciprocal space k vectors
       do ic = 1, ncount
          ssumro(ic,ibox) = ssumr(ic,ibox)
          ssumio(ic,ibox) = ssumi(ic,ibox)
       end do
       if(lrecip .and. (.not. lsolid(ibox)) .or. lrect(ibox)) then
          nN0o(ibox) = nN0(ibox)
          nNlo(:,:,ibox) = nNl(:,:,ibox)
          nNlmo(:,:,:,ibox) = nNlm(:,:,:,ibox)
          k_max_lo(ibox) = k_max_l(ibox)
          k_max_mo(ibox) = k_max_m(ibox)
          k_max_no(ibox) = k_max_n(ibox)
          nnco(:,ibox) = nnc(:,ibox)
       end if
    else if (type .eq. 4) then
       ! restore the reciprocal space k vectors
       do ic = 1, ncount
          ssumr(ic,ibox) = ssumro(ic,ibox)
          ssumi(ic,ibox) = ssumio(ic,ibox)
       end do
       if(lrecip .and. (.not. lsolid(ibox)) .or. lrect(ibox)) then
          nN0(ibox) = nN0o(ibox)
          nNl(:,:,ibox) = nNlo(:,:,ibox)
          nNlm(:,:,:,ibox) = nNlmo(:,:,:,ibox)
          k_max_l(ibox) = k_max_lo(ibox)
          k_max_m(ibox) = k_max_mo(ibox)
          k_max_n(ibox) = k_max_no(ibox)
          nnc(:,ibox) = nnco(:,ibox)
       end if
    end if

#ifdef __DEBUG_KSPACE__
    write(io_output,*) myid,' in recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)
#endif
    return
  end subroutine recip

!> \brief store old k vectors and reciprocal sum
  subroutine save_kvector(ibox)
    integer,intent(in)::ibox
    integer::ic,ncount
    calpo(ibox) = calp(ibox)
    numvecto(ibox) = numvect(ibox)
    ncount = numvect(ibox)
    do ic = 1,ncount
       kxo(ic,ibox) = kx(ic,ibox)
       kyo(ic,ibox) = ky(ic,ibox)
       kzo(ic,ibox) = kz(ic,ibox)
       prefacto(ic,ibox) = prefact(ic,ibox)
    end do
  end subroutine save_kvector

!> \brief restore old k vectors and reciprocal sum and calp
  subroutine restore_kvector(ibox)
    integer,intent(in)::ibox
    integer::ic,ncount

    calp(ibox) = calpo(ibox)
    numvect(ibox) = numvecto(ibox)
    ncount = numvecto(ibox)
    do ic = 1,ncount
       kx(ic,ibox) = kxo(ic,ibox)
       ky(ic,ibox) = kyo(ic,ibox)
       kz(ic,ibox) = kzo(ic,ibox)
       prefact(ic,ibox) = prefacto(ic,ibox)
    end do
  end subroutine restore_kvector

  subroutine recip_atom(ibox,vrecipnew,vrecipold,type,ii)
      integer::ic,izz,ii,imolty,ibox,ncount,type
      real::vrecipnew,vrecipold,sumr(2),sumi(2) ,arg

      ncount = numvect(ibox)

      if ( type .eq. 1 ) then
! recalculate the reciprocal space part for one-particle move, translation,
! rotation, swap, flucq, and swatch.
! old conformation izz = 1 (which is 0 for swap inserted molecule)
! new conformation izz = 2 (which is 0 for swap removed molecule)

#ifdef __DEBUG_KSPACE__
         write(io_output,*) myid,' in recip_atom:',moltion(1),moltion(2)
         do izz = 1,2
            imolty = moltion(izz)
            do ii = 1, nunit(imolty)
               write(io_output,*) rxuion(ii,izz),ryuion(ii,izz),rzuion(ii,izz),qquion(ii,izz)
            end do
         end do
#endif

         do 30 ic = 1, ncount
            do 20 izz = 1,2
! izz = 1: old configuration
! izz = 2: new configuration

               sumr(izz) = 0.0E0_dp
               sumi(izz) = 0.0E0_dp
               imolty = moltion(izz)
                  if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,izz) + ky(ic,ibox)*ryuion(ii,izz) + kz(ic,ibox)*rzuion(ii,izz)
                     sumr(izz) = sumr(izz) +  qquion(ii,izz)*cos(arg)
                     sumi(izz) = sumi(izz) +  qquion(ii,izz)*sin(arg)
                  end if
 20         continue
            ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1) + sumr(2)
            ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1) + sumi(2)
 30      continue
         vrecipnew = 0.0E0_dp
         vrecipold = 0.0E0_dp
         do ic = 1,ncount
            vrecipnew = vrecipnew + (ssumrn(ic,ibox)* ssumrn(ic,ibox) + ssumin(ic,ibox)* ssumin(ic,ibox))*prefact(ic,ibox)
            vrecipold = vrecipold + (ssumr(ic,ibox)* ssumr(ic,ibox) + ssumi(ic,ibox)* ssumi(ic,ibox))*prefact(ic,ibox)
         end do

         vrecipnew = vrecipnew*qqfact
         vrecipold = vrecipold*qqfact

      else if (type .eq. 2) then

! update the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumrn(ic,ibox)
            ssumi(ic,ibox) = ssumin(ic,ibox)
         end do

      else if (type .eq. 3) then

! store the reciprocal space k vectors

         do ic = 1, ncount
            ssumro(ic,ibox) = ssumr(ic,ibox)
            ssumio(ic,ibox) = ssumi(ic,ibox)
         end do

      else if (type .eq. 4) then

! restore the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumro(ic,ibox)
            ssumi(ic,ibox) = ssumio(ic,ibox)
         end do

      end if

#ifdef __DEBUG_KSPACE__
      write(io_output,*) myid,' in recip_atom:',ssumr(100,ibox),ibox,ssumrn(100,ibox)
#endif

      return
  end subroutine recip_atom

  subroutine ee_recip(ibox,vrecipnew,vrecipold,type)
      integer::ic,zzz,ii,imolty,ibox,ncount,type
      real::vrecipnew,vrecipold,sumr(2),sumi(2) ,arg

      ncount = numvect(ibox)

      if ( type .eq. 1 ) then

! recalculate the reciprocal space part for one-particle move, translation,
! rotation, swap, flucq, and swatch.
! old conformation zzz = 1 (which is 0 for swap inserted molecule)
! new conformation zzz = 2 (which is 0 for swap removed molecule)

#ifdef __DEBUG_KSPACE__
         write(io_output,*) myid,' in ee_recip:',moltion(1),moltion(2)
         do zzz = 1,2
            imolty = moltion(zzz)
            do ii = 1, nunit(imolty)
               write(io_output,*) rxuion(ii,zzz),ryuion(ii,zzz),rzuion(ii,zzz),qquion(ii,zzz)
            end do
         end do
#endif

         do 30 ic = 1, ncount
            do 20 zzz = 1,2
! zzz = 1: old configuration
! zzz = 2: new configuration

               sumr(zzz) = 0.0E0_dp
               sumi(zzz) = 0.0E0_dp
               imolty = moltion(zzz)
               do ii = 1, nunit(imolty)
! if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,zzz) + ky(ic,ibox)*ryuion(ii,zzz) + kz(ic,ibox)*rzuion(ii,zzz)
                     sumr(zzz) = sumr(zzz) +  qquion(ii,zzz)*cos(arg)
                     sumi(zzz) = sumi(zzz) +  qquion(ii,zzz)*sin(arg)
! end if
               end do
 20         continue
            ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1) + sumr(2)
            ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1) + sumi(2)
 30      continue
         vrecipnew = 0.0E0_dp
         vrecipold = 0.0E0_dp
         do ic = 1,ncount
            vrecipnew = vrecipnew + (ssumrn(ic,ibox)* ssumrn(ic,ibox) + ssumin(ic,ibox)* ssumin(ic,ibox))*prefact(ic,ibox)
            vrecipold = vrecipold + (ssumr(ic,ibox)* ssumr(ic,ibox) + ssumi(ic,ibox)* ssumi(ic,ibox))*prefact(ic,ibox)
         end do

         vrecipnew = vrecipnew*qqfact
         vrecipold = vrecipold*qqfact

      else if (type .eq. 2) then

! update the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumrn(ic,ibox)
            ssumi(ic,ibox) = ssumin(ic,ibox)
         end do

      else if (type .eq. 3) then

! store the reciprocal space k vectors

         do ic = 1, ncount
            ssumro(ic,ibox) = ssumr(ic,ibox)
            ssumio(ic,ibox) = ssumi(ic,ibox)
         end do

      else if (type .eq. 4) then

! restore the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumro(ic,ibox)
            ssumi(ic,ibox) = ssumio(ic,ibox)
         end do

      end if

#ifdef __DEBUG_KSPACE__
      write(io_output,*) myid,' in ee_recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)
#endif
      return
  end subroutine ee_recip

!> \brief Calculates the reciprocal space contribution to pressure using
!> thermodynamic definition.
!> \see J. Chem. Phys. Vol. 109 P2791.
!> \par History
!> written in 1998 by Bin Chen \n
!> modified to calculate surface tension, 11/24/03 JMS
!> Bai Xue 2016: Enable efficient calculation of cos and sin if lrecip=.true. 
  subroutine recippress(ibox,repress,pxx,pyy,pzz,pxy,pyx,pxz,pzx, pyz,pzy)
    integer::ncount,ibox,i,ii,imolty,l,m,n,m_pre,n_pre,nm,nn
    real::factor,repress,repressx,repressy,repressz,recipintra,piix,piiy,piiz,xcmi,ycmi,zcmi,arg
    real::pxx,pyy,pzz,intraxx,intrayy,intrazz,intraxy,intraxz,intrazy,intrayz,intrayx,intrazx,pxy,pyx,pyz,pzy,pxz,pzx
    real::kx1,ky1,kz1,A(vectormax1,vectormax1,vectormax1),B(vectormax1,vectormax1,vectormax1),hmatik(9),bx1,by1,bz1

    repress  = 0.0E0_dp
    repressx = 0.0E0_dp
    repressy = 0.0E0_dp
    repressz = 0.0E0_dp
    recipintra = 0.0E0_dp
    pxy = 0.0E0_dp
    pxz = 0.0E0_dp
    pyx = 0.0E0_dp
    pyz = 0.0E0_dp
    pzx = 0.0E0_dp
    pzy = 0.0E0_dp

    intraxx = 0.0E0_dp
    intrayy = 0.0E0_dp
    intrazz = 0.0E0_dp
    intraxy = 0.0E0_dp
    intrazy = 0.0E0_dp
    intraxz = 0.0E0_dp
    intrazx = 0.0E0_dp
    intrayz = 0.0E0_dp
    intrayx = 0.0E0_dp

    if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
       bx1 = boxlx(ibox)
       by1 = boxly(ibox)
       bz1 = boxlz(ibox)
       do i = 1,9
          hmatik(i) = 0.0E0_dp
       end do
       hmatik(1) = twopi/bx1
       hmatik(5) = twopi/by1
       hmatik(9) = twopi/bz1
    else
       do i = 1,9
          hmatik(i) = twopi*hmati(ibox,i)
       end do
    end if


    ! RP for MPI
    do ncount = myid+1,numvect(ibox),numprocs
       ! do ncount = 1, numvect(ibox)
       factor = prefact(ncount,ibox)*(ssumr(ncount,ibox)*ssumr(ncount,ibox) + ssumi(ncount,ibox)* ssumi(ncount,ibox))
       repressx = repressx + factor*(1.0E0_dp - (1.0E0_dp/(4.0E0_dp*calp(ibox) *calp(ibox))&
        + 1.0E0_dp/(kx(ncount,ibox)*kx(ncount,ibox)+ ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)* kz(ncount,ibox)))&
        *2.0E0_dp*kx(ncount,ibox)*kx(ncount,ibox))
       repressy = repressy + factor*(1.0E0_dp - (1.0E0_dp/(4.0E0_dp*calp(ibox) *calp(ibox))&
        + 1.0E0_dp/(kx(ncount,ibox)*kx(ncount,ibox)+ ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)* kz(ncount,ibox)))&
        *2.0E0_dp*ky(ncount,ibox)*ky(ncount,ibox))
       repressz = repressz + factor*(1.0E0_dp - (1.0E0_dp/(4.0E0_dp*calp(ibox) *calp(ibox))&
        + 1.0E0_dp/(kx(ncount,ibox)*kx(ncount,ibox)+ ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)* kz(ncount,ibox)))&
        *2.0E0_dp*kz(ncount,ibox)*kz(ncount,ibox))
       pxy = pxy + factor*(0.0E0_dp - (1.0E0_dp/(4.0E0_dp*calp(ibox) *calp(ibox))&
        + 1.0E0_dp/(kx(ncount,ibox)*kx(ncount,ibox)+ ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)* kz(ncount,ibox)))&
        *2.0E0_dp*kx(ncount,ibox)*ky(ncount,ibox))
       pxz = pxz + factor*(0.0E0_dp - (1.0E0_dp/(4.0E0_dp*calp(ibox) *calp(ibox))&
        + 1.0E0_dp/(kx(ncount,ibox)*kx(ncount,ibox)+ ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)* kz(ncount,ibox)))&
        *2.0E0_dp*kx(ncount,ibox)*kz(ncount,ibox))
       pyz = pyz + factor*(0.0E0_dp - (1.0E0_dp/(4.0E0_dp*calp(ibox) *calp(ibox))&
        + 1.0E0_dp/(kx(ncount,ibox)*kx(ncount,ibox)+ ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)* kz(ncount,ibox)))&
        *2.0E0_dp*ky(ncount,ibox)*kz(ncount,ibox))
    end do

    ! RP added for MPI
    call mp_sum(repressx,1,groupid)
    call mp_sum(repressy,1,groupid)
    call mp_sum(repressz,1,groupid)
    call mp_sum(pxy,1,groupid)
    call mp_sum(pxz,1,groupid)
    call mp_sum(pyz,1,groupid)

    repress = repressx + repressy + repressz
    ! keep x,y,z separate for surface tension calculation
    pxx = repressx
    pyy = repressy
    pzz = repressz
    pyx = pxy
    pzx = pxz
    pzy = pyz

    if(lrecip .and. (.not. lsolid(ibox)) .or. lrect(ibox)) then

       do i = myid+1, nchain, numprocs
          ! check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             if ( .not. lelect(imolty) ) cycle
             xcmi = xcm(i)
             ycmi = ycm(i)
             zcmi = zcm(i)

             ! loop over all beads ii of chain i
             do ii = 1, nunit(imolty)
                if(abs(qqu(i,ii)).lt.1e-6) cycle
                ! compute the vector of the bead to the COM (p)
                piix = rxu(i,ii) - xcmi
                piiy = ryu(i,ii) - ycmi
                piiz = rzu(i,ii) - zcmi

                ncount = 0
                kx1 = hmatik(1)*rxu(i,ii)+hmatik(4)*ryu(i,ii)+hmatik(7)*rzu(i,ii)
                ky1 = hmatik(2)*rxu(i,ii)+hmatik(5)*ryu(i,ii)+hmatik(8)*rzu(i,ii)
                kz1 = hmatik(3)*rxu(i,ii)+hmatik(6)*ryu(i,ii)+hmatik(9)*rzu(i,ii)
                A(1,1,1) = 1.0E0_dp
                B(1,1,1) = 0.0E0_dp
                A(2,1,1) = cos(kx1)
                B(2,1,1) = sin(kx1)
                A(1,2,1) = cos(ky1)
                B(1,2,1) = sin(ky1)
                A(1,1,2) = cos(kz1)
                B(1,1,2) = sin(kz1)
                do l = 0,nN0(ibox)

                   if( l .gt. 0) then
                      ncount = ncount + 1
                      A(l+1,1,1) = A(l,1,1)*A(2,1,1)-B(l,1,1)*B(2,1,1)
                      B(l+1,1,1) = B(l,1,1)*A(2,1,1)+A(l,1,1)*B(2,1,1)
                      factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,1,1)+ssumi(ncount,ibox)*A(l+1,1,1))*qqu(i,ii)
                      recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                      ! keep x,y and z separate for surface tension calculation
                      intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                      intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                      intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                      intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                      intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                      intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                      intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                      intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                      intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)
                   end if

                   do m = 0, nNl(l+1,1,ibox)

                      if( m .gt. 0) then
                         ncount = ncount + 1
                         A(l+1,m+1,1) = A(l+1,m,1)*A(1,2,1)-B(l+1,m,1)*B(1,2,1)
                         B(l+1,m+1,1) = B(l+1,m,1)*A(1,2,1)+A(l+1,m,1)*B(1,2,1)
                         factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,m+1,1)+ssumi(ncount,ibox)*A(l+1,m+1,1))*qqu(i,ii)
                         recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                         intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                         intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                         intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                         intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                         intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                         intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                         intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                         intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                         intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)
                      end if

                      do n = 1, nNlm(l+1,m+1,1,ibox)
                         ncount = ncount + 1
                         A(l+1,m+1,n+1) = A(l+1,m+1,n)*A(1,1,2)-B(l+1,m+1,n)*B(1,1,2)
                         B(l+1,m+1,n+1) = B(l+1,m+1,n)*A(1,1,2)+A(l+1,m+1,n)*B(1,1,2)
                         factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,m+1,n+1)+ssumi(ncount,ibox)*A(l+1,m+1,n+1))*qqu(i,ii)
                         recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                         intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                         intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                         intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                         intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                         intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                         intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                         intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                         intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                         intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)
                      end do

                      if (l.ne.0.or.m.ne.0) then

                         do n = 1,nNlm(l+1,m+1,2,ibox)
                            ncount = ncount + 1
                            nn = nNlm(l+1,m+1,2,ibox)
                            if(n .eq. 1) then
                               n_pre = 0
                            else
                               n_pre = nn
                            end if
                            A(l+1,m+1,nn+n+1) = A(l+1,m+1,n_pre+n)*A(1,1,2)+B(l+1,m+1,n_pre+n)*B(1,1,2)
                            B(l+1,m+1,nn+n+1) = B(l+1,m+1,n_pre+n)*A(1,1,2)-A(l+1,m+1,n_pre+n)*B(1,1,2)
                            factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,m+1,nn+n+1)+ssumi(ncount,ibox)*A(l+1,m+1,nn+n+1))*qqu(i,ii)
                            recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                            intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                            intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                            intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                            intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                            intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                            intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                            intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                            intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                            intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)
                         end do
                      end if
                   end do

                   if (l .ne. 0) then
                      do m = 1,nNl(l+1,2,ibox)
                         ncount = ncount + 1
                         nm = nNl(l+1,2,ibox)
                         if(m .eq. 1) then
                            m_pre = 0
                         else
                            m_pre = nm
                         end if
                         A(l+1,nm+m+1,1) = A(l+1,m_pre+m,1)*A(1,2,1)+B(l+1,m_pre+m,1)*B(1,2,1)
                         B(l+1,nm+m+1,1) = B(l+1,m_pre+m,1)*A(1,2,1)-A(l+1,m_pre+m,1)*B(1,2,1)
                         factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,nm+m+1,1)+ssumi(ncount,ibox)*A(l+1,nm+m+1,1))*qqu(i,ii)
                         recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                         intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                         intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                         intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                         intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                         intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                         intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                         intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                         intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                         intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)

                         do n = 1, nNlm(l+1,m+1,3,ibox)
                            ncount = ncount + 1
                            A(l+1,nm+m+1,n+1) = A(l+1,nm+m+1,n)*A(1,1,2)-B(l+1,nm+m+1,n)*B(1,1,2)
                            B(l+1,nm+m+1,n+1) = B(l+1,nm+m+1,n)*A(1,1,2)+A(l+1,nm+m+1,n)*B(1,1,2)
                            factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,nm+m+1,n+1)+ssumi(ncount,ibox)*A(l+1,nm+m+1,n+1))*qqu(i,ii)
                            recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                            intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                            intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                            intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                            intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                            intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                            intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                            intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                            intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                            intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)
                         end do

                         do n = 1,nNlm(l+1,m+1,4,ibox)
                            ncount = ncount + 1
                            nn = nNlm(l+1,m+1,4,ibox)
                            if(n .eq. 1) then
                               n_pre = 0
                            else
                               n_pre = nn
                            end if
                            A(l+1,nm+m+1,nn+n+1) = A(l+1,nm+m+1,n_pre+n)*A(1,1,2)+B(l+1,nm+m+1,n_pre+n)*B(1,1,2)
                            B(l+1,nm+m+1,nn+n+1) = B(l+1,nm+m+1,n_pre+n)*A(1,1,2)-A(l+1,nm+m+1,n_pre+n)*B(1,1,2)
                            factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*B(l+1,nm+m+1,nn+n+1)+ssumi(ncount,ibox)*A(l+1,nm+m+1,nn+n+1))*qqu(i,ii)
                            recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                            intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                            intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                            intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                            intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                            intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                            intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                            intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                            intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                            intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)
                         end do
                      end do
                   end if
                end do
             end do
          end if
       end do

    else

       ! the intramolecular part should be substracted
       ! RP for MPI
       do i = myid+1, nchain, numprocs
          ! check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             if ( .not. lelect(imolty) ) cycle
             xcmi = xcm(i)
             ycmi = ycm(i)
             zcmi = zcm(i)

             ! loop over all beads ii of chain i
             do ii = 1, nunit(imolty)
                ! compute the vector of the bead to the COM (p)
                piix = rxu(i,ii) - xcmi
                piiy = ryu(i,ii) - ycmi
                piiz = rzu(i,ii) - zcmi

                do ncount = 1,numvect(ibox)
                   ! compute the dot product of k and r
                   arg = kx(ncount,ibox)*rxu(i,ii) + ky(ncount,ibox)*ryu(i,ii) + kz(ncount,ibox)*rzu(i,ii)
                   factor = prefact(ncount,ibox)*2.0E0_dp*(-ssumr(ncount,ibox)*sin(arg)+ssumi(ncount,ibox)*cos(arg))*qqu(i,ii)
                   recipintra = recipintra + factor*(kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy+kz(ncount,ibox)*piiz)
                   ! keep x,y and z separate for surface tension calculation
                   intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                   intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                   intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                   intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                   intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                   intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                   intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                   intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                   intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)

                end do
             end do
          end if
       end do

    end if

    call mp_sum(recipintra,1,groupid)
    call mp_sum(intraxx,1,groupid)
    call mp_sum(intrayy,1,groupid)
    call mp_sum(intrazz,1,groupid)
    call mp_sum(intraxy,1,groupid)
    call mp_sum(intraxz,1,groupid)
    call mp_sum(intrayx,1,groupid)
    call mp_sum(intrayz,1,groupid)
    call mp_sum(intrazx,1,groupid)
    call mp_sum(intrazy,1,groupid)

    repress = (repress + recipintra)*qqfact

    pxx = (pxx + intraxx)*qqfact
    pyy = (pyy + intrayy)*qqfact
    pzz = (pzz + intrazz)*qqfact

    pxy = pxy + intraxy
    pyx = pyx + intrayx
    pxz = pxz + intraxz
    pzx = pzx + intrazx
    pyz = pyz + intrayz
    pzy = pzy + intrazy

#ifdef __DEBUG_KSPACE__
    write(io_output,*) myid,' in recippress. Internal part:',intraxx,intrayy,intrazz
#endif
    return
  end subroutine recippress

  subroutine allocate_kspace()
    use sim_system,only:nbxmax
    integer::jerr
    if (allocated(kx)) deallocate(kx,ky,kz,prefact,ssumr,ssumi,ssumrn,ssumin,ssumro,ssumio,kxo,kyo,kzo,prefacto,calpo,calp,numvect,numvecto,&
                                  k_max_l,k_max_m,k_max_n,k_max_lo,k_max_mo,k_max_no,nN0,nNl,nNlm,nN0o,nNlo,nNlmo,nnc,nnco,stat=jerr)
    allocate(kx(vectormax,nbxmax),ky(vectormax,nbxmax),kz(vectormax,nbxmax),prefact(vectormax,nbxmax)&
     ,ssumr(vectormax,nbxmax),ssumi(vectormax,nbxmax),ssumrn(vectormax,nbxmax),ssumin(vectormax,nbxmax)&
     ,ssumro(vectormax,nbxmax),ssumio(vectormax,nbxmax),kxo(vectormax,nbxmax),kyo(vectormax,nbxmax)&
     ,kzo(vectormax,nbxmax),prefacto(vectormax,nbxmax),calpo(nbxmax),calp(nbxmax),numvect(nbxmax),numvecto(nbxmax)&
     ,k_max_l(nbxmax),k_max_m(nbxmax),k_max_n(nbxmax),k_max_lo(nbxmax),k_max_mo(nbxmax),k_max_no(nbxmax),nN0(nbxmax),nNl(vectormax1,2,nbxmax)&
     ,nNlm(vectormax1,vectormax1,4,nbxmax),nN0o(nbxmax),nNlo(vectormax1,2,nbxmax)&
     ,nNlmo(vectormax1,vectormax1,4,nbxmax),nnc(vectormax1,nbxmax),nnco(vectormax1,nbxmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'allocate_kspace: allocation failed',jerr)
    numvect=0
  end subroutine allocate_kspace

  subroutine compute_kmax(ibox)
    integer, intent(in)::ibox
    if ((.not.lsolid(ibox)).or.lrect(ibox)) then
       k_max_l(ibox) = aint(calp(ibox)*calp(ibox)*boxlx(ibox)*rcut(ibox)/onepi)+1
       k_max_m(ibox) = aint(calp(ibox)*calp(ibox)*boxly(ibox)*rcut(ibox)/onepi)+1
       k_max_n(ibox) = aint(calp(ibox)*calp(ibox)*boxlz(ibox)*rcut(ibox)/onepi)+1
    else
       k_max_l(ibox) = aint(calp(ibox)*calp(ibox)*hmat(ibox,1)*rcut(ibox)/onepi)+2
       k_max_m(ibox) = aint(calp(ibox)*calp(ibox)*hmat(ibox,5)*rcut(ibox)/onepi)+2
       k_max_n(ibox) = aint(calp(ibox)*calp(ibox)*hmat(ibox,9)*rcut(ibox)/onepi)+2
    end if
  end subroutine compute_kmax
    
end MODULE energy_kspace
