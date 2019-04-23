!> \brief Feuston-Garofalini force field for SiO2/H2O
!>
!> It is composed of a modified Born-Mayer-Huggins (BMH) pair potential
!> and a three-body potential. The BMH part is the sum of an exponential
!> repulsive term and a screened Coulomb term. (J Phys Chem 1990,94:5351)
!> U2(r) = A*exp(-r/rho) + q_i*q_j/r*erfc(r/beta)
!>       + \sigma_{x=1}^{x=6} a^x/{1+exp[b^x*(r-c^x)]}
!> U3(r_jik) = h3(r_ij,r_ik,theta_jik) + h3(r_jk,r_ji,theta_kji) + h3(r_ki,r_kj,theta_ikj)
!> h3(rij,rik,theta_jik) = lambda*exp[gamma_ij/(r_ij-r_ij^*)+gamma_ik/(r_ik-r_ik^*)]
!>                  *[cos(theta_jik)-cos(theta_jik^*)]^2 for r_ij<r_ij^* and r_ik<r_ik^*
!>                       = 0   otherwise
!> vvdW_1 = A, vvdW_2 = rho, vvdW_3 = q_i*q_j, vvdW_4 = beta
!> \todo Combine with energy_3body for a more general implementation.
!> \author KE ANDERSON
MODULE energy_garofalini
  use var_type,only:dp
  use const_phys,only:qqfact
  use util_math,only:erfunc
  use sim_particle,only:neigh_cnt,neighbor,ndij,nxij,nyij,nzij
  use sim_system
  implicit none
  private
  save
  public::init_garofalini,idx_garofalini,isNeighbor,garofalini,vthreebody,triad,triad_en

  integer,parameter::pair_max=50000
  real::ga(6,3),gb(6,3),gc(6,3),glambda(4),grij(4,2),ggamma(4,2),gtheta(4),grijsq(4,2)
  real::dxij(pair_max),dyij(pair_max),dzij(pair_max),dij(pair_max),dik(pair_max),dxik(pair_max),dyik(pair_max),dzik(pair_max)
  integer::itr1(pair_max),itr2(pair_max),itr3(pair_max),ntr

contains
  subroutine init_garofalini()
    use const_math,only:degrad
    nntype = 3
    vvdW(:,1:6) = 0.0_dp
    ga = 0.0_dp
    gb = 0.0_dp
    gc = 0.0_dp
    glambda = 0.0_dp
    ggamma = 0.0_dp
    grij = 0.0_dp
    grijsq = 0.0_dp
    gtheta = 0.0_dp
    lij(1:6) = .true.
    qelect(1:6) = 0.0_dp
    lqchg(1:6) = .false.

    ! Parameters (galpha,grho,gbeta,ga,gb,gc; lambda,grij,ggamma,gtheta)
    ! Si-Si
    vvdW(1,1) = 13597175.7E0_dp
    vvdW(2,1) = 0.29E0_dp
    vvdW(3,1) = 16.0_dp
    vvdW(4,1) = 2.29E0_dp
    mass(1) = 28.09E0_dp
    chemid(1) = 'Si'

    ! O-O
    vvdW(1,2) = 5251972.5E0_dp
    vvdW(2,2) = 0.29E0_dp
    vvdW(3,2) = 4.0_dp
    vvdW(4,2) = 2.34E0_dp
    mass(2) = 16.00E0_dp
    chemid(2) = 'O  '

    ! H-H
    vvdW(1,3) = 246299.4E0_dp
    vvdW(2,3) = 0.35E0_dp
    vvdW(3,3) = 1.0_dp
    vvdW(4,3) = 2.1E0_dp
    ga(3,1) = -38243.8E0_dp
    gb(3,1) = 6.0E0_dp
    gc(3,1) = 1.51E0_dp
    ga(3,2) = 2515.9E0_dp
    gb(3,2) = 2.0E0_dp
    gc(3,2) = 2.42E0_dp
    mass(3) = 1.0078E0_dp
    chemid(3) = 'H'

    ! Si-O
    vvdW(1,4) =  21457024.2E0_dp
    vvdW(2,4) = 0.29E0_dp
    vvdW(3,4) = -8.0_dp
    vvdW(4,4) = 2.34E0_dp

    ! Si-H
    vvdW(1,5) = 499842.9E0_dp
    vvdW(2,5) = 0.29E0_dp
    vvdW(3,5) = 4.0_dp
    vvdW(4,5) = 2.31E0_dp
    ga(5,1) = -33715.5E0_dp
    gb(5,1) = 6.0E0_dp
    gc(5,1) = 2.2E0_dp

    ! 0-H
    vvdW(1,6) = 2886049.4E0_dp
    vvdW(2,6) = 0.29E0_dp
    vvdW(3,6) = -2.0_dp
    vvdW(4,6) = 2.26E0_dp
    ga(6,1) = -15096.7E0_dp
    gb(6,1) = 15.0E0_dp
    gc(6,1) = 1.05E0_dp
    ga(6,2) = 55353.6E0_dp
    gb(6,2) = 3.2E0_dp
    gc(6,2) = 1.50E0_dp
    ga(6,3) = -6038.7E0_dp
    gb(6,3) = 5.0E0_dp
    gc(6,3) = 2.0E0_dp

    ! Si-O-Si
    glambda(1) = 21732.3E0_dp
    ggamma(1,1) = 2.0E0_dp
    grij(1,1) = 2.6E0_dp
    grijsq(1,1) = grij(1,1)*grij(1,1)
    gtheta(1) = cos(109.5E0_dp*degrad)

    ! O-Si-O
    glambda(2) = 1376379.0E0_dp
    ggamma(2,1) = 2.8E0_dp
    grij(2,1) = 3.0E0_dp
    grijsq(2,1) = grij(2,1)*grij(2,1)
    gtheta(2) = cos(109.5E0_dp*degrad)

    ! H-O-H
    glambda(3) = 2535435.0E0_dp
    ggamma(3,1) = 1.3E0_dp
    grij(3,1) = 1.6E0_dp
    grijsq(3,1) = grij(3,1)*grij(3,1)
    gtheta(3) = cos(104.5E0_dp*degrad)

    ! Si-O-H
    glambda(4) = 362205.0E0_dp
    ggamma(4,1) = 2.0E0_dp
    ggamma(4,2) = 1.2E0_dp
    grij(4,1) = grij(1,1)
    grij(4,2) = 1.5E0_dp
    grijsq(4,1) = grij(4,1)*grij(4,1)
    grijsq(4,2) = grij(4,2)*grij(4,2)
    gtheta(4) = cos(109.5E0_dp*degrad)
  end subroutine init_garofalini

!DEC$ ATTRIBUTES FORCEINLINE :: type_2body
  function idx_garofalini(ntii,ntjj)
    integer::idx_garofalini
    integer,intent(in)::ntii,ntjj

    if (ntii.eq.ntjj) then
       idx_garofalini = ntii
    else
       idx_garofalini = ntii+ntjj+1
    end if
  end function idx_garofalini

  function isNeighbor(rijsq,ntii,ntjj)
    logical::isNeighbor
    real,intent(in)::rijsq
    integer,intent(in)::ntii,ntjj
    integer::ntij

    ntij = idx_garofalini(ntii,ntjj)
    isNeighbor = (ntij.eq.4.and.rijsq.lt.grijsq(2,1)).or.(ntij.eq.6.and.rijsq.lt.grijsq(3,1))
  end function isNeighbor

  function garofalini(rij,ntii,ntjj)
    real::garofalini
    real,intent(in)::rij
    integer,intent(in)::ntii,ntjj
    integer::ntij,i

    ntij = idx_garofalini(ntii,ntjj)

    garofalini = vvdW(1,ntij)*exp(-rij/vvdW(2,ntij)) + vvdW(3,ntij)*qqfact/rij*erfunc(rij/vvdW(4,ntij))

    ! CSF term
    do i=1,3
       garofalini = garofalini + ga(ntij,i)/(1+exp(gb(ntij,i)*(rij-gc(ntij,i))))
    end do
  end function garofalini

  function vthreebody() result(vthree)
    real::vthree
    integer::ntang,nta,ntb,tri,i,j,k
    real::vthreea,thetac,g,p

#ifdef __DEBUG_GAROFALINI__
    write(io_output,*) 'start Vthreebody in ',myid
#endif

    vthree = 0.0E0_dp
    do tri = 1,ntr
       i = itr1(tri)
       j = itr2(tri)
       k = itr3(tri)

       ! skip if i (central atom) is hydrogen
       if(ntype(moltyp(i),1).eq.3) cycle

       ! determine type
       if(ntype(moltyp(j),1).eq.ntype(moltyp(k),1)) then
          ntang = ntype(moltyp(j),1)
#ifdef __DEBUG_GAROFALINI__
             write(30,*) tri,'a:',j,i,k,'(',ntang,')',dij(tri),dik(tri)
#endif
          nta = 1
          ntb = 1
          if ((dij(tri).gt.grij(ntang,1)).or.(dik(tri).gt.grij(ntang,1))) then
#ifdef __DEBUG_GAROFALINI__
             write(30,*) '   catch distance Si-O-Si',dij(tri),dik(tri)
#endif
             cycle
          end if
       else
          ntang = 4
#ifdef __DEBUG_GAROFALINI__
             write(30,*) tri,'b:',j,i,k,'(',ntang,')',dij(tri),dik(tri)
#endif
          if(ntype(moltyp(j),1).eq.1) then
             nta = 1
             ntb = 2
             if((dij(tri).gt.grij(4,nta)).or.(dik(tri).gt.grij(4,ntb))) then
#ifdef __DEBUG_GAROFALINI__
                write(30,*)'  catch distance Si-O-H',dij(tri),dik(tri)
#endif
                cycle
             end if
          else
             nta = 2
             ntb = 1
             if((dij(tri).gt.grij(4,nta)).or.(dik(tri).gt.grij(4,ntb))) then
#ifdef __DEBUG_GAROFALINI__
                write(30,*) '  catch distance H-O-Si',dij(tri),dik(tri)
#endif
                cycle
             end if
          end if
       end if

       ! dij = sqrt(dijsq(tri))
       ! dik = sqrt(diksq(tri))

       thetac = (dxij(tri)*dxik(tri)+dyij(tri)*dyik(tri)+ dzij(tri)*dzik(tri))/(dij(tri)*dik(tri))
       if ( thetac .ge. 1.0E0_dp ) thetac = 1.0E0_dp
       if ( thetac .le. -1.0E0_dp ) thetac = -1.0E0_dp

       p = (thetac-gtheta(ntang))**2
       g = exp((ggamma(ntang,nta)/(dij(tri)-grij(ntang,nta))) + (ggamma(ntang,ntb)/(dik(tri)-grij(ntang,ntb))))

       vthreea = glambda(ntang)*p*g

#ifdef __DEBUG_GAROFALINI__
          write(69,'(A11,I10,I10,I10,F15.6)') 'vthree jik',j,i,k,vthreea
          write(69,'(A11,3F10.6,F25.6)') 'ij',dxij(tri),dyij(tri),dzij(tri),dij(tri)
          write(69,'(A11,3F10.6,F25.6)') 'ik',dxik(tri),dyik(tri),dzik(tri),dik(tri)
#endif

       vthree = vthree + vthreea
    end do

#ifdef __DEBUG_GAROFALINI__
    write(io_output,*) 'end Vthreebody in ',myid
#endif
  end function vthreebody

  subroutine triad()
    integer::i,j,k,ptr,ptr2

    ntr = 0
    do i=1,nchain
       do ptr = 1,neigh_cnt(i)-1
          j = neighbor(ptr,i)
          do ptr2 = ptr+1,neigh_cnt(i)
             k = neighbor(ptr2,i)
             ntr = ntr + 1

             itr1(ntr) = i
             itr2(ntr) = j
             itr3(ntr) = k

             dij(ntr) = ndij(ptr,i)
             dik(ntr) = ndij(ptr2,i)

             dxij(ntr) = nxij(ptr,i)
             dyij(ntr) = nyij(ptr,i)
             dzij(ntr) = nzij(ptr,i)

             dxik(ntr) = nxij(ptr2,i)
             dyik(ntr) = nyij(ptr2,i)
             dzik(ntr) = nzij(ptr2,i)
          end do
       end do
    end do
  end subroutine triad

!> \brief Determine correct 3-body interactions for single particles
!> particle moved: \a i
!> \param cnt number of pairs within cutoff
!> \param ni(cnt) identity of pair molecule; set in energy
  function triad_en(i,cnt,ni,nrij,nxi,nyi,nzi,lupdate) result(vthree)
    real::vthree
    integer::i,j,k,m,imolty,jmolty,kmolty,mmolty,atomj,atomk,atomm
    integer::nta,ntb,ntang,cnt,ni(:),temp_cnt,itemp,number(nmax)
    integer::temp_nei(nmax),itype,jtype,ktype,mtype
    logical::ltemp(nmax),lupdate
    real::vthreea,thetac,p,g,nrij(:),nxi(:),nyi(:),nzi(:)
    real::temp_dist(nmax),temp_x(nmax),temp_y(nmax),temp_z(nmax)

    vthree = 0.0E0_dp
    temp_cnt = 0
    do itemp = 1,nmax
       ltemp(itemp) = .false.
       temp_dist(itemp) = 0.0E0_dp
       temp_x(itemp) = 0.0E0_dp
       temp_y(itemp) = 0.0E0_dp
       temp_z(itemp) = 0.0E0_dp
    end do

    imolty = moltyp(i)
    itype = ntype(imolty,1)
    do 10 j = 1,cnt
       atomj = ni(j)
       jmolty = moltyp(atomj)
       jtype = ntype(jmolty,1)

       ! skip if i=H; go straight to l-j-i loop
       if (itype.ne.3) then
          ! skip if Si-Si or O-O
          if(jtype.eq.itype) goto 10

          ! skip if Si-H
          if(jtype.eq.3.and.itype.eq.1) goto 10

          ! loop over other pairs with i as central atom
          do 20 k = j+1,cnt
             atomk = ni(k)
             kmolty = moltyp(atomk)
             ktype = ntype(kmolty,1)

             if(ktype.eq.itype) goto 20
             if(ktype.eq.3.and.itype.eq.1) goto 20

             ! determine type
             if(jtype.eq.ktype) then
                ntang = jtype
                nta = 1
                ntb = 1
                if((nrij(j).gt.grij(ntang,1)).or. (nrij(k).gt.grij(ntang,1))) then
                   goto 20
                end if
             else
                ntang = 4
                if(jtype.eq.1) then
                   nta = 1
                   ntb = 2
                   if((nrij(j).gt.grij(4,nta)).or. (nrij(k).gt.grij(4,ntb))) then
                      goto 20
                   end if
                else
                   nta = 2
                   ntb = 1
                   if((nrij(j).gt.grij(4,nta)).or. (nrij(k).gt.grij(4,ntb))) then
                      goto 20
                   end if
                end if
             end if

             thetac = (nxi(j)*nxi(k)+nyi(j)*nyi(k)+ nzi(j)*nzi(k))/(nrij(j)*nrij(k))
             if ( thetac .ge. 1.0E0_dp ) thetac = 1.0E0_dp
             if ( thetac .le. -1.0E0_dp ) thetac = -1.0E0_dp

             p = (thetac-gtheta(ntang))**2
             g = exp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta)) )+(ggamma(ntang,ntb)/(nrij(k)-grij(ntang,ntb))))

             vthreea = glambda(ntang)*p*g

#ifdef __DEBUG_GAROFALINI__
             write(io_output,*) 'vthree jik',atomj,i,atomk,vthreea
#endif

             vthree = vthree + vthreea
             ! actual neighbor counter
             if(.not.ltemp(j).and.lupdate) then
                temp_cnt = temp_cnt+1
                number(temp_cnt) = j
                temp_nei(temp_cnt) = ni(j)
                temp_dist(temp_cnt) = nrij(j)
                temp_x(temp_cnt) = nxi(j)
                temp_y(temp_cnt) = nyi(j)
                temp_z(temp_cnt) = nzi(j)
#ifdef __DEBUG_GAROFALINI__
                write(io_output,*) 'temps a',temp_cnt,temp_nei(temp_cnt)
                write(io_output,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                write(io_output,*)
#endif
                ltemp(j) = .true.
             end if
             if(.not.ltemp(k).and.lupdate) then
                temp_cnt = temp_cnt+1
                number(temp_cnt) = k
                temp_nei(temp_cnt) = ni(k)
                temp_dist(temp_cnt) = nrij(k)
                temp_x(temp_cnt) = nxi(k)
                temp_y(temp_cnt) = nyi(k)
                temp_z(temp_cnt) = nzi(k)
#ifdef __DEBUG_GAROFALINI__
                write(io_output,*) 'temps b',temp_cnt,temp_nei(temp_cnt)
                write(io_output,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                write(io_output,*)
#endif
                ltemp(k) = .true.
             end if
20        end do
       end if

       ! now starting loop to check i-j-m pairs
       if(jtype.ne.3.and..not.(itype.eq.3.and.jtype.eq.1)) then
#ifdef __DEBUG_GAROFALINI__
          write(io_output,*)
          write(io_output,*) 'j',atomj,' neigh_cnt(j)',neigh_cnt(atomj),':', (neighbor(m,atomj),m=1,neigh_cnt(atomj))
#endif
          !> \bug why does it use neigh_cnt instead of cnt?
          do 30 m = 1,neigh_cnt(atomj)
             atomm = neighbor(m,atomj)
             mmolty = moltyp(atomm)
             mtype = ntype(mmolty,1)

             if(mtype.eq.jtype) goto 30
             if(mtype.eq.3.and.jtype.eq.1) goto 30
             if(atomm.eq.i) goto 30

             ! determine type
             if(itype.eq.mtype) then
                ntang = mtype
                nta = 1
                ntb = 1
                if((nrij(j).gt.grij(ntang,1)) .or. (ndij(m,atomj).gt.grij(ntang,1))) goto 30
             else
                ntang = 4
                if(itype.eq.1) then
                   nta = 1
                   ntb = 2
                   if((nrij(j).gt.grij(4,nta)).or. (ndij(m,atomj).gt.grij(4,ntb))) then
                      goto 30
                   end if
                else
                   nta = 2
                   ntb = 1
                   if((nrij(j).gt.grij(4,nta)).or. (ndij(m,atomj).gt.grij(4,ntb))) then
                      goto 30
                   end if
                end if
             end if
             ! because the values saved are i-j not j-i must multiply all by -1
             thetac = ((-nxi(j))*nxij(m,atomj)+(-nyi(j))* nyij(m,atomj)+(-nzi(j))*nzij(m,atomj))/ (nrij(j)*ndij(m,atomj))
             if ( thetac .ge. 1.0E0_dp ) thetac = 1.0E0_dp
             if ( thetac .le. -1.0E0_dp ) thetac = -1.0E0_dp

             p = (thetac-gtheta(ntang))**2
             g = exp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta)) )+(ggamma(ntang,ntb)/(ndij(m,atomj)- grij(ntang,ntb))))

             vthreea = glambda(ntang)*p*g

#ifdef __DEBUG_GAROFALINI__
             write(io_output,'(A13,I5,I5,I5,F15.7)') 'vthree ijm',i,atomj ,atomm,vthreea
#endif
             vthree = vthree + vthreea
             ! actual neighbor counter
             if(.not.ltemp(j).and.lupdate) then
                temp_cnt = temp_cnt+1
                number(temp_cnt) = j
                temp_nei(temp_cnt) = ni(j)
                temp_dist(temp_cnt) = nrij(j)
                temp_x(temp_cnt) = nxi(j)
                temp_y(temp_cnt) = nyi(j)
                temp_z(temp_cnt) = nzi(j)
                ltemp(j) = .true.
#ifdef __DEBUG_GAROFALINI__
                write(io_output,*) 'temps ijm',temp_cnt,temp_nei(temp_cnt)
                write(io_output,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                write(io_output,*)
#endif
             end if
30        end do
       end if
10  end do

#ifdef __DEBUG_GAROFALINI__
    write(io_output,*) 'triad_en vthree',vthree
#endif
  end function triad_en
end MODULE energy_garofalini
