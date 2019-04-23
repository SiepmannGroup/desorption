MODULE sim_cell
  use var_type,only:dp,RealPtr
  use util_runtime,only:err_exit
  use sim_system,only:nmax
  implicit none
  private
  save
  public::CellType,CellMaskType,matops,setpbc,mimage,build_linked_cell,update_linked_cell,get_cell_neighbors,allocate_sim_cell

  type CellType
     logical::ortho& !<.true. if the simulation box is orthorhombic
      ,solid !<.true. if the simulation box can change shape
     logical,dimension(3)::pbc
     real::cut,vol,calp,boxl(3),ang(3),height(3),hmat(3,3),hmati(3,3)
  end type CellType

  type CellMaskType
     logical,pointer::ortho& !<.true. if the simulation box is orthorhombic
      ,solid!<.true. if the simulation box can change shape
     logical,dimension(3)::pbc
     real,pointer::cut,vol,calp
     type(RealPtr)::boxl(3),ang(3),height(3),hmat(3,3),hmati(3,3)
  end type CellMaskType

  real,allocatable,target,public::hmat(:,:),hmati(:,:),cell_length(:,:),min_width(:,:),cell_vol(:),cell_ang(:,:)

  real::bx,by,bz,hbx,hby,hbz,bxi,byi,bzi

  integer,parameter::cmax = 1000& !< cmax is the max number of cells for the linkcell list
   ,cmaxa = 512 !< cmaxa is the max number of molecules per cell
  real::dcellx,dcelly,dcellz
  integer::nicell(cmax),iucell(cmax,cmaxa),ncellx,ncelly,ncellz
  integer,allocatable::icell(:)

contains
!> \brief Calculates for non cubic simulation cell
!> \b cell_length: boxlengths \n
!> \b min_length: minimum boxwidths \n
!> \b cell_vol: box volume \n
!> \b hmati: inverse H matrix
!> \author Neeraj Rai (in Merck Apr 2005)
  subroutine matops(ibox)
    use sim_system,only:boxlx,boxly,boxlz
    integer,intent(in)::ibox
    real::abx,aby,abz,bcx,bcy,bcz,cax,cay,caz
    real::elem(9),inv_vol,adj(9),cosa,cosb,cosg
    integer::i

    do i=1,9
       elem(i)=hmat(ibox,i)
    end do

    ! calculating the length of cell vectors
    cell_length(ibox,1)=sqrt(elem(1)*elem(1)+elem(2)*elem(2)+ elem(3)*elem(3))
    cell_length(ibox,2)=sqrt(elem(4)*elem(4)+elem(5)*elem(5)+ elem(6)*elem(6))
    cell_length(ibox,3)=sqrt(elem(7)*elem(7)+elem(8)*elem(8)+ elem(9)*elem(9))

    boxlx(ibox) = cell_length(ibox,1)
    boxly(ibox) = cell_length(ibox,2)
    boxlz(ibox) = cell_length(ibox,3)

    ! calculating cross product of cell vectors
    abx=elem(2)*elem(6)-elem(3)*elem(5)
    aby=elem(3)*elem(4)-elem(1)*elem(6)
    abz=elem(1)*elem(5)-elem(2)*elem(4)
    bcx=elem(5)*elem(9)-elem(6)*elem(8)
    bcy=elem(6)*elem(7)-elem(4)*elem(9)
    bcz=elem(4)*elem(8)-elem(5)*elem(7)
    cax=elem(8)*elem(3)-elem(2)*elem(9)
    cay=elem(1)*elem(9)-elem(3)*elem(7)
    caz=elem(2)*elem(7)-elem(1)*elem(8)

    ! calculating cell volume
    cell_vol(ibox) = elem(1)*bcx+elem(2)*bcy+elem(3)*bcz

    if(abs(cell_vol(ibox)).lt.1E-16_dp) then
       call err_exit(__FILE__,__LINE__,'Volume of cell negligible, check input H matrix',-1)
    end if

    inv_vol = 1.0E0_dp/cell_vol(ibox)

    ! calculating minimum cell widths
    min_width(ibox,1) = cell_vol(ibox)/sqrt(bcx*bcx+bcy*bcy+ bcz*bcz)
    min_width(ibox,2) = cell_vol(ibox)/sqrt(cax*cax+cay*cay+ caz*caz)
    min_width(ibox,3) = cell_vol(ibox)/sqrt(abx*abx+aby*aby+ abz*abz)

    ! calculating adjoint for inverting the h-matrix
    adj(1)=elem(5)*elem(9)-elem(6)*elem(8)
    adj(2)=elem(3)*elem(8)-elem(2)*elem(9)
    adj(3)=elem(2)*elem(6)-elem(3)*elem(5)
    adj(4)=elem(6)*elem(7)-elem(4)*elem(9)
    adj(5)=elem(1)*elem(9)-elem(3)*elem(7)
    adj(6)=elem(3)*elem(4)-elem(1)*elem(6)
    adj(7)=elem(4)*elem(8)-elem(5)*elem(7)
    adj(8)=elem(2)*elem(7)-elem(1)*elem(8)
    adj(9)=elem(1)*elem(5)-elem(2)*elem(4)

    ! inverting the matrix
    hmati(ibox,1) = inv_vol * adj(1)
    hmati(ibox,2) = inv_vol * adj(2)
    hmati(ibox,3) = inv_vol * adj(3)
    hmati(ibox,4) = inv_vol * adj(4)
    hmati(ibox,5) = inv_vol * adj(5)
    hmati(ibox,6) = inv_vol * adj(6)
    hmati(ibox,7) = inv_vol * adj(7)
    hmati(ibox,8) = inv_vol * adj(8)
    hmati(ibox,9) = inv_vol * adj(9)

    ! calculating alpha, beta and gamma using the dot product rule
    ! cell_ang(ibox,1)=alpha;2= beta; 3= gamma
    ! In crystallography Literature angle between b & c is alpha (1),
    ! angle between c & a is Beta (2) and angle between a & b is gamma (3)
    cosa = (elem(4)*elem(7)+elem(5)*elem(8)+elem(6)*elem(9))/ (cell_length(ibox,2)*cell_length(ibox,3))
    cosb = (elem(1)*elem(7)+elem(2)*elem(8)+elem(3)*elem(9))/ (cell_length(ibox,1)*cell_length(ibox,3))
    cosg = (elem(4)*elem(1)+elem(5)*elem(2)+elem(6)*elem(3))/ (cell_length(ibox,2)*cell_length(ibox,1))
    cell_ang(ibox,1) = acos(cosa)
    cell_ang(ibox,2) = acos(cosb)
    cell_ang(ibox,3) = acos(cosg)

    return
  end subroutine matops

  subroutine setpbc(ibox)
    use sim_system,only:boxlx,boxly,boxlz,lpbcx,lpbcy,lpbcz,lfold
    integer,intent(in)::ibox

    if ( lpbcx ) then
       bx = boxlx(ibox)
       if ( lfold ) then
          hbx = 0.5E0_dp * bx
       else
          bxi = 1.0E0_dp / bx
       end if
    end if

    if ( lpbcy ) then
       by = boxly(ibox)
       if ( lfold ) then
          hby = 0.5E0_dp * by
       else
          byi = 1.0E0_dp / by
       end if
    end if

    if ( lpbcz ) then
       bz = boxlz(ibox)
       if ( lfold ) then
          hbz = 0.5E0_dp * bz
       else
          bzi = 1.0E0_dp / bz
       end if
    end if

    return
  end subroutine setpbc

  pure subroutine mimage(rxuij,ryuij,rzuij,ibox)
    use sim_system,only:lsolid,lrect,lpbcx,lpbcy,lpbcz,lfold,ltwice
    integer,intent(in)::ibox
    real,intent(inout)::rxuij,ryuij,rzuij

    real::hsx,hsy,hsz,sx,sy,sz

    if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
       ! non-hexagonal box
       hsx = 0.5E0_dp*hmat(ibox,1)
       hsy = 0.5E0_dp*hmat(ibox,5)
       hsz = 0.5E0_dp*hmat(ibox,9)
       ! sx, sy, sz are the coordinates of vector (rxuij,ryuij,rzuij) in the
       ! basis (n1,n2,n3), which is the transpose of the H matrix, and is the
       ! transforming matrix from basis (n1,n2,n3) to the canonical basis
       ! (e1,e2,e3), where e1, e2, e3 are the three perpendicular unit vector.
       sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox,4) +rzuij*hmati(ibox,7)
       sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox,5) +rzuij*hmati(ibox,8)
       sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox,6) +rzuij*hmati(ibox,9)

       !         if ( sx .gt. 0.5E0_dp ) then
       !            sx = sx-1E0_dp
       !         else if ( sx .lt. -0.5E0_dp ) then
       !            sx = sx+1E0_dp
       !         end if
       !         if ( sy .gt. 0.5E0_dp ) then
       !            sy = sy-1E0_dp
       !         else if ( sy .lt. -0.5E0_dp ) then
       !            sy = sy+1E0_dp
       !         end if
       !         if ( sz .gt. 0.5E0_dp ) then
       !            sz = sz-1E0_dp
       !         else if ( sz .lt. -0.5E0_dp ) then
       !            sz = sz+1E0_dp
       !         end if
       !         sx = sx-nint(sx)
       !         sy = sy-nint(sy)
       !         sz = sz-nint(sz)
       !         rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
       !         ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+sz*hmat(ibox,8)
       !         rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+sz*hmat(ibox,9)

       !         print*,rxuij,ryuij,rzuij

       ! Here it implies that in the H matrix, the first vector must only have
       ! the x component, the second vector the x,y components, and only the
       ! third can have all the x,y,z components
       if ( rzuij .gt. hsz ) then
          rzuij=rzuij-hmat(ibox,9)
          sz=sz-1E0_dp
          if ( rzuij .gt. hsz ) then
             rzuij=rzuij-hmat(ibox,9)
             sz=sz-1E0_dp
          end if
       else if ( rzuij .lt. -hsz ) then
          rzuij=rzuij+hmat(ibox,9)
          sz=sz+1E0_dp
          if ( rzuij .lt. -hsz ) then
             rzuij=rzuij+hmat(ibox,9)
             sz=sz+1E0_dp
          end if
       end if

       ryuij=sy*hmat(ibox,5)+sz*hmat(ibox,8)
       if ( ryuij .gt. hsy ) then
          ryuij=ryuij-hmat(ibox,5)
          sy=sy-1E0_dp
          if ( ryuij .gt. hsy ) then
             ryuij=ryuij-hmat(ibox,5)
             sy=sy-1E0_dp
          end if
       else if ( ryuij .lt. -hsy ) then
          ryuij=ryuij+hmat(ibox,5)
          sy=sy+1E0_dp
          if ( ryuij .lt. -hsy ) then
             ryuij=ryuij+hmat(ibox,5)
             sy=sy+1E0_dp
          end if
       end if

       rxuij=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
       if ( rxuij .gt. hsx ) then
          rxuij=rxuij-hmat(ibox,1)
          if ( rxuij .gt. hsx ) then
             rxuij=rxuij-hmat(ibox,1)
          end if
       else if ( rxuij .lt. -hsx ) then
          rxuij=rxuij+hmat(ibox,1)
          if ( rxuij .lt. -hsx ) then
             rxuij=rxuij+hmat(ibox,1)
          end if
       end if
       !         print*, sqrt(rxuij**2+ryuij**2+rzuij**2)
       !         print*
    else
       ! orthorhombic box
       if ( lpbcx ) then
          if ( lfold ) then
             if ( rxuij .gt. hbx ) then
                rxuij=rxuij-bx
             else
                if (rxuij.lt.-hbx) rxuij=rxuij+bx
             end if
          else
             ! rxuij = rxuij - bx*anint(rxuij*bxi)
             rxuij = rxuij - bx*aint(rxuij*bxi+sign(0.5E0_dp,rxuij))
          end if
       end if

       if ( lpbcy ) then
          if ( lfold ) then
             if ( ryuij .gt. hby ) then
                ryuij=ryuij-by
             else
                if (ryuij.lt.-hby) ryuij=ryuij+by
             end if
          else
             ! ryuij  = ryuij - by*anint(ryuij*byi)
             ryuij = ryuij - by*aint(ryuij*byi+sign(0.5E0_dp,ryuij))
          end if
       end if

       if ( lpbcz ) then
          if ( lfold ) then
             if (rzuij.gt.hbz) then
                rzuij=rzuij-bz
             else
                if (rzuij.lt.-hbz) rzuij=rzuij+bz
             end if
          else
             ! rzuij  = rzuij - bz*anint(rzuij*bzi)
             rzuij = rzuij - bz*aint(rzuij*bzi+sign(0.5E0_dp,rzuij))
          end if
       end if

       ! --- JLR 11-17-09
       ! --- for orginal RPLC setup, mimage must be applied twice or energy errors can result
       ! --- this may be removed in future if larger system size is used!!!
       if (ltwice(ibox)) then !everthing is folded again
          if ( rxuij .gt. hbx ) then
             rxuij=rxuij-bx
          else
             if (rxuij.lt.-hbx) rxuij=rxuij+bx
          end if
          if ( ryuij .gt. hby ) then
             ryuij=ryuij-by
          else
             if (ryuij.lt.-hby) ryuij=ryuij+by
          end if
       end if
       ! --- END JLR 11-17-09
    end if
    return
  end subroutine mimage

!> \brief Decodes x,y,z to a single number
!DEC$ ATTRIBUTES FORCEINLINE :: linkdecode
  pure function linkdecode(i,j,k,ncellx,ncelly,ncellz)
    integer::linkdecode
    integer,intent(in)::i,j,k,ncellx,ncelly,ncellz

    linkdecode = (i-1)*ncelly*ncellz + (j-1)*ncellz + k
  end function linkdecode

!> \brief Set up or update linked cell list
  subroutine build_linked_cell()
    use sim_system,only:xcm,ycm,zcm,boxlink,rintramax,rcut,boxlx,boxly,boxlz,io_output,nchain,nboxi
    integer::ibox,ncell,n,i,j,k,ic
    integer,save::ncello=0

    ibox = boxlink

! determine dcell
! write(io_output,*) 'linkcell used',rcut,rintramax
! rintramax is the maximum distance between endpoints in a molecule
    dcellx = rcut(ibox) + rintramax

! find hypothetical ncell
    ncellx = int( boxlx(ibox) / dcellx )
    ncelly = int( boxly(ibox) / dcellx )
    ncellz = int( boxlz(ibox) / dcellx )

! make dcells larger so each each cell is the same size
    dcellx = boxlx(ibox) / dble(ncellx)
    dcelly = boxly(ibox) / dble(ncelly)
    dcellz = boxlz(ibox) / dble(ncellz)

! now reweight ncell one more time
    ncellx = anint( boxlx(ibox) / dcellx )
    ncelly = anint( boxly(ibox) / dcelly )
    ncellz = anint( boxlz(ibox) / dcellz )

    ncell = ncellx * ncelly * ncellz

    if (ncell .ne. ncello) then
       if (ncello .eq. 0) then
          write(io_output,*) 'number of linkcells set to',ncell
       else
          write(io_output,*) 'number of linkcells changed to',ncell
       end if
    end if

    ncello = ncell

    if (ncell.gt.cmax) then
       write(io_output,*) 'ncell,cmax',ncell,cmax
       call err_exit(__FILE__,__LINE__,'ncell greater than cmax in linkcell',-1)
    end if

    do n = 1, ncell
       nicell(n) = 0
    end do

! assign molecules to cells
    do n = 1, nchain
       if (nboxi(n).eq.boxlink) then
          i = int(xcm(n) / dcellx) + 1
          j = int(ycm(n) / dcelly) + 1
          k = int(zcm(n) / dcellz) + 1
          ic = linkdecode(i,j,k,ncellx,ncelly,ncellz)

          if (ic.gt.cmax) then
             write(io_output,*) 'ic,cmax',ic,cmax
             call err_exit(__FILE__,__LINE__,'ic gt cmax',-1)
          end if

          icell(n) = ic
          nicell(ic) = nicell(ic) + 1

          if (nicell(ic).gt.cmaxa) then
             write(io_output,*) 'nicell,cmaxa',nicell(ic) ,cmaxa
             call err_exit(__FILE__,__LINE__,'nicell gt cmaxa',-1)
          end if

          iucell(ic,nicell(ic)) = n
       else
          icell(n) = 0
       end if
    end do
  end subroutine build_linked_cell

!> \brief Update cell's occupants
  subroutine update_linked_cell(imol)
    use sim_system,only:xcm,ycm,zcm,nboxi,boxlink
    integer,intent(in)::imol
    integer::ic,ico,i,j,k,n

    if (nboxi(imol).eq.boxlink) then
       i = int(xcm(imol) / dcellx) + 1
       j = int(ycm(imol) / dcelly) + 1
       k = int(zcm(imol) / dcellz) + 1
       ic = linkdecode(i,j,k,ncellx,ncelly,ncellz)
    else
       ic=0
    end if

    ico=icell(imol)
    if (ic.ne.ico) then
       if (ico.gt.0) then
! first remove our molecule
          do n = 1, nicell(ico)
             if (iucell(ico,n).eq.imol) then
! replace removed occupant with last occupant and erase last spot
                iucell(ico,n) = iucell(ico,nicell(ico))
                iucell(ico,nicell(ico)) = 0
                icell(imol) = 0
                exit
             end if
          end do

          if (n.gt.nicell(ico)) then
             call err_exit(__FILE__,__LINE__,'screwup for iinit = 2 for linkcell',-1)
          else
             nicell(ico) = nicell(ico) - 1
          end if
       end if

       if (ic.gt.0) then
! now we will add the molecule
          icell(imol) = ic
          nicell(ic) = nicell(ic) + 1
          if (nicell(ic).gt.cmaxa) call err_exit(__FILE__,__LINE__,'nicell too big',-1)
          iucell(ic,nicell(ic)) = imol
       end if
    end if
  end subroutine update_linked_cell

!> \brief Determine the cell neighbors
  subroutine get_cell_neighbors(xcmi,ycmi,zcmi,ibox,jcell,nmole)
    use sim_system,only:boxlx,boxly,boxlz
    integer,intent(in)::ibox
    real,intent(inout)::xcmi,ycmi,zcmi
    integer,intent(out)::jcell(nmax),nmole
    integer::i,j,k,ia,ja,ka,ib,jb,kb,ic,n

! check perodic boundaries
    if (xcmi.gt.boxlx(ibox)) then
       xcmi = xcmi - boxlx(ibox)
    else if (xcmi.lt.0) then
       xcmi = xcmi + boxlx(ibox)
    end if

    if (ycmi.gt.boxly(ibox)) then
       ycmi = ycmi - boxly(ibox)
    else if (ycmi.lt.0) then
       ycmi = ycmi + boxly(ibox)
    end if

    if (zcmi.gt.boxlz(ibox)) then
       zcmi = zcmi - boxlz(ibox)
    else if (zcmi.lt.0) then
       zcmi = zcmi + boxlz(ibox)
    end if

    i = int(xcmi / dcellx) + 1
    j = int(ycmi / dcelly) + 1
    k = int(zcmi / dcellz) + 1

    nmole = 0
    do ia = i-1, i+1
       do ja = j-1, j+1
          do ka = k-1, k+1
             if (ia.gt.ncellx) then
                ib = ia - ncellx
             else if (ia.lt.1) then
                ib = ia + ncellx
             else
                ib = ia
             end if

             if (ja.gt.ncelly) then
                jb = ja - ncelly
             else if (ja.lt.1) then
                jb = ja + ncelly
             else
                jb = ja
             end if

             if (ka.gt.ncellz) then
                kb = ka - ncellz
             else if (ka.lt.1) then
                kb = ka + ncellz
             else
                kb = ka
             end if

             ic = linkdecode(ib,jb,kb,ncellx,ncelly,ncellz)
             do n = 1, nicell(ic)
                nmole = nmole + 1
                jcell(nmole) = iucell(ic,n)
             end do
          end do
       end do
    end do
  end subroutine get_cell_neighbors

  subroutine allocate_sim_cell
    use sim_system,only:nbxmax
    integer::jerr

    if (allocated(hmat)) deallocate(hmat,hmati,cell_length,min_width,cell_vol,cell_ang,stat=jerr)
    allocate(hmat(nbxmax,9),hmati(nbxmax,9),cell_length(nbxmax,3),min_width(nbxmax,3),cell_vol(nbxmax),cell_ang(nbxmax,3),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'allocate_sim_cell: allocation failed',jerr)

    hmat=0.0E0_dp

    if (allocated(icell)) deallocate(icell,stat=jerr)
    allocate(icell(nmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'allocate_sim_cell.linked_cell: allocation failed',jerr)
  end subroutine allocate_sim_cell
end MODULE sim_cell
