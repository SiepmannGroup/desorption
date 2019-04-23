!> \brief Setup Maginns interphase switch
      subroutine ipswsetup
      use util_runtime,only:err_exit
      use sim_system
      use sim_cell
      implicit none

      integer::i,j,k,tnw,ibox

      real::lx,hmata(9),hmatc(9),hm(9)

      logical::lhm,llwell

      ibox = 1
      if (lmipsw.and.(nbox.gt.1)) call err_exit(__FILE__,__LINE__,'ipsw only for 1 box',myid+1)
      if (lmipsw.and.lnpt.and.pmvol.gt.1.0E-7_dp) call err_exit(__FILE__,__LINE__,'ipsw only for NVT',myid+1)
      read(35,*)
      read(35,*) (lwell(i),i=1,nmolty)
      read(35,*)
      tnw = 0
      llwell = .false.
      do i = 1, nmolty
         if (lwell(i)) then
            llwell = .true.
            nwell(i) = temtyp(i)
            tnw = tnw+nwell(i)*nunit(i)
            do j = 1, nunit(i)
               read(35,*) (awell(j,k,i), k = 1,nunit(i))
            end do
         end if
      end do
      if (nw.lt.tnw) then
         write(io_output,*) 'increase nw in ipswpar to ', tnw
         call err_exit(__FILE__,__LINE__,'',myid+1)
      end if
      read(35,*)
      read(35,*) bwell
      read(35,*)
      read(35,*) lstagea, lstageb, lstagec
      if (lmipsw) then
         if ((.not.lstagea).and.(.not.lstageb).and.(.not.lstagec)) call err_exit(__FILE__,__LINE__,'one stage must be true',myid+1)
      end if
      if (llwell.and.lstagea) call err_exit(__FILE__,__LINE__,'ipsw well NOT for stage a',myid+1)
      if (lstagea) then
         if (lstageb.or.lstagec) call err_exit(__FILE__,__LINE__,'only one lstage must be true',myid+1)
      end if
      if (lstageb) then
         if (lstagea.or.lstagec) call err_exit(__FILE__,__LINE__,'only one lstage must be true',myid+1)
      end if
      if (lstagec) then
         if (lstagea.or.lstageb) call err_exit(__FILE__,__LINE__,'only one lstage must be true',myid+1)
      end if
      read(35,*)
      read(35,*) etais, lambdais
      if ((lambdais.le.-1.0E-6_dp).or.(lambdais.ge.1.000001)) call err_exit(__FILE__,__LINE__,'lambdais must be between 0 and 1',myid+1)
      if ((etais.le.-1.0E-6_dp).or.(etais.ge.1.000001)) call err_exit(__FILE__,__LINE__,'etais must be between 0 and 1',myid+1)
      read(35,*)
      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         read(35,*) hmata(1),hmata(2),hmata(3)
         read(35,*) hmata(4),hmata(5),hmata(6)
         read(35,*) hmata(7),hmata(8),hmata(9)
         read(35,*) hmatc(1),hmatc(2),hmatc(3)
         read(35,*) hmatc(4),hmatc(5),hmatc(6)
         read(35,*) hmatc(7),hmatc(8),hmatc(9)
      else if (.not.lsolid(ibox).and.(.not.lrect(ibox))) then
         read(35,*) lena, lenc
      end if
      if (lmipsw.and.lstageb) then
         if (lsolid(ibox).and.(.not.lrect(ibox))) then
            lhm = .false.
            do i = 1, 9
               hm(i) = (1.0E0_dp-lambdais)*hmata(i)+lambdais*hmatc(i)
               if (abs(hmat(ibox,i)-hm(i)).gt.1.0E-6_dp) lhm = .true.
            end do
            if (lhm) then
               write(io_output,*) 'enter correct hmat'
               write(io_output,*) hm(1),hm(2),hm(3)
               write(io_output,*) hm(4),hm(5),hm(6)
               write(io_output,*) hm(7),hm(8),hm(9)
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if
         else if (.not.lsolid(ibox)) then
            lx = (1.0E0_dp-lambdais)*lena+lambdais*lenc
            if (abs(boxlx(1)-lx).gt.1.0E-6_dp) then
               write(io_output,*) 'input correct boxl', lx
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if
         end if
      end if
      if (lmipsw.and.lstagea) then
         if (lsolid(ibox).and.(.not.lrect(ibox))) then
            lhm = .false.
            do i = 1, 9
               if (abs(hmat(ibox,i)-hmata(i)).gt.1.0E-6_dp) lhm = .true.
            end do
            if (lhm) then
               write(io_output,*) 'enter correct hmat'
               write(io_output,*) hmata(1),hmata(2),hmata(3)
               write(io_output,*) hmata(4),hmata(5),hmata(6)
               write(io_output,*) hmata(7),hmata(8),hmata(9)
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if
!write(io_output,*) boxlx(1),lena
         else if (.not.lsolid(ibox)) then
            if (abs(boxlx(1)-lena).gt.1.0E-6_dp) then
               write(io_output,*) 'input correct boxl', lena
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if
         end if
      end if
      if (lmipsw.and.lstagec) then
         if (lsolid(ibox).and.(.not.lrect(ibox))) then
            lhm = .false.
            do i = 1, 9
               if (abs(hmat(ibox,i)-hmatc(i)).gt.1.0E-6_dp) lhm = .true.
            end do
            if (lhm) then
               write(io_output,*) 'enter correct hmat'
               write(io_output,*) hmatc(1),hmatc(2),hmatc(3)
               write(io_output,*) hmatc(4),hmatc(5),hmatc(6)
               write(io_output,*) hmatc(7),hmatc(8),hmatc(9)
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if
         else if (.not.lsolid(ibox)) then
            if (abs(boxlx(1)-lenc).gt.1.0E-6_dp) then
               write(io_output,*) 'input correct boxl', lenc
               call err_exit(__FILE__,__LINE__,'',myid+1)
            end if
         end if
      end if
      dhmat(1,1) = hmatc(1)-hmata(1)
      dhmat(1,2) = hmatc(4)-hmata(4)
      dhmat(1,3) = hmatc(7)-hmata(7)
      dhmat(2,1) = hmatc(2)-hmata(2)
      dhmat(2,2) = hmatc(5)-hmata(5)
      dhmat(2,3) = hmatc(8)-hmata(8)
      dhmat(3,1) = hmatc(3)-hmata(3)
      dhmat(3,2) = hmatc(6)-hmata(6)
      dhmat(3,3) = hmatc(9)-hmata(9)
      read(35,*)
      read(35,*) iratipsw
      if (lmipsw.and.lstageb) then
         if (mod(iratipsw,iratp).ne.0) call err_exit(__FILE__,__LINE__,'iratipsw must be integer multiple of iratp if lstageb is true',myid+1)
      end if
      do i = 1, nmolty
         if (lwell(i)) then
            do j = 1, nwell(i)*nunit(i)
               read(35,*) sxwell(j,i),sywell(j,i),szwell(j,i)
            end do
         end if
      end do
      if (lmipsw) then
         write(io_output,*) '*****************************************'
         write(io_output,*) 'Some parameters used for interphase switch'
         write(io_output,*) 'moltyp, lwell, and nwell'
         do i = 1, nmolty
            if (lwell(i)) then
               write(io_output,*) i, lwell(i), nwell(i)
            else
               write(io_output,*) i, lwell(i), 'not defined'
            end if
         end do
         write(io_output,*) 'lstagea,lstageb,lstagec', lstagea,lstageb,lstagec
         write(io_output,*) 'etais, lambdais', etais,lambdais
         write(io_output,*) '*****************************************'
      else
         lstagea = .false.
         lstageb = .false.
         lstagec = .false.
         do i = 1, nmolty
            lwell(i) = .false.
         end do
         lambdais = 0.0E0_dp
      end if

      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         do i = 1, nmolty
            if (lwell(i)) then
               do j = 1, nwell(i)*nunit(i)
                  rxwell(j,i) = sxwell(j,i)*hmat(ibox,1)+ sywell(j,i)*hmat(ibox,4)+szwell(j,i)*hmat(ibox,7)
                  rywell(j,i) = sxwell(j,i)*hmat(ibox,2)+ sywell(j,i)*hmat(ibox,5)+szwell(j,i)*hmat(ibox,8)
                  rzwell(j,i) = sxwell(j,i)*hmat(ibox,3)+ sywell(j,i)*hmat(ibox,6)+szwell(j,i)*hmat(ibox,9)
               end do
            end if
         end do
      else
         do i = 1, nmolty
            if (lwell(i)) then
               do j = 1, nwell(i)*nunit(i)
                  rxwell(j,i) = sxwell(j,i)*boxlx(ibox)
                  rywell(j,i) = sywell(j,i)*boxly(ibox)
                  rzwell(j,i) = szwell(j,i)*boxlz(ibox)
               end do
            end if
         end do
      end if

      return
      end
