      subroutine dipole(ibox,mtype)
      use sim_system
      implicit none

      integer::ibox,mtype,i,imolty,zzz,ii
      real::dipox(2),dipoy(2),dipoz(2)
      
      if ( mtype .eq. 0 ) then
! in sumup, initiliaze dipole to be 0 and then sum up all the dipoles

         dipolex(ibox) = 0.0E0_dp
         dipoley(ibox) = 0.0E0_dp
         dipolez(ibox) = 0.0E0_dp
         
         do i = 1, nchain
            if ( nboxi(i) .eq. ibox ) then
               imolty = moltyp(i)
               do ii = 1,nunit(imolty)
                  dipolex(ibox) = dipolex(ibox)  + qqu(i,ii)*rxu(i,ii)
                  dipoley(ibox) = dipoley(ibox)  + qqu(i,ii)*ryu(i,ii)
                  dipolez(ibox) = dipolez(ibox)  + qqu(i,ii)*rzu(i,ii)
               end do
            end if
         end do

      else if(mtype .eq. 1 ) then
         
! calculate the dipole moment after the traslation, rotation and 
! charge move

         do zzz = 1,2
            dipox(zzz) = 0.0E0_dp
            dipoy(zzz) = 0.0E0_dp
            dipoz(zzz) = 0.0E0_dp
            imolty = moltion(zzz)
            do i = 1,nunit(imolty)
               dipox(zzz) = dipox(zzz) +  qquion(i,zzz)*rxuion(i,zzz)
               dipoy(zzz) = dipoy(zzz) + qquion(i,zzz)*ryuion(i,zzz)
               dipoz(zzz) = dipoz(zzz) + qquion(i,zzz)*rzuion(i,zzz)
            end do
         end do
         dipolex(ibox) = dipolex(ibox) - dipox(1) + dipox(2) 
         dipoley(ibox) = dipoley(ibox) - dipoy(1) + dipoy(2)
         dipolez(ibox) = dipolez(ibox) - dipoz(1) + dipoz(2)

      else if(mtype .eq. 2) then

! store the old dipole moment

         dipolexo = dipolex(ibox) 
         dipoleyo = dipoley(ibox)
         dipolezo = dipolez(ibox)
         
      else if(mtype .eq. 3) then

! restore the old dipole moment

         dipolex(ibox) = dipolexo 
         dipoley(ibox) = dipoleyo
         dipolez(ibox) = dipolezo

      end if
      
       
      return
      end
