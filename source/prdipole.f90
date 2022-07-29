subroutine prdipole(dip)
      implicit none
      real*8 dip(3)

      real*8 dcal

      write(*,'(''dipole moment  X         Y          Z'')')
      dcal=sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
      write(*,'(3f16.6,''  total     (au/Debye)      :'',2f11.4)')dip(1:3),dcal,dcal*2.5418    

!     dref=sqrt(dipref(1)**2+dipref(2)**2+dipref(3)**2)
!     write(*,'(3f12.6,''  reference (au/Debye): '',2f10.5)')dipref(1:3),dref,dref*2.5418    
!     if(dcal.gt.0.05)then
!     dang = 180d0*acos(sum(dip*dipref)/(1.d-8+dcal*dref))/3.14159d0
!     write(*,'(''angle  (Dcal-Dref) deg: '', f8.3)')dang
!     endif

end

subroutine prsec(sec)
      implicit none
      real*8 sec(3)

      real*8 dcal

      write(*,'(''second moment  X^2       Y^2        Z^2'')')
      dcal=(sec(1)+sec(2)+sec(3))/3.d0
      write(*,'(3f16.6,''  average   (au)            :'',f11.4)')sec(1:3),dcal

end


subroutine prpolar(alp)
      implicit none
      real*8 alp(6)

      real*8 av  

      av=(alp(1)+alp(3)+alp(6))/3d0    
      write(*,'(''dipole polarizability      X         Y           Z'',5x,''total (au/A^3) :'',2f11.4)')av,av*0.148185
      write(*,'(15x,'' X '',3f16.6)') alp(1)
      write(*,'(15x,'' Y '',3f16.6)') alp(2:3)
      write(*,'(15x,'' Z '',3f16.6)') alp(4:6)


end

subroutine prbeta(beta)
      implicit none
      real*8 beta(6,3)
      integer i,j,k
      real*8 tmp(3,3,3),av(3)
      character*3 s1,s2,s3
      character*3 xyz
      xyz='XYZ'

      k=0
      do i=1,3
         do j=1,i
            k=k+1
            tmp(j,i,1:3)=beta(k,1:3)
            tmp(i,j,1:3)=beta(k,1:3)
         enddo
      enddo
      
      av = 0
      do i=1,3
         do j=1,3
            av(i) = av(i) + tmp(i,j,j)
         enddo
      enddo

      write(*,'(''dipole hyperpolarizability (au):'')')
      open(unit=21,file='.data_beta')
      do i=1,3
         do j=1,3
            s1(3:3)=xyz(i:i)
            s2(3:3)=xyz(i:i)
            s3(3:3)=xyz(i:i)
            s1(1:1)=xyz(1:1)
            s2(1:1)=xyz(2:2)
            s3(1:1)=xyz(3:3)
            s1(2:2)=xyz(j:j)
            s2(2:2)=xyz(j:j)
            s3(2:2)=xyz(j:j)
            write(*,'(3(a3,f14.4,5x))') s1,tmp(1,j,i),s2,tmp(2,j,i),s3,tmp(3,j,i)
            do k=1,3
               write(21,'(F14.6)') tmp(k,j,i)
            enddo
         enddo
      enddo
      close(21)
      write(*,'(''sum beta_i  X Y Z : '',3f12.3)')av
            
end
