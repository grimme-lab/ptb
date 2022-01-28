
subroutine calcrab(n,at,xyz,rab)
      implicit none
      integer n
      integer,intent(in)     :: at(n)       
      real*8, intent(in)     :: xyz(3,n)
      real*8, intent(out)    :: rab(n*(n+1)/2)

      integer i,j,k
      real*8 rab2, r

      k=0
      do i=1,n
         do j=1,i-1
            rab2= &
     &              (xyz(1,i)-xyz(1,j))**2+&
     &              (xyz(2,i)-xyz(2,j))**2+&
     &              (xyz(3,i)-xyz(3,j))**2
            k = k + 1
            r = sqrt(rab2)
            rab(k) = r
         enddo
         k = k +1
         rab(k)=1.d-12 ! avoid NANs
      enddo

end
