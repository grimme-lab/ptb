      subroutine setespot(nshell,qsh,jab,ves)
      implicit none 
      integer, intent(in) :: nshell
      real*8, intent(in) ::  qsh(nshell),jab(nshell,nshell)
      real*8, intent(inout) :: ves(nshell) 
      real*8 qshi,vesi
      integer i,j

      ves = 0
      do i=1,nshell
         qshi=qsh(i)
         vesi=0.0d0
         do j=1,i-1
            ves(j)=ves(j)+qshi*jab(j,i) 
            vesi=vesi+qsh(j)*jab(j,i)
         enddo
         vesi=vesi+qshi*jab(i,i)
         ves(i)=ves(i)+vesi
      enddo
      end  
