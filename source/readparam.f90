subroutine rdelparam(fname)
use parcom
      implicit none
      character*80 fname    

      integer ich
      integer iz  
      integer nshell
      integer i,nn  
      real*8 xx(10)
      character*80 atmp     

      ich=11
      open(unit=ich,file=fname)

 10   read(ich,'(a)',end=100) atmp

      if(index(atmp,'Z=').ne.0)then
         call readl(atmp,xx,nn)
!        ordinal number         
         iz    =idint(xx(1))
         ishell=idint(xx(2))
 20      read(ich,'(a)',end=100) atmp
         if(index(atmp,'lev=').ne.0)then
            call readl(atmp,xx,nn)
         endif
         if(index(atmp,'$end').ne.0) goto 10
      endif


100   close(ich)

end
