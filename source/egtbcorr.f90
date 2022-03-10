!--------------------------------------------------------
!--------------------------------------------------------

subroutine egtbcorr(n,at,rab,wbo,e)
      use parcom
      use com
      implicit none
      integer n
      integer,intent(in)     :: at(n)       
      real*8, intent(in)     :: rab(n*(n+1)/2) 
      real*8, intent(in)     :: wbo(n,n)                  ! WBO 
      real*8, intent(out)    :: e                         ! total energy correction

      real*8 wbocut,w                   
      integer i,j,k,ati,atj

      e = 0

      k = 0
      do i=1, n
         ati = at(i)
         do j=1, i - 1
            atj = at(j)
            wbocut = ener_par1(9,atj)+ener_par1(9,ati)
            w = exp(-glob_par(1)*(wbo(j,i)+wbocut)**2)
            e = e + w * (ener_par1(10,atj)+ener_par1(10,ati))
         enddo
      enddo

      e = -0.001*e

      open(unit=112,file='.CP')
      write(112,*) e
      close(112)

      end
