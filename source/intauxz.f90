subroutine intauxz(nat,at,xyz,v1)
      use bascom
      use parcom
      implicit none          
      integer, intent(in)  :: nat
      integer, intent(in)  :: at(nat)     
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(out)  :: v1(naux)

      real*8 alp,pqx,pqy,pqz,ab,eab,eabo,f,fn,pq1,pq2,pq3,dx,dy,dz,cpsq,v00,f02
      real*8 s00

      real*8, parameter ::       PI = 3.1415926535897932384626433832795029d0
      real*8, parameter ::       PI2= 2.0d0/PI
      real*8, parameter ::       TPI= 2.0d0*PI

      integer i,k
      integer iat

      s00 = 0
      do i = 1, naux
         iat = aux_at (i)
         alp = sqrt(aux_exp(i)) ! aux function describes rho so take sqrt when the (s|Z/r|s) formula is used
         pqx =alp*xyz(1,iat)
         pqy =alp*xyz(2,iat)
         pqz =alp*xyz(3,iat)
         ab  =alp**2 
         eab =2d0*alp
         eabo=1d0/eab
         fn  =(pi2*alp)**1.5d0
         f   =fn*tpi*eabo 
!        s00 = s00 + v1(i)*fn*(pi*eabo)**1.5d0  ! norm ok
         pq1 =2d0*pqx*eabo
         pq2 =2d0*pqy*eabo
         pq3 =2d0*pqz*eabo
         v00 =0  
         do k=1,nat
            dx=pq1-xyz(1,k)
            dy=pq2-xyz(2,k)
            dz=pq3-xyz(3,k)
            cpsq=dx*dx+dy*dy+dz*dz
            v00=v00+atom_par(5,at(k))*f02(eab*cpsq)
         enddo   
         v1(i) = v00 * f
      enddo

!     write(*,*) 'auxnorm q',s00

end
