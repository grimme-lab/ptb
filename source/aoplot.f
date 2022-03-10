      real*8 dex2
      real*8 cont(10),expo(10)
      real*8 dex(10)
      character*1 atmp

      PI=4.D0*ATAN(1.D0)
      DO I=-1,10
         DEX(I)=DEX2(I)
      ENDDO

      do ns = 1, 7

      iam = 0
      read(*,*) n, atmp
      if(atmp.eq.'p') iam = 1
      if(atmp.eq.'d') iam = 2
      if(atmp.eq.'f') iam = 3
      write(*,*) iam, n
      tnorm=0
      do i=1,n
         read(*,*) expo(i),cont(i)
         XNORM=(2.D0*EXPO(i)/PI)**0.75D0*(4.D0*EXPO(i))**(IAM/2.D0)/
     .          SQRT(DEX(2*IAM-1))
         tnorm = tnorm + cont(i)
         cont(i)=cont(i)*xnorm
      enddo
      cont = cont /tnorm

      x=0

      write(42,*)
 10   val = 0
      do i=1,n
         val = val + x**dble(iam)*cont(i)*exp(-expo(i)*x*x)
      enddo   
      x = x + 0.02
      write(42,*) x,val
      write(* ,*) x,val
      if(x.lt.3) goto 10

      enddo

      end

      DOUBLE PRECISION FUNCTION DEX2(M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(M .LT. 2) THEN
         DEX2=1
      ELSE
         DEX2=1
         DO 10 I=1,M,2
   10    DEX2=DEX2*I
      ENDIF
      RETURN
      END

