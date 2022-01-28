      SUBROUTINE rdmat(IUOUT,X,M)
      REAL*8 X(M*(M+1)/2)
      REAL*8 R(M,M)
      CHARACTER*128 line
      integer i,j,ij
      N=M
      NKPB=6
      IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      read (IUOUT,'(a)') line    
      DO 90 J=1,N
      read (IUOUT,*)idum,(R(J,II),II=K1,K2)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      read (IUOUT,'(a)') line    
      DO 110 J=1,N
      read (IUOUT,*)idum,(R(J,II),II=K1,K2)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 CONTINUE
      ij = 0
      do i=1,m
         do j=1,i
            ij = ij + 1
            X(ij)= 0.5*(R(j,i)+R(i,j))
         enddo
      enddo
      RETURN
      END

      SUBROUTINE rdmat2(IUOUT,R,M)
      REAL*8 R(M,M)
      CHARACTER*128 line
      CHARACTER*17  a17 
      integer i,j,ij
      N=M
      NKPB=6
      IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      read (IUOUT,'(a)') line    
      read (IUOUT,'(a)') line    
      read (IUOUT,'(a)') line    
      read (IUOUT,'(a)') line    
      DO 90 J=1,N
      read (IUOUT,'(a16,6f10.6)')a17,(R(J,II),II=K1,K2)
!     write(*    ,'(a16,6f10.6)')a17,(R(J,II),II=K1,K2)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      read (IUOUT,'(a)') line    
      read (IUOUT,'(a)') line    
      read (IUOUT,'(a)') line    
      read (IUOUT,'(a)') line    
      DO 110 J=1,N
      read (IUOUT,'(a16,6f10.6)')a17,(R(J,II),II=K1,K2)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 CONTINUE
      RETURN
      END
