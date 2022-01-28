!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dtrf3(s,li,lj)
      implicit none
      real*8, intent(inout) :: s(6,6)
      integer,intent(in)    :: li,lj
      interface
         SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,
     .                    B,LDB,BETA,C,LDC)
         INTENT(IN) ALPHA,BETA,K,LDA,LDB,LDC,M,N,TRANSA,TRANSB,A,B
         INTENT(INOUT) C
*     .. Scalar Arguments ..
         DOUBLE PRECISION ALPHA,BETA
         INTEGER K,LDA,LDB,LDC,M,N
         CHARACTER TRANSA,TRANSB
*     .. Array Arguments ..
         DOUBLE PRECISION A(lda,*),B(ldb,*),C(ldc,*)
         END SUBROUTINE DGEMM
      end interface
c CAO-AO transformation
      real*8 :: trafo(6,6)
      parameter (trafo = reshape((/ ! copied from scf.f, simplyfied
c --- dS
     . sqrt(1.0d0/5.0d0),
     . sqrt(1.0d0/5.0d0),
     . sqrt(1.0d0/5.0d0),
     . 0.0d0,0.0d0,0.0d0,
c --- dx²-y²
     . 0.5d0*sqrt(3.0d0),
     .-0.5d0*sqrt(3.0d0),
     . 0.0d0,0.0d0,0.0d0,0.0d0,
c --- dz²                        %SG: change sign for ORCA convention
     .-0.5d0,-0.5d0,1.0d0,
     . 0.0d0, 0.0d0,0.0d0,
c --- rest
     . 0.0d0,0.0d0,0.0d0,sqrt(3.0d0),0.0d0,0.0d0,
     . 0.0d0,0.0d0,0.0d0,0.0d0,sqrt(3.0d0),0.0d0,
     . 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,sqrt(3.0d0) /), shape(trafo)))

      real*8 s2(6,6),sspher,dum(6,6)
      integer ii,jj,m,n

!     transformation not needed for pure s/p overlap -> do nothing
!     if (li.lt.2.and.lj.lt.2) return

! --- means li.ge.2.or.lj.ge.2, so one of them is a d-shell
!     assuming its on jat ... a wild guess
      select case(li)
      case(0) ! s-d
         do jj=1,6
         sspher=0
         do m=1,6
            sspher=sspher+trafo(m,jj)*s(m,1)
         enddo
         s2(jj,1)=sspher
         enddo
         s(1:5,1) = s2(2:6,1)
         return
      case(1) ! p-d
         do ii=1,3
         do jj=1,6
         sspher=0
         do m=1,6
            sspher=sspher+trafo(m,jj)*s(m,ii)
         enddo
         s2(jj,ii)=sspher
         enddo
         s(1:5,ii) = s2(2:6,ii)
         enddo
         return
      end select
!     wasn't there, then try iat ...
      select case(lj)
      case(0) ! d-s
         do jj=1,6
         sspher=0
         do m=1,6
            sspher=sspher+trafo(m,jj)*s(1,m)
         enddo
         s2(1,jj)=sspher
         enddo
         s(1,1:5) = s2(1,2:6)
         return
      case(1) ! d-p
         do ii=1,3
         do jj=1,6
         sspher=0
         do m=1,6
            sspher=sspher+trafo(m,jj)*s(ii,m)
         enddo
         s2(ii,jj)=sspher
         enddo
         s(ii,1:5) = s2(ii,2:6)
         enddo
         return
      end select
!     if not returned up to here -> d-d
! CB: transposing s in first dgemm is important for integrals other than S
      CALL DGEMM('T','N',6,6,6,1.D0,s,6,trafo,6,0.D0,dum,6)
      CALL DGEMM('T','N',6,6,6,1.D0,dum,6,trafo,6,0.D0,s2,6)
      s(1:5,1:5) = s2(2:6,2:6)
      return

      end

