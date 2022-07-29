!ccccccccccccccccccccccccccccccccccccccccccccc
! compute F (P S) - S (P F)
! (P.Pulay, J.Comp.Chem. 556 (1982)
!ccccccccccccccccccccccccccccccccccccccccccccc

subroutine fpcom(ndim,S,F,P)
      implicit none
      integer ndim
      real*8 S(ndim*(ndim+1)/2)
      real*8 F(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)

      real*8,allocatable ::PF(:,:), FP(:,:)
      real*8,allocatable ::Fdum(:,:),Pdum(:,:), Sdum(:,:)
      real*8,allocatable ::SP(:,:),PS(:,:)
      real*8 norm,maxv,fppf
      integer i,j,ij
              
      allocate(FP(ndim,ndim),PF(ndim,ndim),Fdum(ndim,ndim),Pdum(ndim,ndim))
      allocate(Sdum(ndim,ndim),PS(ndim,ndim))

      call blowsym(ndim,S,Sdum)
      call blowsym(ndim,F,Fdum)
      call blowsym(ndim,P,Pdum)

      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,Pdum,ndim,Fdum,ndim,0.0d0,PF,ndim)
      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,Pdum,ndim,Sdum,ndim,0.0d0,PS,ndim)
      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,Fdum,ndim,PS  ,ndim,0.0d0,FP  ,ndim)
      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,Sdum,ndim,PF  ,ndim,0.0d0,Pdum,ndim)

      norm=0
      maxv=0
      do i=1,ndim
         do j=1,i
            fppf = FP(j,i)-pdum(j,i)
            if(abs(fppf).gt.maxv) maxv=abs(fppf)
            norm = norm + fppf**2
         enddo
      enddo

!     call prmat(6,fppf,ndim,0,'FP-PF')
      write(*,'('' F(3)P(2)S - SP(2)F(3)     : '',2E16.8)') sqrt(norm)/dble(ndim*(ndim+1)/2),maxv

end
