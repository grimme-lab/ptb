!ccccccccccccccccccccccccccccccccccccccccccccc
! compute F (P S) - S (P F)
! (P.Pulay, J.Comp.Chem. 556 (1982)
!ccccccccccccccccccccccccccccccccccccccccccccc

subroutine fpcom(ndim,S,F,P)
      use gtb_la, only : la_gemm
      implicit none
      integer ndim
      real*8 S(ndim*(ndim+1)/2)
      real*8 F(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)

      real*8,allocatable ::PF(:,:), FP(:,:), fppf(:)
      real*8,allocatable ::Fdum(:,:),Pdum(:,:), Sdum(:,:)
      real*8,allocatable ::SP(:,:),PS(:,:)
      real*8 norm,tr
      integer i,j,ij
              
      allocate(FP(ndim,ndim),PF(ndim,ndim),fppf(ndim*(ndim+1)/2),Fdum(ndim,ndim),Pdum(ndim,ndim))
      allocate(Sdum(ndim,ndim),PS(ndim,ndim))

      call blowsym(ndim,S,Sdum)
      call blowsym(ndim,F,Fdum)
      call blowsym(ndim,P,Pdum)

      call la_gemm('N','N',ndim,ndim,ndim,1.0d0,Pdum,ndim,Fdum,ndim,0.0d0,PF,ndim)
      call la_gemm('N','N',ndim,ndim,ndim,1.0d0,Pdum,ndim,Sdum,ndim,0.0d0,PS,ndim)
      call la_gemm('N','N',ndim,ndim,ndim,1.0d0,Fdum,ndim,PS  ,ndim,0.0d0,FP  ,ndim)
      call la_gemm('N','N',ndim,ndim,ndim,1.0d0,Sdum,ndim,PF  ,ndim,0.0d0,Pdum,ndim)

      tr  =0
      norm=0
      ij = 0
      do i=1,ndim
         do j=1,i-1
            ij = ij + 1
            fppf(ij)=FP(j,i)-pdum(j,i)
            norm = norm + fppf(ij)**2
         enddo
         ij = ij + 1
         fppf(ij)=FP(i,i)-pdum(i,i)
         tr = tr + fppf(ij)
      enddo

!     call prmat(6,fppf,ndim,0,'FP-PF')
      write(*,*) 'TR and norm of FPS - SPF ',tr,sqrt(norm)

end
