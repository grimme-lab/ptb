!ccccccccccccccccccccccccccccccccccccccccccccc
! compute F P S - S P F
! tr (FDS-SDF) =! min
! (P.Pulay, J.Comp.Chem. 556 (1982)
!ccccccccccccccccccccccccccccccccccccccccccccc

subroutine fpcom(ndim,S,F,P)
      implicit none
      integer ndim
      real*8 S(ndim*(ndim+1)/2)
      real*8 F(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)

      real*8,allocatable ::fps(:),spf(:),tmp(:)
      real*8 norm, tr
      integer i,j,ij
              
      allocate(fps(ndim*(ndim+1)/2),spf(ndim*(ndim+1)/2),tmp(ndim*(ndim+1)/2))

      fps = 0
      spf = 0
!     call prmat(6,S,ndim,0,'S')
!     call prmat(6,F,ndim,0,'F')
!     call prmat(6,P,ndim,0,'P')

      call ml3dgm(fps,F,P,S,ndim)
      call ml3dgm(spf,S,P,F,ndim)

      norm=0
      ij = 0
      do i=1,ndim
         do j=1,i
            ij = ij + 1
            tmp(ij) = fps(ij)-spf(ij)
            norm = norm + tmp(ij)**2
         enddo
      enddo

      ij = 0
      tr = 0
      do i=1,ndim
         ij = ij + i
         tr = tr + tmp(ij)
      enddo

!     call prmat(6,tmp,ndim,0,'FPS-SPF')
      write(*,*) 'TR and norm of FPS - SPF ',tr,sqrt(norm)

end
