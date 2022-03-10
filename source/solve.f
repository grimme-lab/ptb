
      subroutine solve(ndim,nel,nopen,homo,et,focc,H,S,P,e,ge,fail)
      use gtb_la, only : la_sygvx
      implicit none
      integer ndim,nel,nopen,homo
      real*8 H(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)
      real*8 e(ndim)
      real*8 focc(ndim)
      real*8 ge        
      real*8 et        
      logical fail

      integer i,j,info,lwork,ij,iu
      integer ihomoa,ihomob
      real*8 nfoda,nfodb,ga,gb,efa,efb
      integer,allocatable ::iwork(:),ifail(:)
      real*8 ,allocatable ::D(:,:),hdum(:,:),sdum(:,:),work(:)
      real*8 ,allocatable ::focca(:), foccb(:)

      fail =.false.
      allocate (D(ndim,ndim),hdum(ndim,ndim),sdum(ndim,ndim),
     .          focca(ndim),foccb(ndim))

      call blowsym(ndim,H,hdum)
      call blowsym(ndim,S,sdum)

! for a large basis, taking only the occ. eigenvalues is faster than a full diag
      iu=min(homo+5,ndim)
      allocate(iwork(5*ndim),ifail(ndim),work(1))
      call la_sygvx(1,'V','I','U',ndim, hdum, ndim, sdum, ndim, ga, gb, 
     .              1, IU, 1d-7, ij, e, D, ndim, WORK, -1   , IWORK, 
     .              IFAIL, INFO )
      lwork=idint(work(1))
      deallocate(work)
      allocate(work(lwork))          
      call la_sygvx(1,'V','I','U',ndim, hdum, ndim, sdum, ndim, ga, gb, 
     .              1, IU, 1d-7, ij, e, D, ndim, WORK, LWORK, IWORK, 
     .              IFAIL, INFO )
      if(info.ne.0) fail=.true.
      e(iu+1:ndim)=100d0
      deallocate(hdum)

      ga=0
      gb=0
c Fermi smearing                                          
c     convert restricted occ first to alpha/beta             
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0d0
         foccb=0.0d0
         ihomoa=0
         ihomob=0
      endif
      if(ihomoa+1.le.ndim) then 
         call FERMISMEAR(.false.,ndim,ihomoa,et,e,focca,nfoda,efa,ga)
      endif
      if(ihomob+1.le.ndim.and.nel.gt.1) then
         call FERMISMEAR(.false.,ndim,ihomob,et,e,foccb,nfodb,efb,gb)
      endif
      focc = focca + foccb
      ge = ga + gb

      call dmat(ndim,focc,D,sdum)
      call packsym(ndim,sdum,P)

      end

!ccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solve_spec(ndim,nel,nopen,homo,et,focc,H,X,P,e,ge,fail)
      use gtb_la, only : la_gemm, la_syevx
      implicit none
      integer ndim,nel,nopen,homo
      real*8 H(ndim*(ndim+1)/2)
      real*8 X(ndim,ndim)
      real*8 P(ndim*(ndim+1)/2)
      real*8 e(ndim)
      real*8 focc(ndim)
      real*8 ge        
      real*8 et        
      logical fail

      integer i,j,ij,info,lwork,liwork,iu
      integer ihomoa,ihomob
      real*8 nfoda,nfodb,ga,gb,efa,efb
      integer,allocatable ::iwork(:),ifail(:)
      real*8 ,allocatable ::D(:,:),hdum(:,:),sdum(:,:),work(:)
      real*8 ,allocatable ::focca(:), foccb(:)

      allocate (D(ndim,ndim),hdum(ndim,ndim),focca(ndim),foccb(ndim))

      call blowsym(ndim,H,hdum)
!     go to OAO basis
      call la_gemm('N','N',ndim,ndim,ndim,1d0,hdum,ndim,
     &               X,ndim,0d0,D,ndim)
      call la_gemm('T','N',ndim,ndim,ndim,1d0,X,ndim,
     &               D,ndim,0d0,hdum,ndim)

      iu=min(homo+5,ndim)
! for a large basis, taking only the occ. eigenvalues is faster than a full diag
      lwork  = -1
      allocate(iwork(5*ndim),ifail(ndim),work(1))
      call la_syevx('V','I','U',ndim, hdum, ndim, ga, gb, 1, IU, 1d-6,
     .              ij, e, D, ndim, WORK,LWORK,IWORK,IFAIL,INFO)
      lwork = idint(work(1))
      deallocate(work)
      allocate(work(lwork))
      fail =.false.
      call la_syevx('V','I','U',ndim, hdum, ndim, ga, gb, 1, IU, 1d-6,
     .              ij, e, D, ndim, WORK,LWORK,IWORK,IFAIL,INFO)

      if(info.ne.0.or.ij.ne.iu) fail=.true.
      e(iu+1:ndim)=100d0

      ga=0
      gb=0
c Fermi smearing                                          
c     convert restricted occ first to alpha/beta             
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0d0
         foccb=0.0d0
         ihomoa=0
         ihomob=0
      endif
      if(ihomoa+1.le.ndim) then 
         call FERMISMEAR(.false.,ndim,ihomoa,et,e,focca,nfoda,efa,ga)
      endif
      if(ihomob+1.le.ndim.and.nel.gt.1) then
         call FERMISMEAR(.false.,ndim,ihomob,et,e,foccb,nfodb,efb,gb)
      endif
      focc = focca + foccb
      ge = ga + gb

!     get AO MO coeff. from OAO ones
      call la_gemm('T','N',ndim,ndim,ndim,1d0,X,ndim,
     &               D,ndim,0d0,hdum,ndim)
      call dmat(ndim,focc,hdum,D)
      call packsym(ndim,D,P)

      end

!ccccccccccccccccccccccccccccccccccccccccccccc

      subroutine blowsym(n,matin,matout)
      implicit none
      integer n
      real*8 matin (n*(n+1)/2)
      real*8 matout(n,n)      
      integer i,j,ij

      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            matout(j,i)=matin(ij)
            matout(i,j)=matin(ij)
         enddo
      enddo

      end

      subroutine blowsym84(n,matin,matout)
      implicit none
      integer n
      real*8 matin (n*(n+1)/2)
      real*4 matout(n,n)      
      integer i,j,ij

      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            matout(j,i)=matin(ij)
            matout(i,j)=matin(ij)
         enddo
      enddo

      end

      subroutine packsym(n,matin,matout)
      implicit none
      integer n
      real*8 matin (n,n)      
      real*8 matout(n*(n+1)/2)
      integer i,j,ij

      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            matout(ij)=matin(j,i)
         enddo
      enddo

      end

      subroutine packsym48(n,matin,matout)
      implicit none
      integer n
      real*4 matin (n,n)      
      real*8 matout(n*(n+1)/2)
      integer i,j,ij

      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            matout(ij)=matin(j,i)
         enddo
      enddo

      end

      subroutine addsym(n,f,matin1,matin2,matout)
      implicit none
      integer n
      real*8 f
      real*8 matin1(n*(n+1)/2)
      real*8 matin2(n*(n+1)/2)
      real*8 matout(n*(n+1)/2)
      integer i,j,ij

      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            matout(ij)=matin1(ij)+matin2(ij)*f
         enddo
      enddo

      end

      subroutine adddsym(n,f,matin2,matout)
      implicit none
      integer n
      real*8 f
      real*8 matin2(n*(n+1)/2)
      real*8 matout(n*(n+1)/2)
      integer i,j,ij

      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            matout(ij)=matout(ij)+matin2(ij)*f
         enddo
      enddo

      end

!ccccccccccccccccccccccccccccccccccccccccccccc
! density matrix
! C   : MO coefficient
! focc: occupations      
! P  dmat
!ccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dmat(ndim,focc,C,P)
      use gtb_la, only : la_gemm
      implicit none
      integer ndim
      real*8 focc(*)
      real*8 C(ndim,ndim)
      real*8 P(ndim,ndim)
      integer i,m
      real*8,allocatable ::Ptmp(:,:)
              
      allocate(Ptmp(ndim,ndim))                  
      do m=1,ndim  
         do i=1,ndim
            Ptmp(i,m)=C(i,m)*focc(m)
         enddo
      enddo
      call la_gemm('N','T',ndim,ndim,ndim,1.0d0,C,
     .               ndim,Ptmp,ndim,0.0d0,P,ndim)
      deallocate(Ptmp)

      end

      subroutine dmat4(ndim,focc,C,P)
      use gtb_la, only : la_gemm
      implicit none
      integer ndim
      real*8 focc(*)
      real*4 C(ndim,ndim)
      real*4 P(ndim,ndim)
      integer i,m
      real*4,allocatable ::Ptmp(:,:)
              
      allocate(Ptmp(ndim,ndim))                  
      do m=1,ndim  
         do i=1,ndim
            Ptmp(i,m)=C(i,m)*focc(m)
         enddo
      enddo
      call la_gemm('N','T',ndim,ndim,ndim,1.0e0,C,
     .               ndim,Ptmp,ndim,0.0e0,P,ndim)
      deallocate(Ptmp)

      end
