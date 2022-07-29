!! ------------------------------------------------------------------------
!  compute restricted OS atomic spin populations on qs from MOs U
!  and Fermi occ.
!! ------------------------------------------------------------------------
subroutine spinden(n,ndim,nel,nopen,homo,et,S,e,U,psh)
      implicit none
      integer n,ndim,nsh,nel,nopen,homo
      real*8 et    
      real*8 psh(10,n)
      real*8 S(ndim*(ndim+1)/2)
      real*8 e(ndim)
      real*8 U(ndim,ndim)

      integer i,j
      integer ihomoa,ihomob
      real*8 nfoda,nfodb,ga,gb,efa,efb
      real*8 ,allocatable ::focca(:), foccb(:), qa(:,:), qb(:,:)
      real*8 ,allocatable ::sdum(:,:)            
      real*8 ,allocatable ::P(:)            

      allocate (qa(10,n),qb(10,n),sdum(ndim,ndim),focca(ndim),foccb(ndim),P(ndim*(ndim+1)/2))

      ga=0
      gb=0
! Fermi smearing                                          
!     convert restricted occ first to alpha/beta             
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

      call dmat(ndim,focca,U,sdum)
      call packsym(ndim,sdum,P)
      call mpop2(n,ndim,P,S,qa)

      call dmat(ndim,foccb,U,sdum)
      call packsym(ndim,sdum,P)
      call mpop2(n,ndim,P,S,qb)

      psh = qa - qb

      end

