!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mixed Mulliken/Loewdin population analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mlpop(n, nao, p, s, q, qsh)
      use bascom
      use parcom, only: mull_loew14
      implicit none             
      integer n, nao
      real*8 s (nao*(nao+1)/2)   ! overlap
      real*8 p (nao*(nao+1)/2)   ! density, unchanged 
      real*8 q(n)                ! aomic pops       
      real*8 qsh(10,n)           ! shell pops

      real*8, allocatable :: aux(:),vecs(:,:),e(:),s1(:,:),
     .                       s2(:,:),e1(:),e2(:),cc(:,:),x(:,:)
      integer lwork,info

      integer i,m,ii,ish
      real*8  ps

      allocate (
     .          e(nao),e1(nao),e2(nao),
     .          vecs(nao,nao),cc(nao,nao),x(nao,nao),
     .          s1(nao,nao),s2(nao,nao))

      call blowsym(nao,S,vecs)

      allocate (aux(1))
      lwork=-1
      call dsyev ('V','U',nao,vecs,nao,e,aux,lwork,info)
      lwork=idint(aux(1))
      deallocate(aux)
      allocate (aux(lwork))
      call dsyev ('V','U',nao,vecs,nao,e,aux,lwork,info)
      
      if(e(1).le.0.or.info.ne.0) stop 'sorry, must stop in mlpop!'

      do i=1,nao
         e1(i)=e(i)**(1d0-mull_loew14)
         e2(i)=e(i)**mull_loew14 ! 1/4= more L, 1/8 more M
      enddo

      do m=1,nao 
         do i=1,nao
            cc(i,m)=e1(m)*vecs(i,m)  
         enddo
      enddo
      call dgemm('N','T',nao,nao,nao,1.0d0,vecs,     ! vecs = S^1/2
     .                   nao,cc,nao,0.0d0,s1,nao)
      do m=1,nao 
         do i=1,nao
            cc(i,m)=e2(m)*vecs(i,m)  
         enddo
      enddo
      call dgemm('N','T',nao,nao,nao,1.0d0,vecs,     ! vecs = S^1/2
     .                   nao,cc,nao,0.0d0,s2,nao)

      call blowsym(nao,P,x)

      call dgemm('N','N',nao,nao,nao,1.0d0,x,      ! cc=P*S^1/2
     .                   nao,s1,nao,0.0d0,cc,nao)

      call dgemm('N','N',nao,nao,nao,1.0d0,s2,     ! PL=S^1/2*cc  
     .                   nao,cc,nao,0.0d0,x,nao)

!     call prmat(6,x,nao,nao,'PL')

      qsh= 0
      q  = 0 
      do i=1,nao 
         ii=aoat(i)
         ish=shell2ao(i)
         ps=x(i,i) 
         q(ii)=q(ii)+ps
         qsh(ish,ii)=qsh(ish,ii)+ps
      enddo

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mixed Mulliken/Loewdin population analysis
! stepwise 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mlpop1(nao, s, s34, s14)
      use parcom, only: mull_loew14
      implicit none             
      integer nao
      real*8 s (nao*(nao+1)/2)   ! overlap
!     real*8   X(nao,nao)        ! diag trafo matrix
      real*8 s34(nao,nao)        ! s^3/4
      real*8 s14(nao,nao)        ! s^1/4

      real*8, allocatable :: aux(:),e(:),e1(:),e2(:),cc(:,:),x(:,:)
      integer lwork,info

      integer i,j,m
      real*8  xsum

      allocate(e(nao),e1(nao),e2(nao),cc(nao,nao),x(nao,nao))

      call blowsym(nao,S,X)

      allocate(aux(1))
      lwork=-1
      call dsyev ('V','U',nao,X,nao,e,aux,lwork,info)
      lwork=idint(aux(1))
      deallocate(aux)
      allocate (aux(lwork))
      call dsyev ('V','U',nao,X,nao,e,aux,lwork,info)

      do i=1,nao
         e1(i)=e(i)**(1d0-mull_loew14)
         e2(i)=e(i)**mull_loew14 ! 1/4= more L, 1/8 more M
      enddo

      do m=1,nao 
         do i=1,nao
            cc(i,m)=X(i,m) * e1(m)
         enddo
      enddo
      call dgemm('N','T',nao,nao,nao,1.0d0,X,     ! = S^3/4
     .                   nao,cc,nao, 0.0d0,s34,nao)
      do m=1,nao 
         do i=1,nao
            cc(i,m)=X(i,m) * e2(m)
         enddo
      enddo
      call dgemm('N','T',nao,nao,nao,1.0d0,X,     ! = S^1/4
     .                   nao,cc,nao, 0.0d0,s14,nao)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mixed Mulliken/Loewdin population analysis
! stepwise 1, real*4 version 
! ML trafo output also in real*4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mlpop14(nao, s, s34, s14)
      use parcom, only: mull_loew14
      implicit none             
      integer nao
      real*8 s (nao*(nao+1)/2)   ! overlap
      real*4 s34(nao,nao)        ! s^3/4
      real*4 s14(nao,nao)        ! s^1/4

      real*4, allocatable :: aux(:),e(:),x(:,:)
      real*4, allocatable :: e1(:),e2(:),cc(:,:)
      integer lwork,info

      integer i,m

      allocate(e(nao),e1(nao),e2(nao),x(nao,nao),cc(nao,nao))

      call blowsym84(nao,S,X)

      allocate(aux(1))
      lwork=-1
      call ssyev ('V','U',nao,X,nao,e,aux,lwork,info)
      lwork=int(aux(1))
      deallocate(aux)
      allocate (aux(lwork))
      call ssyev ('V','U',nao,X,nao,e,aux,lwork,info)

      do i=1,nao
         e1(i)=e(i)**(1e0-mull_loew14)       
         e2(i)=e(i)**mull_loew14 ! 1/4= more L, 1/8 more M
      enddo

      do m=1,nao 
         do i=1,nao
            cc(i,m)=X(i,m) * e1(m)
         enddo
      enddo
      call sgemm('N','T',nao,nao,nao,1.0e0,X,     ! = S^3/4
     .                   nao,cc,nao, 0.0e0,s34,nao)

      do m=1,nao 
         do i=1,nao
            cc(i,m)=X(i,m) * e2(m)
         enddo
      enddo
      call sgemm('N','T',nao,nao,nao,1.0e0,X,     ! = S^1/4
     .                   nao,cc,nao, 0.0e0,s14,nao)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mixed Mulliken/Loewdin population analysis
! stepwise 2 for shells and atoms
! ML input in real*4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mlpop2(n, nao, p, s1, s2, q, qsh)
      use bascom
      implicit none             
      integer n, nao
      real*8 p (nao*(nao+1)/2)   ! density, unchanged 
      real*4 s1(nao,nao)         ! s^3/4
      real*4 s2(nao,nao)         ! s^1/4
      real*8 q(n)                ! atomic pops       
      real*8 qsh(10,n)           ! shell pops

      real*4, allocatable :: cc(:,:),x(:,:)

      integer i,ii,ish
      real*8  ps

      allocate (cc(nao,nao),x(nao,nao)) 

      call blowsym84(nao,P,x)

!     call sgemm('N','N',nao,nao,nao,1.0e0,x,      ! cc=P*S^3/4
!    .                   nao,s1,nao,0.0e0,cc,nao)
      call ssymm('L','L',nao,nao,1e0,x,nao,S1,nao,0e0,cc,nao)   

      call sgemm('N','N',nao,nao,nao,1.0e0,s2,     ! PL=S^1/4*cc  
     .                   nao,cc,nao,0.0e0,x,nao)

      qsh= 0
      q  = 0 
      do i=1,nao 
         ps=x(i,i) 
         ii=aoat(i)
         q(ii)=q(ii)+ps
         ish=shell2ao(i)
         qsh(ish,ii)=qsh(ish,ii)+ps
      enddo

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mixed Mulliken/Loewdin population analysis
! stepwise 2 for atoms onl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mlpop3(n, nao, p, s1, s2, q)
      use bascom
      implicit none             
      integer n, nao
      real*8 p (nao*(nao+1)/2)   ! density, unchanged 
      real*8 s1(nao,nao)         ! s^3/4
      real*8 s2(nao,nao)         ! s^1/4
      real*8 q(n)                ! atomic pops       

      real*8, allocatable :: cc(:,:),x(:,:)

      integer i,ii
      real*8  ps

      allocate (cc(nao,nao),x(nao,nao)) 

      call blowsym(nao,P,x)

      call dgemm('N','N',nao,nao,nao,1.0d0,x,      ! cc=P*S^3/4
     .                   nao,s1,nao,0.0d0,cc,nao)

      call dgemm('N','N',nao,nao,nao,1.0d0,s2,     ! PL=S^1/4*cc  
     .                   nao,cc,nao,0.0d0,x,nao)

      q  = 0 
      do i=1,nao 
         ii=aoat(i)
         ps=x(i,i) 
         q(ii)=q(ii)+ps
      enddo

      end
