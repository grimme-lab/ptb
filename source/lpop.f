!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loewdin population analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lpop(n, nao, p, s, q, qsh)
      use bascom
      implicit none             
      integer n, nao
      real*8 s (nao*(nao+1)/2)   ! overlap
      real*8 p (nao*(nao+1)/2)   ! density, on output P_Loewdin basis
      real*8 q(n)                ! aomic
      real*8 qsh(10,n)

      real*8, allocatable ::aux(:),vecs(:,:),e(:),cc(:,:),x(:,:)
      integer lwork,info

      integer i,j,k,m,ii,jj,ish,jsh,lin
      real*8  ps

      lwork  = 1 + 6*nao + 2*nao**2
      allocate (vecs(nao,nao),e(nao),aux(lwork),cc(nao,nao),x(nao,nao))

      call blowsym(nao,S,vecs)

      call dsyev ('V','U',nao,vecs,nao,e,aux,lwork,info)
      
      if(e(1).le.0.or.info.ne.0) stop 'sorry, must stop in S^1/2!'

      do i=1,nao
         e(i)=dsqrt(e(i))
      enddo

      do m=1,nao 
         do i=1,nao
         x (i,m)=     vecs(i,m)
         cc(i,m)=e(m)*vecs(i,m)  
         enddo
      enddo

      call dgemm('N','T',nao,nao,nao,1.0d0,x,        ! vecs = S^1/2
     .                   nao,cc,nao,0.0d0,vecs,nao)

      call blowsym(nao,P,x)

      call dgemm('N','N',nao,nao,nao,1.0d0,x,        ! cc=P*S^1/2
     .                   nao,vecs,nao,0.0d0,cc,nao)

      call dgemm('N','N',nao,nao,nao,1.0d0,vecs,     ! PL=S^1/2*cc  
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

