
***********************************************************************
*                                                                     *
*  make the loewdin orthogonalization matrix x = u' s1/2 u            *
*                                                                     *
***********************************************************************

      subroutine makel(nao, s, x)
      implicit real*8 (a-h,o-z)
      dimension s(*)
      dimension x(nao,nao)
      real*8, allocatable ::aux(:),vecs(:,:),e(:),cc(:,:)

      lwork  = 1 + 6*nao + 2*nao**2
      allocate (vecs(nao,nao),e(nao),aux(lwork),cc(nao,nao))

      k=0
      do i=1,nao
         do j=1,i
            k=k+1
            vecs(j,i)=s(k)
            vecs(i,j)=s(k)
         enddo
      enddo

      call dsyev ('V','U',nao,vecs,nao,e,aux,lwork,info)
      
c     call dHQRII(s,nao,nao,e,vecs)

      do i=1,nao
         if(e(i).lt.0) stop 'sorry, must stop in S^1/2!'
         e(i)=dsqrt(e(i))
      enddo

      do m=1,nao 
         do i=1,nao
         x (i,m)=     vecs(i,m)
         cc(i,m)=e(m)*vecs(i,m)
         enddo
      enddo

      call dgemm('N','T',nao,nao,nao,1.0d0,x,
     .                   nao,cc,nao,0.0d0,vecs,nao)

      x = vecs
      deallocate(e,aux,cc,vecs)

      return
      end

