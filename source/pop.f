! atoms and shells

      subroutine mpop(n,ndim,P,S,q,qsh)
      use bascom
      implicit none
      integer n,ndim
      real*8 q(n),qsh(10,n)
      real*8 P(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)

      integer i,j,ij,ii,jj,ish,jsh
      real*8 ps

      qsh= 0
      q  = 0 
      ij = 0
      do i=1,ndim
         ii=aoat(i)
         ish=shell2ao(i)
         do j=1,i-1
            ij=ij+1
            jj=aoat(j)
            jsh=shell2ao(j)
            ps=p(ij)*s(ij)
            q(ii)=q(ii)+ps
            q(jj)=q(jj)+ps
            qsh(ish,ii)=qsh(ish,ii)+ps
            qsh(jsh,jj)=qsh(jsh,jj)+ps
         enddo
         ij=ij+1
         ps=p(ij)*s(ij)
         q(ii)=q(ii)+ps    
         qsh(ish,ii)=qsh(ish,ii)+ps    
      enddo

      end

! just atoms

      subroutine mpop3(n,ndim,P,S,q)
      use bascom
      implicit none
      integer n,ndim
      real*8 q(n)
      real*8 P(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)

      integer i,j,ij,ii,jj
      real*8 ps

      q  = 0 
      ij = 0
      do i=1,ndim
         ii=aoat(i)
         do j=1,i-1
            ij=ij+1
            jj=aoat(j)
            ps=p(ij)*s(ij)
            q(ii)=q(ii)+ps
            q(jj)=q(jj)+ps
         enddo
         ij=ij+1
         ps=p(ij)*s(ij)
         q(ii)=q(ii)+ps    
      enddo

      end



!--------------------------------------------------

      subroutine mpopao(ndim,P,S,qao)
      use bascom
      implicit none
      integer ndim
      real*8 qao(ndim)
      real*8 P(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)

      integer i,j,ij,ii,jj,ish,jsh
      real*8 ps

      qao= 0 
      ij = 0
      do i=1,ndim
         do j=1,i-1
            ij=ij+1
            ps=p(ij)*s(ij)
            qao(i)=qao(i)+ps
            qao(j)=qao(j)+ps
         enddo
         ij=ij+1
         ps=p(ij)*s(ij)
         qao(i)=qao(i)+ps    
      enddo

      end

cccccccccccccccccccccccccccccccccccccccccccccc      
c            Wiberg BOs      
cccccccccccccccccccccccccccccccccccccccccccccc      

      subroutine wiberg(n,ndim,at,rab,P,S,wbo)
      use bascom
      implicit none
      integer n,ndim,at(n)
      real*8 rab(n*(n+1)/2)
      real*8 P(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)
      real*8 wbo(n,n)

      real*8,allocatable ::Ptmp(:,:)
      real*8,allocatable ::si  (:,:)
      real*8,allocatable ::pi  (:,:)
      real*8 xsum
      integer i,j,k,l,m
      integer llao(4)
      data llao /1,3,5,7 /
      integer aose(2,n)

      allocate(Ptmp(ndim,ndim),pi(ndim,ndim),si(ndim,ndim))

      call blowsym(ndim,P,pi)
      call blowsym(ndim,S,si)
      call DGEMM('N','N',ndim,ndim,ndim,1.0d0,pi,
     .                   ndim,si,ndim,0.0d0,Ptmp,ndim)

      m = 1     
      do i=1,n 
         aose(1,i)=m 
         do k=1,bas_nsh(at(i))
            m= m + llao(bas_lsh(k,at(i))+1)
         enddo
         aose(2,i)=m-1
      enddo

      wbo=0
      l  =0
      do i=1,n 
         do j=1,i-1
         l = l + 1
         xsum=0.0
         if(rab(l).lt.20.0)then
            do k=aose(1,i),aose(2,i)     ! AOs on atom i
               do m=aose(1,j),aose(2,j)  ! AOs on atom j
                  xsum=xsum+Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wbo(i,j)=xsum
         wbo(j,i)=xsum
         enddo
         l = l + 1
      enddo

      end

cccccccccccccccccccccccccccccccccccccccccccccc      
c      print Wiberg BOs      
cccccccccccccccccccccccccccccccccccccccccccccc      

      subroutine prwbo(n,at,wbo)
      use bascom
      implicit none
      logical pr
      integer n,at(n)
      real*8 wbo(n,n)

      real*8,allocatable ::wb  (:,:)
      real*8 xsum
      integer i,j,k,l,m,ibmax,imem(n),lin
      character*2 asym

      allocate(wb(n,n))

      wb = wbo

      write(*,*)'largest Wiberg (=(PS)*(SP)) bond orders for each atom'
      write(*,*)'          total WBO             WBO to atom > 0.05 ...'
      do i=1,n
         do j=1,n
            imem(j)=j
         enddo
         call wibsort(n,i,imem,wb)
         ibmax=0
         xsum =0
         do j=1,n
            if(wb(j,i).gt.0.05)ibmax=j
            xsum=xsum+wb(j,i)
         enddo
         write(*,'(i6,a4,1x,f6.3,9(4x,a2,i4,f6.3))')
     .   i,asym(at(i)),xsum,
     .   (asym(at(imem(j))),imem(j),wb(j,i),j=1,ibmax)
      enddo

      end

cccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE wibsort(ncent,imo,imem,qmo)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension qmo(ncent,ncent)
      dimension imem(ncent)

      do 140   ii = 2,ncent
         i = ii - 1
         k = i
         pp= qmo(i,imo)
         do 120   j = ii, ncent
            if (qmo(j,imo) .lt. pp) go to 120
            k = j
            pp=qmo(j,imo)
  120    continue
         if (k .eq. i) go to 140
         qmo(k,imo) = qmo(i,imo)
         qmo(i,imo) = pp

         ihilf=imem(i)
         imem(i)=imem(k)
         imem(k)=ihilf
  140 continue

      end

