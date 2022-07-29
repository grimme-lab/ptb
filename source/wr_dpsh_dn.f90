      implicit none
      real*8 psh1(10)
      real*8 psh2(10)
      real*8 psh3(10)
      real*8 n1,n2,x,xmax,xmin,nn,thr
      real*8 q(1000)                    
      integer i,n,ns,norm

      xmax = -1000.
      xmin =  1000.
      norm = 0
      nn   = 0
      psh3 = 0
      thr  = 3d-3
      open(unit=1,file='.datadpdn')
 10   read(1,*,end=100) n,ns,psh1(1:ns)
      read(1,*)         n,ns,psh2(1:ns)
      n1=sum(psh1(1:ns))
      n2=sum(psh2(1:ns))
      if(abs(n1-n2).lt.n1*thr.or.abs(n1-n2).gt.5d0)then
         write(*,*) norm,'dq too small:', n1,n2, ' skipping'
         goto 10
      endif
      norm = norm + 1
      psh3(1:ns)=psh3(1:ns)+(psh2(1:ns)-psh1(1:ns))/(n2-n1)
      nn   = nn   + n1
      q(norm) = n1
      do i=1,ns
         x=(psh2(i)-psh1(i))/(n2-n1)
         if(x.gt.xmax) xmax = x
         if(x.lt.xmin) xmin = x
      enddo
      goto 10
100   continue
      psh3(1:ns)=psh3(1:ns)/dble(norm)
      if(abs(1d0-sum(psh3(1:ns))).gt.1d-6) stop 'norm error'
      write(*,*) '# molecules ', norm
      nn = nn / dble(norm)
      q(1:norm)=q(1:norm)-nn
      write(*,*) 'average N, max. change   ', nn, maxval(abs(q(1:norm)))
      write(*,*) 'min/max. dp/dq  ', xmin,xmax     
      do i=1,ns
        write(*,'(6x,''dpshdq('',i2,'','',i2,'')='',F16.12)')i,n,psh3(i)
      enddo
      end
