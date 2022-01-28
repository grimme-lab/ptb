      subroutine CalcGammaF(val,x,m)

!     computation of function Fm(x) described above
!     for a batch of m values: val(i)=Fi(x), i=0,m
!     relative accuracy better than 5.e-15
!     RA 1/04   modified 6/05 : better accuracy
!     and better speed by initializing array ai

      implicit real(8)(a-h,o-z)
      real*8 val(0:m)

      data pi/3.1415926535897932384626433832795029d0/
      data switch,accur/18.d0,1.d-8/
      data a1,a2,ah/1.d0,2.d0,0.5d0/
      real(8),dimension(75) ::     ! modified 6/05
     *           ai=(/(1.d0/(2.d0*i-1.d0),i=1,75)/)

      if (x.gt.switch) go to 305
!     get Fm
      xx=x+x
      t=ai(1+m)
      s=t
      do i=1,(75-m-1)
        t=t*xx*ai(1+m+i)
        s=s+t
        if (abs(t).lt.s*accur) go to 225
      enddo
      write (*,*) 'gammav fails for x= ',x,'  m=',m
      stop 'no convergence in <gammav>'
 225  continue
!     now Fi(x), i=m,..,0 by downward recursion
      ee=exp(-x)
      s=s*ee
      do i=m,1,-1
        val(i)=s
        s=(xx*s+ee)*ai(i)
      enddo
      val(0)=s
      return

 305  continue
!     get F0 from asymptotic expansion
      xi=a1/x
      s = sqrt(pi*xi)*ah
      ee=0.d0
      if (x.lt.37.0d0+m*2) ee = exp(-x) ! modified 6/05
      xdi = -ah*xi
      t=ee*xdi
      an=-a1
      do i=1,50
        an=an+a2
        hh=t*xdi*an
        if(abs(hh).gt.abs(t)) stop 'error in <gammav>'
        s=s+t
        if (abs(t).lt.accur*s) go to 425
        t=hh
      enddo
      write (*,*) 'gammav fails for x= ',x,'  m=',m
      stop 'no convergence in <gammaf>'
 425  val(0) = s
!     now Fi(x) by upward recursion: i=1,..,m
      n=1
      xd=xi*ah
      an=a1
      do i=1,m
        val(i) = (an*val(i-1)-ee)*xd
        an=an+a2
      enddo
      return

      end
