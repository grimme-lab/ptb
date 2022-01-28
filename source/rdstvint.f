      subroutine rdstvint
     . (fname,n,ndim,xyz,S,F,h,T,P,e,ekin,enuc,dip,pnt,alp)
      use bascom, only: aoat
      use mocom                

      implicit none
      integer n,ndim
      real*8 S(ndim*(ndim+1)/2)
      real*8 F(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)
      real*8 h(ndim*(ndim+1)/2)
      real*8 T(ndim*(ndim+1)/2)
      real*8 xyz(3,n) 
      real*8 e,ekin,enuc,dip(3),pnt(3),alp(6)

      real*8 xx(10)
      real*8 xyzrd(3),rmsd
      character*128 line
      character*(*) fname
      integer nn,i,ich,iat

      pnt = 0
      alp = 0
      alp(1)=-99

      ich=142
      open(unit=ich,file=fname)
 100  read(ich,'(a)',end=200)line
      if(index(line,'MOLECULAR ORBITALS').ne.0) then
         read(ich,'(a)',end=200)line
         call rdmat2(ich,cmo_ref,ndim)
!        call prmat(6,cmo_ref,ndim,ndim,'Cref')
!        stop
      endif
      if(index(line,'FINAL SINGLE POINT ENERGY').ne.0) then
         call readl(line,xx,nn)
         e = xx(nn)
      endif
      if(index(line,'KINETIC ENERGY INTEGRALS').ne.0) then
         call rdmat(ich,T,ndim)
      endif
      if(index(line,'OVERLAP INTEGRALS ').ne.0) then
         call rdmat(ich,S,ndim)
      endif

      if(index(line,'CORE HAMILTONIAN ').ne.0) then
         call rdmat(ich,h,ndim)
         h=h-T
      endif

      if(index(line,'FOCK MATRIX ').ne.0) then
         call rdmat(ich,F,ndim)
      endif

      if(index(line,'DENSITY').ne.0.and.index(line,'*').eq.0) then
         read(ich,'(a)',end=200)line
         call rdmat(ich,P,ndim)
      endif
      if(index(line,'ORBITAL EN').ne.0)then
         read(ich,'(a)',end=200)line
         read(ich,'(a)',end=200)line
         read(ich,'(a)',end=200)line
         do i=1,ndim
            read(ich,'(a)',end=200)line
            call readl(line,xx,nn)
            epsref(i)=xx(3)
         enddo
      endif
      if(index(line,'NO LB      ZA    FRAG').ne.0)then
         do i=1,n
            read(ich,'(a)',end=200)line
            call readl(line,xx,nn)
            xyzrd(1:3)=xx(nn-2:nn)
            rmsd=sqrt(sum((xyzrd(1:3)-xyz(1:3,i))**2))
            if(rmsd.gt.1.d-4) then 
               write(*,*) 'coords differ at atom ',i
               stop 'input coords and ORCA coords differ'
            endif
         enddo
      endif
      if(index(line,'The uncorrected charge').ne.0)then
         do i=1,ndim
            read(ich,'(a)',end=200)line
!           if(index(line,'s').ne.0) lao(i)=0
!           if(index(line,'p').ne.0) lao(i)=1
!           if(index(line,'d').ne.0) lao(i)=2
!           if(index(line,'f').ne.0) lao(i)=3
            call readl(line,xx,nn)
            iat = idint(xx(2)) + 1
            aoat(i) = iat
         enddo
      endif
!     if(index(line,'E(XC)              :').ne.0)then
!        call readl(line,xx,nn)
!        exc=xx(1)
!     endif
      if(index(line,'Nuclear Repulsion  :').ne.0)then
         call readl(line,xx,nn)
         enuc=xx(1)
      endif
      if(index(line,'Kinetic Energy     :').ne.0)then
         call readl(line,xx,nn)
         ekin=xx(1)
      endif
      if(index(line,'Total Dipole Moment    :').ne.0)then
         call readl(line,xx,nn)
         dip(1:3)=xx(1:3)
      endif
      if(index(line,'The origin for moment calc').ne.0)then
         call readl(line,xx,nn)
         pnt(1:3)=xx(1:3)
      endif
      if(index(line,'The raw cartesian tensor').ne.0)then
         read(ich,'(a)',end=200)line
         call readl(line,xx,nn)
         alp(1)=xx(1)
         read(ich,'(a)',end=200)line
         call readl(line,xx,nn)
         alp(2)=xx(1)
         alp(3)=xx(2)
         read(ich,'(a)',end=200)line
         call readl(line,xx,nn)
         alp(4)=xx(1)
         alp(5)=xx(2)
         alp(6)=xx(3)
      endif
      goto 100
 200  continue

      close(ich)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine rdstvint0(fname,ndim,ex)
      implicit none
      character*(*) fname
      integer ndim
      logical ex

      character*128 line
      real*8 xx(10)
      integer nn,ich

      ich=141
      inquire(file=fname,exist=ex)
      if(.not.ex) then
         write(*,*) 'orca output file not found'
         return
      endif
         
      open(unit=ich,file=fname)
 100  read(ich,'(a)',end=200)line
      if(index(line,'Basis Dimension        Dim').ne.0) then
         call readl(line,xx,nn)
         ndim = idint(xx(nn))
         return
      endif
      goto 100
 200  continue

      close(ich)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine rdstvintp
     .          (fname,n,ndim,xyz,S,P)
      use bascom, only: aoat
      implicit none
      integer n,ndim
      real*8 S(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)
      real*8 xyz(3,n) 

      real*8 xx(10)
      real*8 xyzrd(3),rmsd
      character*128 line
      character*(*) fname
      integer nn,i,ich,iat

      ich=142
      open(unit=ich,file=fname)
 100  read(ich,'(a)',end=200)line
      if(index(line,'OVERLAP INTEGRALS ').ne.0) then
         call rdmat(ich,S,ndim)
      endif
      if(index(line,'DENSITY').ne.0.and.index(line,'*').eq.0) then
         read(ich,'(a)',end=200)line
         call rdmat(ich,P,ndim)
      endif
      if(index(line,'NO LB      ZA    FRAG').ne.0)then
         do i=1,n
            read(ich,'(a)',end=200)line
            call readl(line,xx,nn)
            xyzrd(1:3)=xx(nn-2:nn)
            rmsd=sqrt(sum((xyzrd(1:3)-xyz(1:3,i))**2))
            if(rmsd.gt.1.d-4) then 
               write(*,*) 'coords differ at atom ',i
!              stop 'input coords and ORCA coords differ'
            endif
         enddo
      endif
      if(index(line,'The uncorrected charge').ne.0)then
         do i=1,ndim
            read(ich,'(a)',end=200)line
!           if(index(line,'s').ne.0) lao(i)=0
!           if(index(line,'p').ne.0) lao(i)=1
!           if(index(line,'d').ne.0) lao(i)=2
!           if(index(line,'f').ne.0) lao(i)=3
            call readl(line,xx,nn)
            iat = idint(xx(2)) + 1
            aoat(i) = iat
         enddo
      endif
      goto 100
 200  continue

      close(ich)

      end
