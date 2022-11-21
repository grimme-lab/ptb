!---------------------------------------------------
! numerical gradient (and SP energy) of gTB-RPE-D4
! obtained from call to egtb script
! symmetry (recongnition) included      
!      
! SG, Feb 22
!---------------------------------------------------      
      program gngtb
      implicit none
      real*8  ,allocatable :: xyz(:,:)
      real*8  ,allocatable :: g  (:,:)
      real*8  ,allocatable :: dum(:,:)
      integer ,allocatable :: at(:)
      integer ,allocatable :: dgen(:)
      integer ,allocatable :: ict(:,:)
      real*8  trans(9,120)
      real*8  step, er, el, e
      integer i, j, n, ii
      integer ntrans  
      character*2 asym
      logical highsym

      step=0.02d0 ! num. diff. step in Bohr

!--------------------------- read coord
      call rd0('coord',n)
      allocate(at(n),xyz(3,n),g(3,n),dum(3,n),ict(n,120),dgen(n))
      call rd(.false.,'coord',n,xyz,at)

      e = 0
      g = 0

!--------------------------- energy at input point
      call runenergy     
      call rdenergy(e)   
!---------------------------      
      call system('mv energy energy.tmp')
      if(abs(e).lt.1d-6) stop 'no energy'

      if(n.eq.2.and.abs(xyz(1,1)).lt.1d-6.and.
     .              abs(xyz(2,1)).lt.1d-6.and.
     .              abs(xyz(1,2)).lt.1d-6.and.
     .              abs(xyz(2,2)).lt.1d-6)
     .then ! case diatomic on z
         xyz(3,1)=xyz(3,1)+step
!---------------------------    right step
         call wrcoord(n,at,xyz)
         call runenergy     
         call rdenergy(er)   
!---------------------------      
         xyz(3,1)=xyz(3,1)-step*2d0
!---------------------------    left one
         call wrcoord(n,at,xyz)
         call runenergy     
         call rdenergy(el)   
!---------------------------    back
         xyz(3,1)=xyz(3,1)+step
         call wrcoord(n,at,xyz)
         g(3,1)= (er-el)/(2d0*step)
         g(3,2)=-(er-el)/(2d0*step)
      else
         call getsymmetry2(.true.,n,at,xyz,0.01d0,ntrans,ict,trans) 
         if(ntrans.gt.1) then ! symmetric
            dgen = 0
            do i=1,n
               do j=1,ntrans
                  if(ict(i,j).eq.i) dgen(i)=dgen(i)+1   
               enddo
               dgen(i)=ntrans / dgen(i)
            enddo
         else
            dgen = 1
         endif
         do i=2,n ! first atom added later
            if(ntrans.gt.1.and.sum(abs(g(:,i))).gt.1d-12) cycle  !  symmetry related atom already done
            do j=1,3
               if(ntrans.gt.1.and.abs(xyz(j,i)).lt.1d-12) cycle  ! =0 by symmetry
               write(*,*) 'atom ',i,' dir',j,
     .                    ' degen. by symmetry ',dgen(i)
               xyz(j,i)=xyz(j,i)+step
!---------------------------    right step
               call wrcoord(n,at,xyz)
               call runenergy     
               call rdenergy(er)   
               xyz(j,i)=xyz(j,i)-step*2d0
!---------------------------    left one
               call wrcoord(n,at,xyz)
               call runenergy     
               call rdenergy(el)   
!---------------------------    back
               xyz(j,i)=xyz(j,i)+step
               call wrcoord(n,at,xyz)
               g(j,i)= (er-el)/(2d0*step)
            enddo
            if(ntrans.gt.1) then
               dum = 0
               dum(1:3,i) = g(1:3,i)
               g(1:3,i) = 0
               call grdsym(dum,n,ntrans,ict,trans) ! symmetrize
               g = g + dum * dble(dgen(i))
            endif
         enddo
         g(1,1)=-sum(g(1,2:n)) ! use other forces to compute first atom
         g(2,1)=-sum(g(2,2:n))
         g(3,1)=-sum(g(3,2:n))
      endif

!---------------------------  
!     TM output style      
      i=1
      open(unit=43,file='gradient')
      write(43,'(''$grad'')')
      write(43,'(''  cycle = '',i6,4x,''SCF energy ='',F18.11,3x,
     .              ''|dE/dxyz| ='',F10.6)')i,e,sum(abs(g))    
      do i=1,n
         write(43,'(3(F20.14,2x),4x,a2)')xyz(1,i),xyz(2,i),xyz(3,i),
     .                                   asym(at(i))
      enddo
      do i=1,n
         write(43,'(3D22.13)')g(1,i),g(2,i),g(3,i)
      enddo
      write(43,'(''$end'')')
      close(43)
      call system('mv energy.tmp energy')

      end

!-----------------------------------------------------------      

      subroutine runenergy

!     call system('egtb')
      call system(
     .'rm -rf ptb_dump gradient; ptb coord -energy -clean > tmp')

      end

!-----------------------------------------------------------      

      subroutine rdenergy(e)
      implicit none
      real*8 e
      character*80 atmp
      integer i

      open(unit=2,file='energy')
      read(2,*) atmp
      read(2,*) i,e 
      close(2)

      end

!-----------------------------------------------------------      

      subroutine wrcoord(n,at,xyz)
      implicit none
      integer at(*),n,i
      real*8 xyz(3,*)
      character*3  asym

      open(unit=2,file='coord')
      write(2,'(''$coord'')')
      do i=1,n
         write(2,'(3F22.14,5x,a2)')xyz(1:3,i),asym(at(i))
      enddo
      write(2,'(''$end'')')
      close(2)

      end

