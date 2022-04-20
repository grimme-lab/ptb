!---------------------------------------------------
! SG, Feb 22
!---------------------------------------------------      
      program eh2  
      implicit none
      real*8  ,allocatable :: xyz(:,:)
      integer ,allocatable :: at(:)
      integer i, j, n, iat, hleft, hright
      integer ielem(86)   
      integer nh   (86)   
      real*8  eref (86)   
      real*8  etb  (86)   
      real*8  eref_mol, rkref, rk, e_mol, hfak
      character*2 asym
      character*128 atmp
      logical ex

      nh    = 0
      nh(1) = 2
      nh(3) = 1
      nh(4) = 2
      nh(5) = 3
      nh(6) = 4
      nh(7) = 3
      nh(8) = 2
      nh(9) = 1
      nh(12)= 2
      nh(13)= 3
      nh(14)= 4
      nh(15)= 3
      nh(16)= 2
      nh(17)= 1
      nh(35)= 1
      nh(53)= 1
      nh(85)= 1
      nh(34)= 2
      nh(52)= 2
      nh(84)= 2
      nh(33)= 3
      nh(51)= 3
      nh(53)= 3
      nh(32)= 4
      nh(50)= 4
      nh(82)= 4
      nh(31)= 3
      nh(49)= 3
      nh(81)= 3
      nh(19)= 1
      nh(37)= 1
      nh(55)= 1
      nh(20)= 2
      nh(38)= 2
      nh(56)= 2

      nh(30)= 2
      nh(48)= 2
      nh(80)= 2

      nh(29)= 1
      nh(47)= 1
      nh(79)= 1

      nh(28)= 2
      nh(46)= 2
      nh(78)= 2

      nh(27)= 5
      nh(45)= 5
      nh(77)= 5

      nh(26)= 4
      nh(44)= 4
      nh(76)= 4

      nh(25)= 3
      nh(43)= 3
      nh(75)= 3

      nh(24)= 4
      nh(42)= 4
      nh(74)= 4

      nh(23)= 5
      nh(41)= 5
      nh(73)= 5

      nh(22)= 4
      nh(40)= 4
      nh(72)= 4

      nh(21)= 3
      nh(39)= 3
      nh(57)= 3


      open(unit=1,file='/tmp1/grimme/HYDRIDE/.eref')
 10   read(1,*,end=11) i,eref(i)
      goto 10
 11   close(1)
      open(unit=1,file='/tmp1/grimme/HYDRIDE/.etb')
 20   read(1,*,end=21) i,etb(i)
      goto 20
 21   close(1)

!--------------------------- read coord
      call rd0('coord',n)
      allocate(at(n),xyz(3,n))
      call rd(.false.,'coord',n,xyz,at)

!--------------------------- formula
      hright= 0
      ielem = 0
      do i=1,n 
         ielem(at(i)) = ielem(at(i)) + 1
         hright = hright + nh(at(i))
      enddo
      hleft = hright - ielem(1)

!--------------------------- energy at input point

      inquire(file='.eh2',exist=ex)
      if(ex) then
             open(unit=1,file='.eh2')
             read(1,*) eref_mol
             close(1)
      else
 atmp='rmec; cefine -bas def2-QZVP -gf -func wb97x-v -noopt -d4 -senex; ridft>tmp'
             call system(atmp)
             call system("eiger --holumogap | grep 'eV'")
             call rdenergy(eref_mol)
             open(unit=1,file='.eh2')
             write(1,*) eref_mol
             close(1)
      endif

      call system('rmec; egtb')
      call rdenergy(e_mol)   

!-----------------------------------------------------------      

!     write(*,*) '         element      E(ref)        E(TB)'
      rkref = 0
      rk    = 0
      do i=1,n
         iat = at(i)
         rkref = rkref + eref(iat)
         rk    = rk    + etb (iat)
!        write(*,*) iat,eref(iat),etb(iat)
      enddo

      hfak = 0.5d0* float(hleft)

!     write(*,*)  'n H2 ',hfak   
!     write(*,*) '      E H2 (ref)       E mol (ref)'
!     write(*,*)  eref(1), eref_mol
!     write(*,*) '      E H2 (TB)        E mol (TB) '
!     write(*,*)  etb (1), e_mol

      rkref = rkref - eref_mol - hfak * eref(1)
      rk    = rk    - e_mol    - hfak * etb (1)

      write(*,'(''reference '',f10.6,f9.2,3x,''TB '',f10.6,f9.2,3x,&
   &  ''error '',f7.2,3x,'' % '',f7.2,3x,'' per atom '',f7.2)') &
   &  rkref,rkref*627.51,rk,rk*627.51,(rkref-rk)*627.51,100d0*((rkref-rk)/rkref),&
   &  627.41*(rkref-rk)/dble(n)

!-----------------------------------------------------------      
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

