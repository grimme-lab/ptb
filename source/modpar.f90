program modpar 
      use iso_fortran_env, only : wp => real64
   
      implicit none
      integer i,j,n,io
      integer itabrow6 
      real*8 tmp    
   real(wp),parameter :: rcov(118) = 1.889725949_wp * [ & 
   & 0.32_wp,0.46_wp, & ! H,He
!  & 0.24_wp,0.46_wp, & ! H,He H changed
   & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
   & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
   & 1.76_wp,1.54_wp, & ! K,Ca
   &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
   &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
   &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
   & 1.89_wp,1.67_wp, & ! Rb,Sr
   &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
   &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
   &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
   & 2.09_wp,1.76_wp, & ! Cs,Ba
   &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
   &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
   &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
   &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
   &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
   & 2.01_wp,1.81_wp, & ! Fr,Ra
   &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og


      real*8     glob_par  (20)     
      real*8     fock_par  (10,86)
      real*8     fock_lev  (10,86)
      real*8     fock_lev2 (10,86)
      real*8     ener_par1 (10,86)
      real*8     ener_par2 (10,86)
      real*8     expscal (4,10,86)
      real*8   shell_xi  (10,86)
      real*8   shell_cnf1(10,86)
      real*8   shell_cnf2(10,86)
      real*8   shell_cnf3(10,86)
      real*8   shell_cnf4(10,86)
      real*8   shell_resp(10,86,2)

      io=6
      n = 86            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read parameter file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      glob_par =0
      fock_par =0
      fock_lev =0
      fock_lev2=0
      expscal  =0
      shell_resp=0

      open(unit=1,file='~/.atompara')
      read(1,*) glob_par (1:10)
      read(1,*) glob_par(11:20)
      do i=1,n 
         read(1,*) j
         read(1,*) fock_par  (1:10,j)    ! 1-10   fock
         read(1,*) fock_lev  (1:10,j)    ! 11-20  fock
         read(1,*) fock_lev2 (1:10,j)    ! 31-40  fock 
         read(1,*) expscal  (1,1:10,j)   ! 21-30  fock
         read(1,*) expscal  (2,1:10,j)   ! 41-50  E
         read(1,*) ener_par1 (1:10,j)    ! 51-60  E        
         read(1,*) ener_par2 (1:10,j)    ! 61-70  E   
         read(1,*) shell_xi  (1:10,j)    ! 71-80  shellQ
         read(1,*) shell_cnf1(1:10,j)    ! 81-90    "
         read(1,*) shell_cnf2(1:10,j)    ! 91-100   "
         read(1,*) shell_cnf3(1:10,j)    ! 101-110  "
         read(1,*) expscal  (3,1:10,j)   ! 111-120  "
         read(1,*) shell_cnf4(1:10,j)    ! 121-130  "
         read(1,*) shell_resp(1:10,j,1)  ! 121-130  "
         read(1,*) shell_resp(1:10,j,2)  ! 121-130  "
         write(*,*) 'read',j
      enddo
      close(1)

!     do j=1,86
!        if(fock_par(1,j).lt.1d-6) cycle
!        write(42,*) j,shell_cnf4(7,j)
!     enddo
!     stop
!     do i=1,86
!        if(i.lt.21.or.i.gt.31)then
!        shell_respa(1:5,i)=shell_cnf1(6:10,i)
!        shell_respa(10 ,i)=shell_xi  (   9,i)
!        shell_xi  (   9,i)=0
!        shell_cnf1(6:10,i)=0
!        endif
!     enddo
!     do i=1,86
!        shell_resp(10,i,1:2)=shell_resp(10,i,1:2)*2d0
!     enddo
!     do i=1,86
!        shell_cnf1(9,i)=shell_cnf4(7,i)  ! 127 -> 89
!        shell_cnf2(9,i)=shell_cnf4(8,i)  ! 128 -> 99
!        shell_cnf4(7:8,i)=0
!     enddo
      do i=1,86
         expscal(1,8,i) = shell_cnf4(3,i)
      enddo

      open(unit=io,file='~/atompara')
      write(io,111) glob_par (1:10)
      write(io,111) glob_par(11:20)
      do j=1,86
!        if(abs(fock_par(1,j)).lt.1d-8) cycle
         write(*,*) 'write',j
         write(io,*  ) j
         write(io,111) fock_par  (1:10,j)    ! 1-10   fock
         write(io,111) fock_lev  (1:10,j)    ! 11-20  fock
         write(io,111) fock_lev2 (1:10,j)    ! 31-40  fock 
         write(io,111) expscal  (1,1:10,j)   ! 21-30  fock
         write(io,111) expscal  (2,1:10,j)   ! 41-50  E
         write(io,111) ener_par1 (1:10,j)    ! 51-60  E        
         write(io,111) ener_par2 (1:10,j)    ! 61-70  E   
         write(io,111) shell_xi  (1:10,j)    ! 71-80  shellQ
         write(io,111) shell_cnf1(1:10,j)    ! 81-90    "
         write(io,111) shell_cnf2(1:10,j)    ! 91-100   "
         write(io,111) shell_cnf3(1:10,j)    ! 101-110  "
         write(io,111) expscal  (3,1:10,j)   ! 111-120  "
         write(io,111) shell_cnf4(1:10,j)    ! 121-130  "
         write(io,111) shell_resp(1:10,j,1)  ! 121-130  "
         write(io,111) shell_resp(1:10,j,2)  ! 121-130  "
      enddo
      close(1)
111   format(20F15.10)
      end

INTEGER FUNCTION iTabRow6(i)
      implicit none
      INTEGER i
 
      iTabRow6=0
      If (i.gt. 0 .and. i.le. 2) Then
         iTabRow6=1
      Else If (i.gt. 2 .and. i.le.10) Then
         iTabRow6=2
      Else If (i.gt.10 .and. i.le.18) Then
         iTabRow6=3
      Else If (i.gt.18 .and. i.le.36) Then
         iTabRow6=4
      Else If (i.gt.36 .and. i.le.54) Then
         iTabRow6=5
      Else If (i.gt.54) Then
         iTabRow6=6
      End If
 
End

