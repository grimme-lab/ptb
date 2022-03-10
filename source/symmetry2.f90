Subroutine getsymmetry2(pr,n, iat, xyz, symthr, ntrans, ict, trans)
      Use iso_c_binding, Only: C_CHAR, C_NULL_CHAR
      Implicit None
      Integer n, iat(n)
      real*8  xyz(3,n)
      real*8 symthr
      Integer ict(n,120), ntrans
      real*8  trans(9,120)

      Character*4 sfsym
      logical pr, highsym
      Character*4 atmp    

      Real * 8 :: paramar (11)  !parameter array for get_schoenflies_
      Real * 8 :: xx(2)                                                 
      integer  :: nn                                                   

      if(pr)write(*,*)
!parameters for symmetry recognition:
      paramar (1) = - 1         ! verbose, increase for more detailed output (to stdout)
      paramar (2) = 10          ! MaxAxisOrder
      paramar (3) = 100         ! MaxOptCycles
      paramar (4) = 0.001d0     ! ToleranceSame
      paramar (5) = 0.5d0       ! TolerancePrimary
      paramar (6) = symthr      ! ToleranceFinal, THIS IS THE IMPORTANT VALUE
      paramar (7) = 0.5d0       ! MaxOptStep
      paramar (8) = 1.0D-7      ! MinOptStep
      paramar (9) = 1.0D-7      ! GradientStep
      paramar (10) = 1.0D-8     ! OptChangeThreshold
      paramar (11) = 5          ! OptChangeHits

      atmp='    '
      Call get_schoenflies (n, iat, xyz, atmp, paramar)

!TM stuff (trafo table)
      sfsym(1:3)=atmp(2:4)
      if(sfsym(1:1).eq.'D') sfsym(1:1)='d'
      if(sfsym(1:1).eq.'C') sfsym(1:1)='c'
      if(sfsym(1:1).eq.'T') sfsym(1:1)='t'
      if(sfsym(1:1).eq.'O') sfsym(1:1)='o'
      if(sfsym(1:1).eq.'I') sfsym(1:1)='i'
      if(sfsym.eq.'dih') sfsym='d6h'
      if(sfsym.eq.'civ') sfsym='c6v'
      if(sfsym(3:3).gt.'v'.or.sfsym(3:3).lt.'a') sfsym(3:3)=' '
      sfsym(4:4)=' '

      call symtrans(sfsym,1.0d0*symthr,xyz,n,ntrans,ict,trans)

      if(pr.and.ntrans.gt.1) then
         write(*,'(a3,'' symmetry found (for desy threshold: '',e9.2,'')'')')sfsym, symthr
         write(*,'(''# symmetry operations : '',i4)') ntrans
      endif

End
