Subroutine get_schoenflies (n, iat, xyz, sfsym, paramar)
      Use iso_c_binding
      Implicit None

      Interface c_interface
    !Interface to c routine for symmetry recognition
    !attypes are Atom types as integers (e.g 6 for Carbon etc...)
    !x,y,z must be one dimensional sequential(!) arrays of real*8
    !symbol is the recognized schoenflies symbol
         Subroutine schoenflies (natoms, attypes, x, y, z, symbol, &
        & paramar) bind (C, name="schoenflies")
            Use iso_c_binding
            import
            Implicit None
            Integer (c_int), Intent (In) :: natoms
            Type (c_ptr), value :: attypes
            Type (c_ptr), value :: x
            Type (c_ptr), value :: y
            Type (c_ptr), value :: z
            Type (c_ptr), Intent (InOut) :: symbol
            Type (c_ptr), value :: paramar
        !Character (kind=c_char), Intent (out)  :: symbol(*)
         End Subroutine schoenflies
      End Interface c_interface

     !Dummy Arguments
      Character (Len=*)     :: sfsym
      Integer, Intent (In)  :: n
      Integer, Intent (In)  :: iat (n)
      Real * 8, Intent (In) :: xyz (3, n)
      Real * 8, Intent (In) :: paramar (11)


    !local variables for passing to c routine:
      Integer (c_int) :: natoms
      Integer (c_int), Allocatable, Target, Dimension (:) :: attypes
      Real (c_double), Allocatable, Target, Dimension (:) :: x
      Real (c_double), Allocatable, Target, Dimension (:) :: y
      Real (c_double), Allocatable, Target, Dimension (:) :: z
      Real (c_double), Allocatable, Target, Dimension (:) :: c_paramar

     !passing strings to c and get back values requires some nasty pointer
     !stuff, as ISO-C treats strings as pointer-to-char
      Type (c_ptr) :: csptr			 !C-Compatible Pointer
      Character, Pointer, Dimension (:) :: fsptr !Fortran Pointer
      Integer :: sflen (1)			 !Needed for correct call of c_f_pointer

    !local stack:
      Integer :: i

    !Allocate memory for copies of iat, xyz...
    !As they will be passed to C as Pointers, they will be actually the working
    !memory for the C code. They must be allocated with "target" attribute.

      Allocate (attypes(n))
      Allocate (x(n))
      Allocate (y(n))
      Allocate (z(n))
      Allocate (c_paramar(11))

    !now, copy contents
      natoms = n
      attypes (:) = iat (:)
      x (:) = xyz (1, :)
      y (:) = xyz (2, :)
      z (:) = xyz (3, :)
      c_paramar = paramar

      sfsym = sfsym // C_NULL_CHAR
      sflen (1) = 3 !3 Characters soule be enough for schoenflies symbol
      csptr = c_loc (sfsym)! Get address of dummy string
      
      Call schoenflies (natoms, c_loc(attypes), c_loc(x), c_loc(y), &
     & c_loc(z), csptr, c_loc(c_paramar))
     
     Call c_f_pointer (csptr, fsptr, sflen)
      Write (sfsym,*) fsptr !Write symbol back to dummy argument
      Nullify (fsptr)

    !deallocate arrays:
      Deallocate (attypes, x, y, z,c_paramar)
End Subroutine get_schoenflies

!Program test_call
!      Use iso_c_binding, Only: C_CHAR, C_NULL_CHAR
!      Implicit None
!      !local stack:
!      Character (Len=10) :: sfsym
!      Integer :: natoms
!      Real * 8 :: paramar (11)  !parameter array for get_schoenflies_
!      
!      !local heap:
!      Real * 8, Allocatable :: xyz (:, :)
!      Integer, Allocatable :: iat (:)
!      
!
!
!      !parameters for symmetry recognition:
!      paramar (1) = - 1         ! verbose, increase for more detailed output (to stdout)
!      paramar (2) = 20          ! MaxAxisOrder
!      paramar (3) = 200         ! MaxOptCycles
!      paramar (4) = 0.001d0     ! ToleranceSame
!      paramar (5) = 0.5d0       ! TolerancePrimary
!      paramar (6) = 0.1d0       ! ToleranceFinal, THIS IS THE IMPORTANT VALUE
!      paramar (7) = 0.5d0       ! MaxOptStep
!      paramar (8) = 1.0D-7      ! MinOptStep
!      paramar (9) = 1.0D-7      ! GradientStep
!      paramar (10) = 1.0D-10    ! OptChangeThreshold
!      paramar (11) = 5          ! OptChangeHits
!
!     !Read coord file
!      Call rd0 ("coord", natoms)
!      Allocate (xyz(3, natoms))
!      Allocate (iat(natoms))
!      Call rd (.True., "coord", natoms, xyz, iat)
!      sfsym = ""
!      
!     !determine Point_group
!      Call get_schoenflies (natoms, iat, xyz, sfsym, paramar)
!      Write (*,*) "get_schoenflies returned Point Group: ", sfsym
!
!      Deallocate (xyz, iat)
!End Program test_call
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine getsymmetry (pr,n, iat, xyz, symthr, maxatdesy, highsym)
      Use iso_c_binding, Only: C_CHAR, C_NULL_CHAR
      Implicit None
      Integer n, iat(n), maxatdesy, extcode
      real*8  xyz(3,n)
      real*8 symthr
      Character*3 sfsym
      logical pr, highsym
      Character*4 atmp    

      Real * 8 :: paramar (11)  !parameter array for get_schoenflies_
      Real * 8 :: xx(2)                                                 
      integer  :: nn                                                   

      highsym=.false.
      if(n.gt.maxatdesy) return

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

      if(pr)write(*,'(a3,'' symmetry found (for desy threshold: '',e9.2,'')'')')sfsym, symthr

      call readl(sfsym,xx,nn)
 
      if(idint(xx(1)).gt.2.or.n.le.2) highsym=.true. ! atom or diatomic
      if(index(sfsym,'in').ne.0     ) highsym=.true. ! linear
      if(index(sfsym,'oh').ne.0     ) highsym=.true. ! Oh    
      if(index(sfsym,'ih').ne.0     ) highsym=.true. ! Ih    

End

