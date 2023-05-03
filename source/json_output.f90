module json_output

   use iso_fortran_env, only : wp => real64
   implicit none

   private

   public :: write_json

contains

   subroutine write_json(filename, nat, atomtype, z, q, wbo, dipmom)
      ! Writes the charges to the JSON output file.
      ! nat: number of atoms
      ! atomtype: atom types
      ! q: charges
      ! i: loop index
      character(len=*),   intent(in)  :: filename
      integer,            intent(in)  :: nat
      integer,            intent(in)  :: atomtype(nat)
      real(wp),           intent(in)  :: q(nat)
      real(wp),           intent(in)  :: z(nat)
      real(wp),           intent(in)  :: wbo(nat,nat)
      real(wp),           intent(in)  :: dipmom(3)
      integer                         :: i,j,myunit,ierr
      logical                         :: first = .false.

      open(newunit=myunit, file=filename, status='replace', action='write', &
         access='sequential', form='formatted', iostat=ierr)
      if (ierr /= 0) then
         write(*, '(a, i0, a)', advance='no') 'Error opening file ', myunit, &
            ' for writing.'
         stop
      end if

      ! Write the JSON object
      ! Opening curly bracket of the JSON object
      write(myunit, '(a)') '{'

      ! First object: "number of atoms"
      write(myunit, '(3x,a,i0,a)') '"Number of atoms": ', nat, ','


      ! Second object: "atomic charge"
      write(myunit, '(3x,a)') '"Atomic charge"  : ['
      do i = 1, nat-1
         write(myunit, '(f8.4,a,/)', advance='no') z(i)-q(i),","
      end do
      write(myunit, '(f8.4)') z(nat)-q(nat)
      write(myunit, '(a)', advance='no') '  ]'

      ! Separating comma with advanced formatting:
      write(myunit, '(a)') ','

      ! Third object: "atomic number"
      write(myunit, '(3x,a)') '"Atomic number"  : ['
      do i = 1, nat-1
         write(myunit, '(2x,i0,a,/)', advance='no') atomtype(i),","
      end do
      write(myunit, '(2x,i0)') atomtype(nat)
      write(myunit, '(a)', advance='no') '  ]'

      ! Separating comma with advanced formatting:
      write(myunit, '(a)') ','

      ! Fourth object: Wiberg bond order
      write(myunit, '(3x,a)') '"Wiberg bond order"  : ['
      do i=1,nat-1
         do j=i,nat
            if ( wbo(j,i) .gt. 0.01 ) then
               if (first) then
                  write(myunit, '(a)') ','
               end if
               first=.true.
               write(myunit, '(2x,a,i0,a,i0,a,f8.4,a)', advance='no') "[ ",i,", ",j,",",wbo(j,i)," ]"
            end if
         enddo
      enddo
      write(myunit, '(/,a)', advance='no') '  ]'

      ! Separating comma with advanced formatting:
      write(myunit, '(a)') ','

      ! Fifth object: dipole moment
      write(myunit, '(3x,a)') '"Dipole moment (x,y,z)"  : ['
      write(myunit, '(f8.4,a,/,f8.4,a,/,f8.4)', advance='no') dipmom(1),", ",dipmom(2),", ",dipmom(3)
      write(myunit, '(/,a)', advance='no') '  ]'

      ! Enclosing curly bracket of the JSON object
      write(myunit, '(/,a)') '}'

      close(myunit)

   end subroutine write_json

end module json_output
