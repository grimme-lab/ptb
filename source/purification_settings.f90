module purification_settings
   use gtb_accuracy, only: wp
   implicit none

   !> Different calculation types
   integer, parameter :: mcweeny = 1 ! CLassical McWeeny purification

   !> S metric for purification
   integer, parameter :: inv = 1
   integer, parameter :: inv_sqrt = 2 
   integer, parameter :: sqrt_ = 3

   type :: tChempotSet

      !> Number of cycles during the chemical potetial search
      integer :: cycles = 40

      !> Chempot guess
      real(wp) :: guess = -1.5_wp

      !> Increment for search
      real(wp) :: increment = 0.1_wp

   end type

   type :: tMetricSet
 
      !> S metric
      integer :: type = inv_sqrt

      !> Number of cycles during iterative calculation of X
      integer :: cycles = 50

      !> Numerical calculation of the metric
      logical :: iterative = .false.

   end type

   type :: tPurificationSet

      !> Purification type
      integer :: type = mcweeny
      
      !> Purification cycles in itertive methods
      integer :: cycles = 40
      
      !> Printout level during calculation (0 = less, 1 = normal, 2 = verbose)
      integer :: prlvl = 2

      !> Development mode (calculate solve2 as well and divergence between results)
      logical :: dev = .true.

      !> (internal) Number of electrons
      integer :: nel

      type(tMetricSet) :: metric
      type(tChempotSet) :: chempot
   contains
      procedure :: settings => initialize_purification
      procedure :: print => print_settings
   end type tPurificationSet

contains
   
   !> Setup purification setting
   subroutine initialize_purification(self, unit)

      character(len=*), parameter :: eq = "="
      character(len=*), parameter :: colon = ":"
      
      !> Structure
      class(tPurificationSet), intent(inout) :: self

      !> I/O unit number
      integer, intent(in) :: unit
      
      !> Local variables 
      integer :: iostat, ns, nf, nl
      character(len=:), allocatable :: arg1, arg2, arg3
      character(len=80) :: line
      character(len=80) :: str(10) ! Tokens 
      real(wp) :: floats(10) ! Floating numbers
      logical :: logicals(10) ! Booleans

      read_settings: do
         read(unit, '(a)', iostat=iostat) line
         if(iostat .ne. 0) exit read_settings
         call readline(line, floats, str, logicals, ns, nf, nl)

         arg1 = trim(adjustl(str(1)))
         if (arg1(1:1) .eq. '!') cycle ! comment line
         
         arg2 = trim(adjustl(str(2)))
         if (arg2 == eq) then
           
            arg3 = trim(adjustl(str(3)))
            
            select case(arg1)
   
            ! get purification type !
            case('type')
               select case(arg3)
               case('mcweeny')
                  self%type = mcweeny
               case default
                  error stop 'Error: .PUR contains invalid type definition'
               end select
            
            case('metric')
               select case(arg3)
               case('inverse')
                  self%metric%type = inv
               case('inverse_sqrt')
                  self%metric%type = inv_sqrt
               case('sqrt')
                  self%metric%type = sqrt_
               case default
                  error stop 'Error: .PUR contains invalid metric definition.'
               end select
            
            case('cycles')
               if (nf > 0) &
                  self%cycles = int(floats(1))

            case('chempot_cycles')
               if (nf > 0) &
                  self%chempot%cycles = int(floats(1))

            case('chempot_increment')
               if (nf > 0) &
                  self%chempot%increment = int(floats(1))

            end select
         else

            ! logicals !
            select case(arg1)
            case('dev')
               self%dev = .true.
            case('iterative')
               self%metric%iterative= .true.
            endselect         
         
         endif
      enddo read_settings
   
   end subroutine initialize_purification

   subroutine print_settings(self, out)

      !> Purification settings holder
      class(tPurificationSet), intent(in) :: self

      !> I/O unit
      integer,  intent(in) :: out

      ! HEADER !
      write(out,'(/,a)') repeat('*',72)
      write(out,'(a,1x,a,1x,a)') repeat('*',26), "PURIFICATION MODE", repeat('*',27)
      write(out,'(a)') repeat('*',72)

      write(out,'(2x,a)') "__SETTINGS__" 
      write(out,'(2x,a)') "__general__"
      write(out,'(2x,a,6x)',advance='no') "Purification type:            "
      selectcase(self%type)
      case(mcweeny)
         write(out,'(a)') 'McWeeny Grand Canonical Purification'
         if (self%prlvl > 0) &
         write(out,'(2x,a, 8x, i0)') 'Iteration Cycles:           ', self%cycles
      endselect
      
      if (self%prlvl > 0) then
         write(out,'(/,2x,a)') "__metric__"
         write(out,'(2x,a,5x)',advance='no') "S power:                      "
         selectcase(self%metric%type)
         case(inv)
            write(out,'(a)') '-1'
         case(sqrt_)
            write(out,'(a)') '0.5'
         case(inv_sqrt)
            write(out,'(a)') '-0.5'
         endselect
         write(out,'(2x,a,6x,L1)')           "Iterative metric evaluation:  ", self%metric%iterative
         if (self%metric%iterative .and. self%prlvl > 1) &
         write(out,'(2x,a,5x,i0)')           "Number of cycles:             ", self%metric%cycles
      
         write(out,'(/,2x,a)') "__chempot__"
         write(out,'(2x,a,3x,f13.8)')        "Initial Chemical Potential:   ", self%chempot%guess
         if (self%prlvl > 1) &
         write(out,'(2x,a,5x,i0)')           "Number of cycles:             ", self%chempot%cycles

      endif
         

      if(self%prlvl > 1) then
         write(out,'(/,2x,a)') "__chempot__"
         write(out,'(2x,a,5x,L1)') "Development Mode:              ", self%dev
      endif
      write(out,'(a)') repeat('*',72)
      write(out,'()')

   end subroutine print_settings
end module purification_settings