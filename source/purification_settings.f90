module purification_settings
   use gtb_accuracy, only: wp
   implicit none

   !> Different calculation types
   integer, parameter :: mcweeny = 1 ! CLassical McWeeny purification

   !> S metric for purification
   integer, parameter :: inv = 1
   integer, parameter :: inv_sqrt = 2 
   integer, parameter :: sqrt_ = 3

   type :: tPurificationSet
   
      !> S metric
      integer :: metric 

      !> Purification type
      integer :: method 
      
      !> Printout level during calculation
      integer :: plvl

      !> Development mode (calculate solve2 as well and divergence between results)
      logical :: dev = .true.

   contains
      procedure :: settings => initialize_purification
   end type tPurificationSet

contains
   
   !> Setup purification setting
   subroutine initialize_purification(self, unit)

      !> Structure
      class(tPurificationSet), intent(inout) :: self

      !> I/O unit number
      integer, intent(in) :: unit

   end subroutine initialize_purification

end module purification_settings