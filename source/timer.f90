module timer
   use gtb_accuracy, only: wp

   implicit none
   private 
   public :: tTimer

   !> timer class 
   type :: tTimer
      private

      !> Number of timers
      integer, public :: n
      
      !> Wall-time for each timer
      real(wp), allocatable :: wall_time(:)
      
      !> CPU time for each timer
      real(wp), allocatable :: cpu_time(:)

      !> Distinguish running timers
      logical, private,allocatable :: running(:)
      
      !> Timer Name
      character(len=40), allocatable :: tag(:)

      !> Total elapsed wall-time
      real(wp) :: total_wall_time

      !> Total elapsed cpu time
      real(wp) :: total_cpu_time
   contains
      procedure :: new => allocate_timer
      procedure :: deallocate => deallocate_timer

   end type tTimer
contains

   !> initialize timer instance
   subroutine allocate_timer(self,n,verbose)
      
      implicit none
      
      !> instance of timer
      class(tTimer),intent(inout) :: self

      !> number of timers
      integer, intent(in)           :: n
      
      !> if verbose
      logical, intent(in), optional :: verbose

      real(wp) :: time_cpu
      real(wp) :: time_wall

      ! capture negative values !
      if (n < 1) return
      
      self%n = n
      
      allocate( self%wall_time(0:n), source = 0.0_wp )
      allocate( self%cpu_time(0:n),  source = 0.0_wp )
      allocate( self%running(n), source =.false. )
      allocate( self%tag(n) ); self%tag = ' '

   end subroutine allocate_timer

   !> To deallocate memory
   subroutine deallocate_timer(self)
      
      implicit none
      
      !> instance of timer
      class(tTimer),intent(inout) :: self
      
      self%n = 0
      self%total_wall_time = 0
      self%total_cpu_time  = 0
      if (allocated(self%wall_time))   deallocate(self%wall_time)
      if (allocated(self%cpu_time))    deallocate(self%cpu_time)
      if (allocated(self%running)) deallocate(self%running)
      if (allocated(self%tag)) deallocate(self%tag)
   end subroutine deallocate_timer

end module timer