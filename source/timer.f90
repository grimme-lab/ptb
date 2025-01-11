module timer
   use gtb_accuracy, only: wp, i8

   implicit none
   private 
   public :: tTimer

   intrinsic :: system_clock, cpu_time
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
      logical, allocatable, public :: running(:)
      
      !> Timer Name
      character(len=40), allocatable :: tag(:)

      !> Total elapsed wall-time
      real(wp) :: total_wall_time

      !> Total elapsed cpu time
      real(wp) :: total_cpu_time
   contains
      procedure :: new => allocate_timer
      procedure :: deallocate => deallocate_timer
      procedure :: finalize => write_results
      procedure :: click => timer_click
      
      procedure, private :: start_timing
      procedure, private :: stop_timing

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

      call self%deallocate
      ! capture negative values !
      if (n < 1) return
      
      self%n = n
      
      allocate( self%wall_time(0:n), source = 0.0_wp )
      allocate( self%cpu_time(0:n),  source = 0.0_wp )
      allocate( self%running(0:n), source =.false. )
      allocate( self%tag(0:n) ); self%tag = ' '

      ! start total timer !
      call self%start_timing(0)
   
   end subroutine allocate_timer


   !> Start/stop button for each individual timer member
   subroutine timer_click(self, i)
      
      implicit none
      
      !> instance of timer
      class(tTimer),intent(inout) :: self

      !> index
      integer,intent(in) :: i
      
      ! check if appropriate index is given !
      if (i > self%n .or. i < 1) return
      
      ! switcher between start/stop status !
      if (self%running(i)) then
         call self%stop_timing(i)
      else
         call self%start_timing(i)
      endif
      
   end subroutine timer_click

   !> To retrieve the current CPU and wall time
   subroutine timing(time_cpu,time_wall)
      
      implicit none
   
      real(wp),intent(out) :: time_cpu
      real(wp),intent(out) :: time_wall
      
      !> current value of system clock (time passed from arbitary point)
      integer(i8) :: time_count
      
      !> number of clock ticks per second (conversion factor b/n ticks and seconds)
      integer(i8) :: time_rate
      integer(i8) :: time_max
      
      call system_clock(time_count,time_rate,time_max)
      call cpu_time(time_cpu)
   
      ! elapsed time in seconds !
      time_wall = real(time_count,wp)/real(time_rate,wp)
   
   end subroutine timing


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

   !> Stop all the timers
   subroutine write_results(self)

      intrinsic :: pack
      class(tTimer), intent(inout) :: self

      ! Find the indices of the .false. elements
      integer, allocatable :: indices(:)
      integer :: i
      ! Stop timers if they are still running
      if (any(self%running)) then
         indices = pack([(i, i=0, size(self%running)-1)], self%running)
      endif
      do i= 1, size(indices)
         call self%stop_timing(indices(i))
      enddo

   end subroutine write_results

   !> Start timer
   subroutine start_timing(self,i)
      
      implicit none

      !> instance of timer
      class(tTimer),intent(inout) :: self
      
      !> index 
      integer,intent(in) :: i
      
      real(wp) :: time_cpu
      real(wp) :: time_wall

      ! sanity check !
      if (self%running(i))&
         & error stop 'timer already running, cannot start it twice'

      if (i >= 0 .and. i <= self%n) then 
         self%running(i) = .true. 
      else
         error stop 'False timer index'
      endif

      call timing(time_cpu,time_wall)
      self%cpu_time(i) =  self%cpu_time(i) - time_cpu
      self%wall_time(i) = self%wall_time(i) - time_wall

   end subroutine start_timing


   !> To stop counting
   subroutine stop_timing(self,i)
      
      implicit none
   
      !> instance of timer
      class(tTimer),intent(inout) :: self
   
      !> index
      integer,intent(in) :: i
      
      real(wp) :: time_cpu
      real(wp) :: time_wall

      ! sanity check !
      if (.not. self%running(i))&
         & error stop 'timer already running, cannot start it twice'

      if (i >= 0 .and. i <= self%n) then 
         self%running(i) = .false. 
      else
         error stop 'False timer index'
      endif

      
      call timing(time_cpu,time_wall)
      self%cpu_time(i) =  self%cpu_time(i) + time_cpu
      self%wall_time(i) = self%wall_time(i) + time_wall

   end subroutine stop_timing


end module timer