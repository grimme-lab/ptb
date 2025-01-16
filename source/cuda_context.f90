module cuda_
   use gtb_accuracy, only: i4
   use accel_lib
   implicit none
   public :: initialize_ctx, ctx
   private
   type(cuda_context), allocatable :: ctx
contains
   subroutine initialize_ctx()

      integer(i4) :: err

      allocate(ctx)
      err = load_cuda()
      if (err /= 0) then
         stop "CUDA acceleration library could not be loaded."
      end if
      
      ctx = cuda_init(err)
      if (err /= 0) then
         stop "Could not initialize GPU acceleration library."
      end if

   end subroutine initialize_ctx
end module cuda_