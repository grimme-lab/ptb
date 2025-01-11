module matops
   use gtb_accuracy, only: wp

contains

   !> Calculate trace of a matrix
   real(wp) pure function trace(matrix) result(tr)

      !> Input matrix
      real(wp), intent(in) :: matrix(:,:)
      
      !> Counter
      integer :: i

      tr = 0.0_wp
      do i = 1, size(matrix,2)
         tr = tr + matrix(i,i)
      enddo
   end function


end module matops