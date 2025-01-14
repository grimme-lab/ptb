module norms
   use gtb_accuracy, only: wp
contains
   !> Gershgorin norm 
   real(wp) pure function gershgorin(matrix) result(gersh)

      !> Input matrix
      real(wp), intent(in) :: matrix(:,:)
      
      !> Locals
      integer :: i

      gersh=0.0_wp
      do i = 1, size(matrix,2)
         gersh=max(gersh,sum(abs(matrix(:,i))))
      enddo

   end function gershgorin

   !> Frobenius norm
   real(wp) pure function frobenius(matrix) result(frob)

      !> Input matrix
      real(wp), intent(in) :: matrix(:,:)

      frob=sqrt(sum(matrix**2))
      
   end function frobenius


   !> Infinity norm
   real(wp) pure function infinity(matrix) result(out)

      !> Input matrix
      real(wp), intent(in) :: matrix(:,:)

      !> Locals 
      integer :: i
      real(wp), allocatable :: inf(:)

      allocate(inf(size(matrix,2)))
      
      do i = 1, size(matrix, 2)
         inf(i) = sum(abs(matrix(i,:))) 
      enddo

      out = maxval(inf)

   end function infinity

end module norms