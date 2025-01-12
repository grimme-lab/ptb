module matops
   use gtb_accuracy, only: wp
   use iso_fortran_env, only: default => output_unit
   implicit none
   intrinsic :: merge
   interface print_matrix
      module procedure :: print_packed_matrix
      module procedure :: print_full_matrix
   endinterface print_matrix

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

   end function trace

   !> Print packed matrix (ndim*(ndim+1)/2)
   subroutine print_packed_matrix(ndim, matrix, name, unit)
      
      !> Number of BFs
      integer, intent(in) :: ndim

      !> Input matrix 
      real(wp), intent(in) :: matrix(ndim*(ndim+1)/2)
      
      !> Matrix name
      character(len=*), intent(in) :: name

      !> Optionally, I/O unit
      integer, intent(in), optional :: unit

      !> Local variables
      integer :: i, iunit
      real(wp) :: tmp(ndim,ndim)
      logical :: verbose

      ! get unit !
      if (present(unit)) then
         iunit = unit
      else
         iunit = default
      end if
      
      call blowsym(ndim,matrix,tmp)
      
      ! Sanity check if matrix is printable !
      verbose = merge(.true., .false., ndim < 15) 
      
      write(iunit,'(a)') repeat("-", 195)
      write(iunit,'(  1x, a, ":" )') name
      if (verbose) then
         do i = 1, ndim
            write(iunit,'(*(f12.4,1x))') tmp(:,i)
         enddo
      else
         write(iunit,'(*(f12.4,1x))') tmp(:15,1)
      endif

      write(iunit,'(a)') repeat("-", 195)


   end subroutine print_packed_matrix


   !> Print full matrix (ndimndim)
   subroutine print_full_matrix(ndim, matrix, name, unit)
      
      !> Number of BFs
      integer, intent(in) :: ndim

      !> Input matrix 
      real(wp), intent(in) :: matrix(ndim,ndim)
      
      !> Matrix name
      character(len=*), intent(in) :: name

      !> Optionally, I/O unit
      integer, intent(in), optional :: unit

      !> Local variables
      integer :: i,  iunit
      logical :: verbose

      ! get unit !
      if (present(unit)) then
         iunit = unit
      else
         iunit = default
      end if
      
      ! Sanity check if matrix is printable !
      verbose = merge(.true., .false., ndim < 15) 
      
      write(iunit,'(a)') repeat("-", 195)
      write(iunit,'(  1x, a, ":" )') name
      if (verbose) then
         do i = 1, ndim
            write(iunit,'(*(f12.4,1x))') matrix(:,i)
         enddo
      else
         write(iunit,'(*(f12.4,1x))') matrix(:15,1)
      endif

      write(iunit,'(a)') repeat("-", 195)


   end subroutine print_full_matrix

end module matops