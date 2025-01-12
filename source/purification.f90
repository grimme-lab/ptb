module purification_
   use gtb_accuracy, only : wp
   use purification_settings, only : tPurificationSet
   use timer, only : tTimer
   implicit none

   public :: purification
   private
   real(wp), allocatable :: identity(:,:)

contains
   subroutine purification(pur, ndim, H, S, P)
      
      !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Hamiltonian matrix
      real(wp), intent(in) :: H(ndim*(ndim+1)/2)

      !> Overlap matrix
      real(wp), intent(in) :: S(ndim*(ndim+1)/2)

      !> Density matrix
      real(wp), intent(out) :: P(ndim*(ndim+1)/2)

      !> Local variables
      type(tTimer) :: timer_purification
      real(wp), dimension(ndim, ndim) :: Hmat, Smat, metric
      integer :: i, j

      ! start timer !
      call timer_purification%new(1)
      call timer_purification%click(1, 'placeholder')

      call blowsym(ndim, H, Hmat)  
      call blowsym(ndim, S, Smat)

      ! create identity matrix !
      allocate(identity(ndim,ndim), source=0.0_wp)
      do i=1,ndim
         identity(i,i) = 1.0_wp
      enddo

      metric = get_transformation_matrix(pur, ndim, Smat) ! different powers of S
      call timer_purification%click(1)

      deallocate(identity)
      call timer_purification%finalize('Total purification time')


   end subroutine purification

   function get_transformation_matrix(pur, ndim, S) result(X)

       !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Overlap matrix
      real(wp), intent(in) :: S(ndim,ndim)

      !> Transformation matrix
      real(wp) :: X(ndim,ndim)

   end function get_transformation_matrix

end module purification_