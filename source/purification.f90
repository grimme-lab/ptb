module purification_
   use gtb_accuracy, only : wp
   use purification_settings, only : tPurificationSet
   use timer, only : tTimer
   implicit none

   public :: purification
   private
contains
   subroutine purification(ndim, H, S, P)
      
      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Hamiltonian matrix
      real(wp), intent(in) :: H(ndim*(ndim+1)/2)

      !> Overlap matrix
      real(wp), intent(in) :: S(ndim*(ndim+1)/2)

      !> Density matrix
      real(wp), intent(out) :: P(ndim*(ndim+1)/2)
   end subroutine purification

end module purification_