module metrics
   use gtb_accuracy, only : wp
   use gtb_la, only: la_gemm
   use matops, only: trace
   implicit none
   private
   public :: get_nel, thrs

   !> Type for storing all thresholds
   type :: tThreshold
      real(wp) :: general
   endtype tThreshold

   type(tThreshold), parameter :: thrs = &
      & tThreshold(general = 1.0e-8_wp)

   interface get_nel
      module procedure :: get_nel_full
      module procedure :: get_nel_packed
   end interface get_nel
contains

!> Calculate number of electrons
real(wp) function get_nel_full(ndim, P, S) result(nel)

   implicit none

   !> Number of basis functions
   integer, intent(in) :: ndim
   
   !> Overlap matrix
   real(wp), intent(in) :: S(ndim,ndim) 
   
   !> Density Matrix
   real(wp), intent(in) :: P(ndim,ndim)
   
   !> Product of P*S
   real(wp) :: PS(ndim,ndim) 
   
   call la_gemm(P, S, PS) 
   nel = trace(PS)

end function get_nel_full

!> Calculate number of electrons (from packed sym matrices)
real(wp) function get_nel_packed(ndim, Psym, Ssym) result(nel)

   implicit none

   !> Number of basis functions
   integer, intent(in) :: ndim
   
   !> Overlap matrix
   real(wp), intent(in) :: Ssym(ndim*(ndim+1)/2) 
   
   !> Density Matrix
   real(wp), intent(in) :: Psym(ndim*(ndim+1)/2)
   
   !> Overlap matrix
   real(wp) :: S(ndim,ndim) 
   
   !> Density Matrix
   real(wp) :: P(ndim,ndim)
   
   !> Product of P*S
   real(wp) :: PS(ndim,ndim) 

   call blowsym(ndim, Psym, P)
   call blowsym(ndim, Ssym, S)
   
   call la_gemm(P, S, PS) 
   nel = trace(PS)

end function get_nel_packed

end module metrics