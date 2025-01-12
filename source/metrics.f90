module metrics
   use gtb_accuracy, only : wp
   use gtb_la, only: la_gemm, la_symm
   use matops, only: trace, print_matrix
   implicit none
   private
   public :: get_nel, thrs, idempotent, check_density

   !> Type for storing all thresholds
   type :: tThreshold
      real(wp) :: general
   endtype tThreshold

   type(tThreshold), parameter :: thrs = &
      & tThreshold(general = 1.0e-8_wp)

   interface check_density
      module procedure :: check_density_full
      module procedure :: check_density_packed
   end interface check_density

   interface get_nel
      module procedure :: get_nel_full
      module procedure :: get_nel_packed
   end interface get_nel

   interface idempotent
      module procedure :: idempotent_full
      module procedure :: idempotent_packed
   end interface idempotent
   
contains

   !> Check if the density matrix is valid
   subroutine check_density_full(ndim, P, S, nel)

      !> Number of basis functions
      integer, intent(in) :: ndim
      
      !> Density Matrix
      real(wp), intent(in) :: P(ndim,ndim)
      
      !> Overlap matrix
      real(wp), intent(in) :: S(ndim,ndim) 
      
      !> Number of electrons
      integer, intent(in) :: nel

      ! Nel check !
      if (abs(get_nel(ndim,P,S)) - nel > thrs%general) &
         error stop "wrong Nel, check P"
      
      ! Idempotency check! 
      if (.not. idempotent(ndim,P,S)) &
         error stop "P is not idempotent"

   end subroutine check_density_full


   !> Check if the density matrix is valid
   subroutine check_density_packed(ndim, P, S, nel)

      !> Number of basis functions
      integer, intent(in) :: ndim
      
      !> Density Matrix
      real(wp), intent(in) :: P(ndim*(ndim+1)/2)
      
      !> Overlap matrix
      real(wp), intent(in) :: S(ndim*(ndim+1)/2) 
      
      !> Number of electrons
      integer, intent(in) :: nel

      ! Nel check !
      if (abs(get_nel(ndim,P,S)) - nel > thrs%general) &
         error stop "wrong Nel, check P"
      
      ! Idempotency check! 
      if (.not. idempotent(ndim,P,S)) &
         error stop "P is not idempotent"

   end subroutine check_density_packed

   !> Calculate number of electrons
   real(wp) function get_nel_full(ndim, P, S) result(nel)

      !> Number of basis functions
      integer, intent(in) :: ndim
      
      !> Density Matrix
      real(wp), intent(in) :: P(ndim,ndim)
      
      !> Overlap matrix
      real(wp), intent(in) :: S(ndim,ndim) 
      
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
      
      !> Density Matrix
      real(wp), intent(in) :: Psym(ndim*(ndim+1)/2)
      
      !> Overlap matrix
      real(wp), intent(in) :: Ssym(ndim*(ndim+1)/2) 
      
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

   !> Idempotency check for full matrices
   logical function idempotent_full(ndim, mat, S) result(idempot)
      
      !> Number of BFs 
      integer, intent(in) :: ndim
      
      !> Input matrix (P)
      real(wp), intent(in) :: mat(ndim,ndim)
      
      !> Overlap (for orthogonal case)
      real(wp), intent(in), optional :: S(ndim,ndim)

      !> Local variables
      real(wp), dimension(ndim,ndim) :: matmat, mm2
      real(wp) :: frob_diff
      
      matmat = 0.0_wp 
      mm2 = 0.0_wp

      ! P*S*P !
      if (present(S)) then
         call la_symm(mat, S, matmat)
         call la_symm(matmat, mat, mm2, alpha=0.5_wp)
      ! P*P !
      else
         call la_symm(mat, mat, matmat)
      endif

      ! check the divergence !
      frob_diff = sqrt(sum(mm2-mat)**2)
      idempot = merge(.true., .false., frob_diff < thrs%general)
      
   end function idempotent_full

   !> idempotency check for packed matrices
   logical function idempotent_packed(ndim, mat, S) result(idempot)
      
      !> Number of BFs 
      integer, intent(in) :: ndim
      
      !> Input matrix (P)
      real(wp), intent(in) :: mat(ndim*(ndim+1)/2)
      
      !> Overlap (for orthogonal case)
      real(wp), intent(in), optional :: S(ndim*(ndim+1)/2)

      !> Buffer variables
      real(wp), dimension(ndim, ndim) :: Sdum, matblowed, matmat, mm2
      real(wp) :: frob_diff
      real(wp), external :: dlange
      real(wp), allocatable :: work(:)


      allocate(work(ndim))
      
      if(present(S)) &
         & call blowsym(ndim, S, Sdum)
      call blowsym(ndim, mat, matblowed)
      matmat=0.0_wp 
      
      ! dim check !
      if (size(matmat, 2) /= size(matblowed, 1)) then
         error stop "Error: Dimensions do not match for multiplication."
      endif
      
      ! P*S*P !
      if (present(S)) then
         call la_symm(matblowed, Sdum, matmat)
         call la_gemm(matmat, matblowed, mm2, alpha=0.5_wp)
         
       ! P*P !
      else
         call la_symm(matblowed, matblowed, mm2)
      endif
      
      ! check the divergence !
      frob_diff = sqrt(sum(mm2-matblowed)**2)
      idempot = merge(.true., .false., frob_diff < thrs%general)
         
   end function idempotent_packed

end module metrics