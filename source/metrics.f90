module metrics
   use iso_fortran_env, only: stdout => output_unit
   use gtb_accuracy, only : wp
   use gtb_la, only: la_gemm, la_symm
   use matops, only: trace, print_matrix
   implicit none
   private
   public :: get_nel, thrs, idempotent, check_density, analyze_results

   !> Type for storing all thresholds
   type :: tThreshold
      real(wp) :: high8
      real(wp) :: normal
      real(wp) :: low
      real(wp) :: low4
      real(wp) :: low2
      real(wp) :: low1
   endtype tThreshold

   type(tThreshold), parameter :: thrs = &
      & tThreshold(high8=1.0e-8_wp, normal = 1.0e-7_wp, low = 1.0e-6_wp, low4 = 1.0e-4_wp, low1=0.1, &
      & low2=1.0e-2_wp)

   interface band_structure
      module procedure :: get_band_structure_E_full
      module procedure :: get_band_structure_E_packed
   end interface band_structure

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
      if (abs(get_nel(ndim,P,S)) - nel > thrs%normal) &
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

      !> Locals
      logical :: debug = .false.

      if (debug) &
         write(stdout,'(a, 1x, i0, a, 1x, f18.8, 2x)') &
            'nel:', nel, ', nel_calc:', get_nel(ndim,P,S)

      ! Nel check !
      if (abs(get_nel(ndim,P,S)) - nel > thrs%low2) then
         error stop "Error: density matrix gives wrong N_el."
      endif
      
      ! Idempotency check! 
      if (.not. idempotent(ndim,P,S)) &
         error stop "Erorr: density matrix is not idempotent."

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


   real(wp) function get_band_structure_E_full(ndim, P, H) result(E)

      !> Number of basis functions
      integer, intent(in) :: ndim
      
      !> Density Matrix
      real(wp), intent(in) :: P(ndim,ndim)
      
      !> Hamiltonian matrix
      real(wp), intent(in) :: H(ndim,ndim) 
      
      !> Product of P*S
      real(wp) :: PH(ndim,ndim) 
      

      call la_gemm(P, H, PH) 
      E = trace(PH)

   end function get_band_structure_E_full

   real(wp) function get_band_structure_E_packed(ndim, Psym, Hsym) result(E)

      implicit none

      !> Number of basis functions
      integer, intent(in) :: ndim
      
      !> Density Matrix
      real(wp), intent(in) :: Psym(ndim*(ndim+1)/2)
      
      !> Hamiltonian matrix
      real(wp), intent(in) :: Hsym(ndim*(ndim+1)/2) 
      
      !> Locals
      real(wp) :: H(ndim,ndim) 
      real(wp) :: P(ndim,ndim)
      real(wp) :: PH(ndim,ndim) 

      call blowsym(ndim, Psym, P)
      call blowsym(ndim, Hsym, H)
      
      call la_gemm(P, H, PH) 
      E = trace(PH)

   end function get_band_structure_E_packed


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
      idempot = merge(.true., .false., frob_diff < thrs%normal)
      
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
      idempot = merge(.true., .false., frob_diff < thrs%low2)
         
   end function idempotent_packed

   subroutine analyze_results(ndim, P, P2, S, H, n)

       !> Number of BFs 
      integer, intent(in) :: ndim
      
      !> PTB P
      real(wp), intent(in) :: P(ndim*(ndim+1)/2)
      
      !> Purified P
      real(wp), intent(in) :: P2(ndim*(ndim+1)/2)

      !> Hamiltonian matrix
      real(wp), intent(in) :: H(ndim*(ndim+1)/2)

      !> Overlap matrix
      real(wp), intent(in) :: S(ndim*(ndim+1)/2)

      !> Number of atoms
      integer, intent(in) :: n

      ! ΔN_el !
      write(stdout,'(a,1x,f18.12)') &
         'ΔNel =   ', abs(get_nel(ndim, P, S) - get_nel(ndim, P2, S))
      
      write(stdout,'(a,1x,f18.12)') &
         'ΔE_band =', (abs(band_structure(ndim, P, H) - band_structure(ndim, P2, H)))/real(N)

   end subroutine analyze_results

end module metrics