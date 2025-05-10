module metrics
   use iso_fortran_env, only: stdout => output_unit
   use gtb_accuracy, only : dp, sp,wp
   use gtb_la, only: la_gemm, la_symm
   use matops, only: trace, print_matrix
   implicit none
   private
   public :: get_nel, thrs, idempotent, check_density, analyze_results
   public :: adjust_thresholds, band_structure, check_sparsity

   !> Type for storing all thresholds
   type :: tThreshold
      real(wp) :: high
      real(wp) :: normal
      real(wp) :: low
   endtype tThreshold

   type(tThreshold), allocatable :: thrs

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

   interface adjust_thresholds
      module procedure :: adjust_thresholds_sp
      module procedure :: adjust_thresholds_dp
   end interface adjust_thresholds
   
contains
   !> adjust thresholds dynamically
   subroutine adjust_thresholds_sp(ndim, Hmat)


      !> Calculation domain
      integer, intent(in) :: ndim

      !> mat
      real(sp), intent(in) :: Hmat(ndim*(ndim+1)/2)

      if (.not. allocated(thrs)) then
         allocate(thrs)
         if (ndim < 300) then
            thrs%high = 1e-7
            thrs%normal = 1e-6
            thrs%low = 1e-5
         elseif (ndim < 2000) then
            thrs%high = 1e-5
            thrs%normal = 1e-4
            thrs%low = 1e-3
         elseif (ndim < 15000) then
            thrs%high = 1e-4
            thrs%normal = 1e-3
            thrs%low = 1e-2
         else 
            thrs%high = 1e-2
            thrs%normal = 1e-1
            thrs%low = 1e-0
         endif
         
      endif


   end subroutine adjust_thresholds_sp

   !> adjust thresholds dynamically
   subroutine adjust_thresholds_dp(ndim, Hmat)

      !> Calculation domain
      integer, intent(in) :: ndim

      !> mat
      real(dp), intent(in) :: Hmat(ndim*(ndim+1)/2)

      if (.not. allocated(thrs)) then
         allocate(thrs)
         if (ndim < 300) then
            thrs%high = 1e-10
            thrs%normal = 1e-9
            thrs%low = 1e-8
         elseif (ndim < 2000) then
            thrs%high = 1e-7
            thrs%normal = 1e-6
            thrs%low = 1e-5
         elseif (ndim < 15000) then
            thrs%high = 1e-5
            thrs%normal = 1e-4
            thrs%low = 1e-3
         else 
            thrs%high = 1e-3
            thrs%normal = 1e-2
            thrs%low = 1e-1
         endif
         
      endif


   end subroutine adjust_thresholds_dp

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
      logical :: debug = .true.
      real(wp) :: nelP

      if (debug) &
         write(stdout,'(a, 1x, i0, a, 1x, f18.8, 2x)') &
            'nel:', nel, ', nel_calc:', get_nel(ndim,P,S)

      ! Nel check !
      nelP=get_nel(ndim,P,S)
      print*,"nel in P=",nelP,"required nel=",nel
      if (abs(nelP - nel) > thrs%normal) then
         error stop "Error: density matrix gives wrong N_el."
      endif
      
      ! Idempotency check! 
      if (.not. idempotent(ndim,P,S)) &
         error stop "Error: density matrix is not idempotent."

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
      print*,"idempotency deviation",frob_diff,thrs%normal
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
      print*,"idempotency deviation",frob_diff,thrs%normal
      idempot = merge(.true., .false., frob_diff < thrs%low)
         
   end function idempotent_packed

   subroutine analyze_results(ndim, P, P2, S, H, n, pr)

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

      !> printlevel
      integer, intent(in) :: pr

      !> Number of atoms
      integer, intent(in) :: n

      !> Locals
      real(wp) :: d_nel, d_band_str
      logical :: debug

      debug = .false.

      if (pr > 1) &
         write(stdout, '(/,a,/)') 'Enter: analyze_results' 

      if (debug) then
         call print_matrix(ndim, P, 'PTB Density matrix')
         call print_matrix(ndim, P2, 'Purified Density matrix')
      endif
      d_nel =  abs(get_nel(ndim, P, S) - get_nel(ndim, P2, S))
      d_band_str =(abs(band_structure(ndim, P, H) - band_structure(ndim, P2, H)))/real(N)

      write(stdout,'(/)')
      
      ! ΔN_el !
      write(stdout,*) &
         'ΔNel =   ', d_nel
      
      write(stdout,*) &
         'ΔE_band =', d_band_str 
      
      write(stdout,'(/)')

      if (pr > 1) &
         write(stdout, '(/, a, /)') 'Exit: analyze_results' 

   end subroutine analyze_results


   subroutine check_sparsity(pr, ndim, mat, sparse)
      
      !> special truncation parameter for sparsity
      real(wp), parameter :: trunc = 1e-5_wp

      integer, intent(in) :: pr
      integer, intent(in) :: ndim
      real(wp), intent(in) :: mat(ndim,ndim)
      logical, intent(out), optional :: sparse

      !> locals
      integer :: nzeros
      real(wp) :: ratio


      ! check general sparsity !
      nzeros = count(abs(mat) < trunc)
      ratio = real(nzeros)/real(ndim*ndim)
      if (pr > 0) &
         write(stdout, '(a, 2x, f5.2, a)') 'Sparisty:',ratio * 100, '%'
      
      if (present(sparse)) &
         sparse = ratio > 0.9_wp
      
   end subroutine check_sparsity

end module metrics
