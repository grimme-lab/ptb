module purification_
   use ieee_arithmetic, only: ieee_is_NaN
   use iso_fortran_env, only: stdout => output_unit
   use gtb_accuracy, only : wp, ik
   use gtb_lapack_eig, only : la_syevd
   use gtb_la, only : la_gemm
   use matops, only : print_matrix
   use metrics, only: thrs, check_density
   use norms
   use purification_settings
   use timer, only : tTimer
   implicit none

   public :: purification
   private
   real(wp), allocatable :: identity(:,:)
   interface purification
      module procedure :: purification_wrapper
      module procedure :: purification_search
   end interface purification

contains
   subroutine purification_wrapper(pur, ndim, H, S, P)
      
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
      real(wp), dimension(ndim, ndim) :: Hmat, Smat, X
      integer :: i, j

      ! start timer !
      call timer_purification%new(2)
      
      ! Transformation matrix calculation !
      call timer_purification%click(1, 'transformation matrix')

      call blowsym(ndim, H, Hmat)  
      call blowsym(ndim, S, Smat)

      ! create identity matrix !
      allocate(identity(ndim,ndim), source=0.0_wp)
      do i=1,ndim
         identity(i,i) = 1.0_wp
      enddo

      X = get_transformation_matrix(pur, ndim, Smat) ! different powers of S
      call timer_purification%click(1)

      ! purification !
      call timer_purification%click(2, 'Purification procedure')
      call purification(pur, ndim, Hmat, Smat, X, P)      

      call timer_purification%click(2)

      deallocate(identity)
      call timer_purification%finalize('Total purification time')


   end subroutine purification_wrapper

   !> get different powers of S
   function get_transformation_matrix(pur, ndim, S) result(X)


       !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Overlap matrix
      real(wp), intent(in) :: S(ndim,ndim)

      !> Transformation matrix
      real(wp) :: X(ndim,ndim)

      !> Local variables
      logical :: debug  ! debugging mode
      real(wp), dimension(ndim, ndim) :: V, atmp, D ! buffer matrix
      real(wp) :: w(ndim) ! eigenvalues
      integer(ik) :: info ! execution status
      integer(ik) :: i
      real(wp), dimension(ndim) :: s_infinity, s_2 ! different norms
      real(wp), dimension(ndim,ndim) :: XsX, check ! buffer for iterations
      integer :: cycles, type

      debug = .false.
      V = S ! duplicate S
      D = 0.0_wp
      type = pur%metric%type
      cycles = pur%metric%cycles
      
      ! iterative solution !
      if (pur%metric%iterative) then
         selectcase(type)
         case(inv)

            ! Initial guess !
            X = S/((infinity(S))**2)
            if (debug) &
               & call print_matrix(ndim, X, 'Initial guess for itertive method')

            ! Newton Schulz iterations,  !
            do i = 1, cycles
               call la_gemm(S, X, atmp)
               call la_gemm(X, atmp, XsX)

               X = 2.0_wp * X - XsX
               call la_gemm(X, S, check)
               if (norm2(identity - check) < thrs%low) then
                  write(stdout,'(1x, a, 1x, I0)') 'Itertive inversion (s^-1) converged in:', i
                  exit
               endif
            enddo

         end select

      ! Lapack diagonalization !
      else

         call la_syevd(V, w, info)
         
         ! Construct X !
         select case(type)
         case(inv)
            do i = 1, ndim
               D(i,i) = 1.0_wp/w(i)
            enddo
         end select

         call la_gemm(V, D, atmp)
         call la_gemm(atmp, V, X, transb='T')
       
      endif

      ! if debug, check inversion !
      if (debug) then
         selectcase(type)
         case(inv)
            call la_gemm(X, S, atmp)
            call print_matrix(ndim, atmp, 'identity')
         endselect
      endif

   end function get_transformation_matrix


   subroutine purification_search(pur, ndim, Hmat, Smat, X, P)
 
      !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Hamiltonian matrix
      real(wp), intent(in) :: Hmat(ndim,ndim)

      !> Overlap matrix
      real(wp), intent(in) :: Smat(ndim,ndim)
      
      !> Transformation matrix
      real(wp), intent(in) :: X(ndim,ndim)

      !> Density matrix
      real(wp), intent(out) :: P(ndim*(ndim+1)/2)

      !> Locals
      real(wp) :: chempot ! chemical potential
      integer :: i ! counters
      integer :: cycles
      type(tTimer) :: timer_chempot_iteration
      real(wp), dimension(ndim,ndim) :: guess, purified
      integer :: type

      ! Intialization !
      chempot = pur%chempot%guess
      cycles = pur%chempot%cycles
      type = pur%type 
      call timer_chempot_iteration%new(2)

      ! Chemical potential search iterations !
      do i = 1, cycles
         
         ! Initial guess calculation !
         call timer_chempot_iteration%click(1, 'Inital guess')
         call build_guess(pur, ndim, Hmat, X, chempot, guess)
         call timer_chempot_iteration%click(1)
         
         call timer_chempot_iteration%click(2, 'Purification routine')
         select case(type)
         case(mcweeny)
            call mcweeny_(pur, ndim, guess, Smat, purified)
         endselect
         call timer_chempot_iteration%click(2)

      enddo
      call timer_chempot_iteration%finalize('Purification iterator')

   end subroutine purification_search

   subroutine build_guess(pur, ndim, Hmat, X, chempot, guess)

      !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Hamiltonian matrix
      real(wp), intent(in) :: Hmat(ndim,ndim)

      !> Transformation matrix
      real(wp), intent(in) :: X(ndim,ndim)

      !> Chemical potential
      real(wp), intent(in) :: chempot

      !> Guess
      real(wp), intent(out) :: guess(ndim, ndim)

      !> Locals
      integer :: type_
      real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, term5
      logical :: debug = .false.

      type_ = pur%type

      select case(type_)
      ! McWeeny purification !
      case(mcweeny)
         term1 = chempot * X !  μ*s^-1
         call la_gemm(Hmat, X, term2) ! H*s^-1
         call la_gemm(X, term2, term3) ! S^-1*H*S^-1
         term4 = term1 - term3  !  μ*s^-1 - S^-1*H*S^-1
         term5 = (1.0_wp/(min(frobenius(term4), gershgorin(term4)))) * term4 ! α(μ*s^-1 - S^-1*H*S^-1)
         guess = 0.5_wp * (term5 + X) 
         if(debug) &
            call print_matrix(ndim, guess, 'Guess')
      endselect
      
   end subroutine build_guess

   subroutine mcweeny_(pur, ndim, guess, Smat, P)

      !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Guess
      real(wp), intent(in) :: guess(ndim, ndim)

      !> Transformation matrix
      real(wp), intent(in) :: Smat(ndim,ndim)

      !> Purified Density 
      real(wp), intent(out) :: P(ndim,ndim)

      !> Locals
      integer :: cycles
      integer :: i
      real(wp), dimension(ndim,ndim) :: term1, term2, term3, guess2
      real(wp) :: norm
      logical :: debug
      
      debug = .false.
      cycles = pur%cycles
      guess2 = guess
      norm = frobenius(guess)

      ! Check variables !
      if (debug) then
         write(stdout, '(a,1x,a,1x,a)') repeat('!', 40),'DEBUG',repeat('!',40)
         write(stdout, '(1x, a, 1x, f14.7)'), "Norm: ", norm
         call print_matrix(ndim, guess2, 'Guess')
         call print_matrix(ndim, Smat, 'S')
         write(stdout, '(a)') repeat(':',87)
      endif

      do i = 1, 60

         call la_gemm(Smat, guess2, term1) ! S*P
         call la_gemm(guess2, term1, term2) ! P*S*P
         call la_gemm(term2, term1, term3) ! P*S*P*S*P

         P = 3 * term2 - 2 * term3 ! 3* P*S*P - * P*S*P*S*P
         if (debug) &
            call print_matrix(ndim, P*2, 'Purification P')

         ! check for NaN !
         if (any(ieee_is_NaN(P))) &
            error stop 'Error: NaN is encountered during purification'

         ! convergence !
         if (abs(norm2(P)-norm) < thrs%low4) then 
            write(stdout,'(a, 1x, i0, 1x, a)' ) 'McWeeny converged in:', i, 'cycles'
            exit
         endif
         guess2 = P
         norm = norm2(P)

      enddo
      P = P * 2

   end subroutine mcweeny_
end module purification_