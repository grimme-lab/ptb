module purification_
   use ieee_arithmetic, only: ieee_is_NaN
   use iso_fortran_env, only: stdout => output_unit
   use gtb_accuracy, only : wp, ik
   use gtb_lapack_eig, only : la_syevd
   use gtb_la, only : la_gemm
   use matops, only : print_matrix
   use metrics, only: thrs, get_nel
   use norms
   use purification_settings
   use metrics, only: thrs
   use timer, only : tTimer
   implicit none

   public :: purification
   private
   real(wp), allocatable :: identity(:,:)
   integer :: pr

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

      real(wp), dimension(ndim, ndim) :: Hmat, Smat, X, purified
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
      call timer_purification%click(2, 'Purification iterator')
      call purification(pur, ndim, Hmat, Smat, X, purified)    
      call packsym(ndim, purified, P)
      call timer_purification%click(2)

      call timer_purification%finalize('Total purification time')
      deallocate(identity)


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
      real(wp), dimension(ndim, ndim) :: V, atmp, D, atmp2 ! buffer matrix
      real(wp) :: w(ndim) ! eigenvalues
      integer(ik) :: info ! execution status
      integer(ik) :: i
      real(wp), dimension(ndim) :: s_infinity, s_2 ! different norms

      real(wp), dimension(ndim,ndim) :: XsX, check ! buffer for iterations
      integer :: cycles, type, pr

      
      debug = .false.
      V = S ! duplicate S
      D = 0.0_wp
      type = pur%metric%type
      cycles = pur%metric%cycles
      pr = pur%prlvl
      
      ! print section header !
      if (pr > 1) then
         write(stdout, '(/, a, /)') '--> get_transformation_matrix'
         select case(type)
         case(inv)
            write(stdout, '(a)') 'S^-1 Construction'
         case(inv_sqrt)
            write(stdout, '(a)') 'S^-0.5 Construction'
         endselect
      endif

      ! iterative solution !
      if (pur%metric%iterative) then

         selectcase(type)

         case(inv)

            ! Initial guess !
            X = S/((infinity(S))**2)
            if (debug) &
               & call print_matrix(ndim, X, 'Initial guess for itertive method')

            ! Newton-Schulz iterations !
            do i = 1, cycles
               call la_gemm(S, X, atmp, pr=pr)
               call la_gemm(X, atmp, XsX, pr=pr)

               X = 2.0_wp * X - XsX
               call la_gemm(X, S, check, pr=pr)
               if (norm2(identity - check) < thrs%low) then
                  if (pr > 0) &
                     write(stdout,'(1x, a, 1x, I0)') 'Itertive inversion (s^-1) converged in:', i
                  exit
               endif
            enddo

         case(inv_sqrt)

            ! Guess !
            X = S/((frobenius(S))**2)

            if (debug) &
               & call print_matrix(ndim, X, 'Iterative guess: ')
            
            do i = 1, cycles

               call la_gemm(X, X, atmp, pr=pr) ! X^2
               call la_gemm(atmp, S, XsX, pr=pr) ! X^2*S
               call la_gemm(X, XsX, atmp2, pr=pr) ! X*X^2*S
               
               X = 0.5_wp * (3 * X - atmp2) 

               if (debug) then
                  call print_matrix(ndim,  XsX, 'Check:')
               endif

               if (frobenius(identity - XsX) < thrs%low) then
                  if (pr > 0) &
                     write(stdout,'(1x, a, 1x, I0)') 'Itertive inversion (s^0.5) converged in:', i
                  exit
               elseif (any(ieee_is_NaN(XsX))) then
                  error stop 'Error, NaN encountered'
               endif

            enddo


         end select

      ! Lapack diagonalization !
      else

         call la_syevd(V, w, info, pr=pr)
         
         ! Construct X !
         select case(type)
         case(inv)
            do i = 1, ndim
               D(i,i) = 1.0_wp/w(i)
            enddo
         case(inv_sqrt)
            do i = 1, ndim
               D(i,i) = 1.0_wp/sqrt(w(i))
            enddo 
         end select

         call la_gemm(V, D, atmp, pr=pr)
         call la_gemm(atmp, V, X, transb='T', pr=pr)
       
      endif

      ! if debug, check inversion !
      if (debug) then
         write(stdout, '(a,1x,a,1x,a)') repeat(':', 20),'DEBUG',repeat(':',20)
         selectcase(type)
         case(inv)
            call la_gemm(X, S, atmp, pr=pr)
            call print_matrix(ndim, atmp, 'identity')
         case(inv_sqrt)
            call la_gemm(X, X, atmp, pr=pr)
            call la_gemm(S, atmp, atmp2, pr=pr)
            call print_matrix(ndim, atmp2, 'identity')
         endselect
         write(stdout, '(a)') repeat(':', 47)
      endif
      
      ! print trailer !
      if (pr > 1) &
         write(stdout, '(/, a, /)') '<-- get_transformation_matrix'
       
   end function get_transformation_matrix


   subroutine purification_search(pur, ndim, Hmat, Smat, X, purified)
 
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
      real(wp), intent(out) :: purified(ndim,ndim)

      !> Locals
      real(wp) :: chempot ! chemical potential
      real(wp) :: nel, nel_calc ! actual number of electrons
      integer :: i ! counters
      integer :: cycles
      type(tTimer) :: timer_chempot_iteration
      real(wp), dimension(ndim,ndim) :: guess
      integer :: type, pr
      logical :: debug
      
      real(wp) :: lower_bound, upper_bound, incr
      logical :: is_lower_bound, is_upper_bound, bisection

      ! Intialization !
      type = pur%type 
      nel = pur%nel 
      pr = pur%prlvl
      chempot = pur%chempot%guess
      cycles = pur%chempot%cycles
      incr = pur%chempot%increment 
      nel_calc = 0.0_wp

      debug = .true.
      bisection = .false.
      is_lower_bound = .false.
      is_upper_bound = .false.

      if (pr > 1) &
         write(stdout, '(a, /)') '--> purification_search'

      call timer_chempot_iteration%new(2)

      ! Chemical potential search iterations !
      search: do i = 1, cycles

         if (bisection) then
            if (debug) &
               write(stdout,'(a, 1x, f18.8, 1x, a, 1x, f18.8)' ) &
                  'Lower bound:', lower_bound, '| Upper bound:', upper_bound

            chempot = (lower_bound + upper_bound) / 2.0_wp
            if (abs(upper_bound-lower_bound) < thrs%high8) then
               write(stdout,'(a, 1x, i0, 1x, a)' ) &
                  'Chemical potential search. Too short/big interval to bisect.', i, 'cycles'
               exit
            endif
         endif      
         
         ! Initial guess calculation !
         call timer_chempot_iteration%click(1, 'Inital guess')
         call build_guess(pur, ndim, Hmat, X, chempot, guess)
         call timer_chempot_iteration%click(1)
         ! Run Purification !
         call timer_chempot_iteration%click(2, 'Purification routine')
         select case(type)
         case(mcweeny)
            call mcweeny_(pur, ndim, guess, Smat, purified)
         case(sign_iter_pade, sign_iter_newton)
            call sign_iterative(pur, ndim, guess, X, purified)
         case(sign_diagonalization)
            call sign_diagonalization_(pur, ndim, chempot, guess, X, purified)
         endselect
         call timer_chempot_iteration%click(2)
         
         nel_calc = get_nel(ndim, purified, Smat) 

         ! Printout !
         if (pr > 0 .and. type .ne. sign_diagonalization ) then
            write(stdout,'(a, 9x, i0, /, a, 1x, f18.8, /, a, 1x, f18.8)' ) &
               'Search number:             ', i, &
               'Current chemical potential:', chempot, &
               'Number of electrons:       ', nel_calc
         endif

         
         ! Exit condition: compare number of electrons !
         if (abs(nel_calc - real(nel)) < thrs%low4) then
            if (pr > 0 .and. type .ne. sign_diagonalization) &
               write(stdout,'(a, 1x, i0, 1x, a)' ) 'Chemical potential found after:', i, 'cycles'
            exit search
         endif
         
         ! adjust search boundaries !
         if (nel_calc < nel) then
            is_lower_bound = .true.
            lower_bound = chempot
            chempot = chempot + incr
         else                                 
            is_upper_bound = .true.
            upper_bound = chempot
            chempot = chempot - incr
         endif

         bisection = is_lower_bound .and. is_upper_bound
         incr = incr * 2.0_wp

      enddo search
      call timer_chempot_iteration%finalize('Purification iterator')

      if (pr > 1) &
         write(stdout, '(/, a, /)') '<-- purification_search'
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
      integer :: type_, pr, metric
      real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, term5
      logical :: debug = .true.

      type_ = pur%type
      metric = pur%metric%type
      pr = pur%prlvl

      if (pr > 1) &
         write(stdout, '(a, /)') '--> build_guess'

      select case(type_)

      ! McWeeny purification !
      case(mcweeny)
         
         if (pr > 1 .or. debug) &
            write(stdout, '(a, g0)') 'Mcweeny prufication with S^1 with chempot of ', chempot
         term1 = chempot * X !  μ*s^-1
         call la_gemm(Hmat, X, term2, pr=pr) ! H*s^-1
         call la_gemm(X, term2, term3, pr=pr) ! S^-1*H*S^-1
         term4 = term1 - term3  !  μ*s^-1 - S^-1*H*S^-1
         term5 = (1.0_wp/(min(frobenius(term4), gershgorin(term4)))) * term4 ! α(μ*s^-1 - S^-1*H*S^-1)
         guess = 0.5_wp * (term5 + X) 
      
      ! Sign Iteration !
      case(sign_iter_pade, sign_iter_newton, sign_diagonalization)

         select case(metric)
         ! S^-0.5 !
         case(inv_sqrt)

            if (pr > 1 .or. debug) &
               write(stdout, '(a, g0)') 'Sign purification with S^0.5 with chempot of ', chempot
            term1 = chempot * identity ! μ*I
            call la_gemm(Hmat, X, term2, pr=pr) ! H*S^-0.5 
            call la_gemm(X, term2, term3, pr=pr) ! S^-0.5*H*S^-0.5 
            term4 = term3 - term1
            guess = (1.0_wp/(min(frobenius(term4), gershgorin(term4)))) * term4 ! α(S^-0.5*H*S^-0.5 - μI)
         
         ! S^1 !
         case(inv)

            if (pr > 1 .or. debug) &
               write(stdout, '(a, g0)') 'Sign purification with S^-1 woth chempot of ', chempot
            term1 = chempot * identity
            call la_gemm(X, Hmat, term2, pr=pr)
            term3 = term2 - term1
            guess = (1.0_wp/(min(frobenius(term3), gershgorin(term3)))) * term3 ! α(S^-0.5*H*S^-0.5 - μI)

         endselect

      endselect

      if (debug) &
         call print_matrix(ndim, guess, 'Guess')
      
      if (pr > 1) &
         write(stdout, '(/, a, /)') '<-- build_guess'
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
      integer :: i, pr
      real(wp), dimension(ndim,ndim) :: term1, term2, term3, guess2
      real(wp) :: norm
      logical :: debug
      
      debug = .false.
      cycles = pur%cycles
      guess2 = guess
      norm = frobenius(guess)
      pr = pur%prlvl 

      if (pr > 1) &
         write(stdout, '(/, a, /)') '--> mcweeny_'

      ! Check variables !
      if (debug) then
         write(stdout, '(a,1x,a,1x,a)') repeat('!', 40),'DEBUG',repeat('!',40)
         write(stdout, '(1x, a, 1x, f14.7)'), "Norm: ", norm
         call print_matrix(ndim, guess2, 'Guess')
         call print_matrix(ndim, Smat, 'S')
         write(stdout, '(a)') repeat(':',87)
      endif

      do i = 1, cycles

         call la_gemm(Smat, guess2, term1, pr=pr) ! S*P
         call la_gemm(guess2, term1, term2, pr=pr) ! P*S*P
         call la_gemm(term2, term1, term3, pr=pr) ! P*S*P*S*P

         P = 3 * term2 - 2 * term3 ! 3* P*S*P - * P*S*P*S*P
         if (debug) &
            call print_matrix(ndim, P*2, 'Purification P')

         ! check for NaN !
         if (any(ieee_is_NaN(P))) &
            error stop 'Error: NaN is encountered during purification'

         ! convergence !
         if (abs(norm2(P)-norm) < thrs%low4) then 
            if (pr > 0) &
               write(stdout,'(a, 1x, i0, 1x, a)' ) 'McWeeny converged in:', i, 'cycles'
            exit
         endif
         guess2 = P
         norm = norm2(P)

      enddo

      P = P * 2

      if (pr > 1) &
         write(stdout, '(/, a, /)') '<-- mcweeny_'

   end subroutine mcweeny_

   subroutine sign_iterative(pur, ndim, guess, metric, purified)

      !> Purification settings
      type(tPurificationSet), intent(in) :: pur

      !> Number of BFs
      integer, intent(in) :: ndim

      !> Iteration matrix
      real(wp), intent(in) :: guess(ndim,ndim)

      !> Transformation matrix
      real(wp), intent(in) :: metric(ndim,ndim)

      !> Purified Density Matrix
      real(wp), intent(out) :: purified(ndim,ndim)

      !> Locals
      integer :: i, cycles, iter_type
      real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, X
      real(wp) :: norm
      logical :: debug = .true.

      cycles = pur%cycles
      pr = pur%prlvl
      iter_type = pur%type
      X = guess

      if (pr > 1) &
         write(stdout, '(a, /)') '--> sign_iterative'

      selectcase(iter_type)
      case(sign_iter_pade)
         ! pade iterations !
         do i = 1, cycles
            norm = norm2(X)         

            call la_gemm(X, X, term1, pr=pr) ! X^2
            call la_gemm(term1, term1, term2, pr=pr) ! X^4

            term3 = (15.0_wp * identity) - (10.0_wp * term1) + (3.0_wp * term2) ! 15 - 10*X^2 + 3*X^4 
            call la_gemm(X, term3, term4, alpha=0.125_wp, pr=pr)
            X = term4

            if (debug) &
               call print_matrix(ndim, X, 'X')

            if (any(ieee_is_NaN(X))) &
               error stop 'Error: NaN is encountered during purification'

            if (abs(norm2(X)-norm) < thrs%low4) then 
               if (pr > 0) &
                  write(stdout,'(/, a, 1x, i0, 1x, a)' ) 'Pade sign iterative converged in:', i, 'cycles'
               exit
            endif

         enddo
      case(sign_iter_newton)
         ! Newton Schulz iterations !
         do i = 1, cycles
            norm = norm2(X)         
            call la_gemm(X, X, term1, pr=pr) ! X^2
            call la_gemm(X, term1, term2, pr=pr) ! X^3
            term3 = 0.5_wp * (3 * X - term2) ! 1/2 * ( 3*X - X^3)
            X = term3
            if (debug) &
               call print_matrix(ndim, X, 'X')
            
            if (any(ieee_is_NaN(X))) &
               error stop 'Error: NaN is encountered during purification'

            if (abs(norm2(X)-norm) < thrs%low4) then 
               if (pr > 0) &
                  write(stdout,'(/, a, 1x, i0, 1x, a)' ) 'Newton-Schuz sign iterative converged in:', i, 'cycles'
               exit
            endif


         enddo
      endselect

      ! convert to P !
      call sign_to_density_matrix(pur, ndim, X, metric, purified)

      purified = purified * 2

      if (pr > 1) &
         write(stdout, '(/, a, /)') '<-- sign_iterative'

   end subroutine sign_iterative

   subroutine sign_to_density_matrix(pur, ndim, sign_, metric, P)

      !> Purification settings
      type(tPurificationSet), intent(in) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Sign matrix
      real(wp), intent(in) :: sign_(ndim, ndim)
      
      !> Transformation matrix
      real(wp), intent(in) :: metric(ndim, ndim)

      !> Density matrix
      real(wp), intent(out) :: P(ndim,ndim)

      !> Locals
      logical :: debug 
      integer :: metric_type, pr
      real(wp), dimension(ndim,ndim) :: tmp

      debug = .true.
      pr = pur%prlvl
      
      select case(pur%metric%type)
      case(inv_sqrt)
         if (pr > 1) &
            write(stdout, '(/, a, /)') 'Density matrix via S^-0.5'
         call la_gemm(identity-sign_, metric, tmp, pr=pr)
         call la_gemm(metric, tmp, P, alpha=0.5_wp, pr=pr)
      case(inv)
         if (pr > 1) &
            write(stdout, '(/, a, /)') 'Density matrix via S^-1'
         call la_gemm(identity-sign_, metric, P, alpha=0.5_wp, pr=pr)
      end select

      if (debug) then
         call print_matrix(ndim, P, 'Purified density')
      endif
 
      if (any(ieee_is_NaN(P))) &
         error stop 'Error: NaN is encountered during purification'
      
   end subroutine sign_to_density_matrix


   subroutine sign_diagonalization_(pur, ndim, chempot, guess, metric, purified)

      !> Purification settings
      type(tPurificationSet), intent(in) :: pur

      !> Number of BFs
      integer, intent(in) :: ndim

      !> Chemical potential
      real(wp), intent(inout) :: chempot

      !> Iteration matrix
      real(wp), intent(in) :: guess(ndim,ndim)

      !> Transformation matrix
      real(wp), intent(in) :: metric(ndim,ndim)

      !> Purified Density Matrix
      real(wp), intent(out) :: purified(ndim,ndim)

      !> Locals
      logical, save :: first_time =.true.
      real(wp), dimension(:, :), allocatable, save :: eigenvectors
      real(wp), dimension(:), allocatable, save :: eigenvalues
      integer(ik) :: pr, info, i, j, k, cycles
      real(wp) :: nelcalc, sum_, lower_bound, upper_bound, incr, nel
      real(wp), dimension(ndim) :: eigguess
      real(wp), dimension(ndim, ndim) :: eigval, sign_, tmp
      logical :: debug, is_lower_bound, is_upper_bound, bisection

      pr = pur%prlvl
      eigval = 0.0_wp
      nel = real(pur%nel)
      cycles = pur%chempot%cycles
      incr = pur%chempot%increment
      bisection = .false.; is_lower_bound = .false.; is_upper_bound = .false.
      debug =.true.

      if (pr > 1) &
         write(stdout, '(a, /)') '--> sign_diagonalization'

      ! diagonalize only once !
      if (first_time) then
         if (.not. allocated(eigenvalues)) allocate(eigenvectors(ndim,ndim), eigenvalues(ndim))
         eigenvectors = guess
         call la_syevd(eigenvectors, eigenvalues, info, pr=pr)
      endif

      
      readjust_chempot: do k = 1, cycles
         
         nelcalc = 0.0_wp
         eigguess = eigenvalues - chempot

         ! signum () !
         do i = 1, ndim

            if (abs(eigguess(i)) < thrs%normal) then
               eigval(i,i) = 0.0_wp
            else 
               eigval(i,i) = merge(1.0_wp, -1.0_wp, eigguess(i) > 0.0_wp)
            endif
            
            sum_ = 0.0_wp
            ! calculate number of electron !
            do j = 1, ndim
               sum_ = sum_ + eigenvectors(j,i) * eigenvectors(j,i)
            enddo
            nelcalc = nelcalc + (1.0_wp - (sum_ * eigval(i,i)))
         enddo

         ! Printout !
         if (pr > 0) then
            write(stdout,'(a, 9x, i0, /, a, 1x, f18.8, /, a, 1x, g0, /, a, 1x, g0)' ) &
               'Search number:             ', k, &
               'Current chemical potential:', chempot, &
               'Current increment:         ', incr, &
               'Number of electrons:       ', nelcalc
         endif

         
         ! Exit condition: compare number of electrons !
         if (abs(nelcalc - nel) < thrs%normal) then
            if (pr > 0) &
               write(stdout,'(a, 1x, i0, 1x, a)' ) 'Internal chemical potential found after:', k, 'cycles'
            exit readjust_chempot
         endif

         ! readjust chempot internally: bisection !
         if (nelcalc < nel) then
            is_lower_bound = .true.
            lower_bound = chempot
            chempot = chempot + incr
         else                                 
            is_upper_bound = .true.
            upper_bound = chempot
            chempot = chempot - incr
         endif

         bisection = is_lower_bound .and. is_upper_bound
         incr = incr * 2.0_wp

         if (bisection) then
            if (debug) &
               write(stdout,'(a, 1x, f18.8, 1x, a, 1x, f18.8)' ) &
                  'Lower bound:', lower_bound, '| Upper bound:', upper_bound

            chempot = (lower_bound + upper_bound) / 2.0_wp
            if (abs(upper_bound-lower_bound) < thrs%high8) then
               write(stdout,'(a, 1x, i0, 1x, a)' ) &
                  'Internal chemical potential search. Too short/big interval to bisect.', i, 'cycles'
               exit readjust_chempot
            endif
         endif      
 
      enddo readjust_chempot


      call la_gemm(eigval, eigenvectors, tmp, transb='T', pr=pr) ! signum(Λ) * Q^T
      call la_gemm(eigenvectors, tmp, sign_, pr=pr) ! !Q * signum(Λ) * Q^T

      if (debug) &
         call print_matrix(ndim, sign_, 'Sign')

      call sign_to_density_matrix(pur, ndim, sign_, metric, purified)
      
      purified = purified * 2

      if (pr > 1) &
         write(stdout, '(a, /)') '<-- sign_diagonalization'


   end subroutine sign_diagonalization_


end module purification_