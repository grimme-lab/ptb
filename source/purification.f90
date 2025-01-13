module purification_
   use gtb_accuracy, only : wp, ik
   use gtb_lapack_eig, only : la_syevd
   use gtb_la, only : la_gemm
   use matops, only : print_matrix
   use metrics, only: thrs
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

   !> get different powers of S
   function get_transformation_matrix(pur, ndim, S) result(X)

      use purification_settings

       !> Purification settings
      type(tPurificationSet) :: pur

      !> Number of basis functions
      integer, intent(in) :: ndim

      !> Overlap matrix
      real(wp), intent(in) :: S(ndim,ndim)

      !> Transformation matrix
      real(wp) :: X(ndim,ndim)

      !> Local variables
      logical :: debug = .true. ! debugging mode
      real(wp), dimension(ndim, ndim) :: V, atmp, D ! buffer matrix
      real(wp) :: w(ndim) ! eigenvalues
      integer(ik) :: info ! execution status
      integer(ik) :: i
      real(wp), dimension(ndim) :: s_infinity, s_2 ! different norms

      V = S ! duplicate S
      D = 0.0_wp
      
      ! iterative solution !
      if (pur%metric_num) then
         selectcase(pur%metric)
         case(inv)
            ! caculate scaling !
            do i = 1, ndim
               S_2(i) = sum(abs(S(:,i))) 
               S_infinity(i) = sum(abs(S(i,:))) 
            enddo
            print*, s_2
            print*, s_infinity
            call print_matrix(ndim, s, 'S')
            X= (1.0_wp / (maxval(s_2) * maxval(s_infinity))) * S
            
            ! Newton Schulz iterations !
            do i = 1, 10

               call la_gemm(X, X, atmp)
               call la_gemm(atmp, S, atmp)

               X = 2.0_wp * X - atmp
               if (norm2(identity - matmul(X,s)) < thrs%general) exit
            enddo

         end select

      ! Lapack diagonalization !
      else
         call la_syevd(V, w, info)
         
         ! Construct X !
         select case(pur%metric)
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
         selectcase(pur%metric)
         case(inv)
            call la_gemm(X, S, atmp)
            call print_matrix(ndim, atmp, 'identity')
         endselect
      endif
      stop 'inverse'


   end function get_transformation_matrix

end module purification_