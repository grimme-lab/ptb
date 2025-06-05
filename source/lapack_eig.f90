module gtb_lapack_eig
!   use accel_lib
!   use cuda_, only: ctx
   use gtb_accuracy, only: ik, sp, dp, i4
   use iso_fortran_env, only: stdout =>  output_unit
  implicit none
  private

  public :: la_syev, la_syevx, la_sygvx, la_sygvd, la_syevd

  !> Computes all eigenvalues and, optionally, eigenvectors of a
  !> real symmetric matrix A.
  interface la_syev
    pure subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: w(*)
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
    end subroutine ssyev
    pure subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: w(*)
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
    end subroutine dsyev
  end interface la_syev

  !> Computes selected eigenvalues and, optionally, eigenvectors
  !> of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
  !> selected by specifying either a range of values or a range of indices
  !> for the desired eigenvalues.
  interface la_syevx
    pure subroutine ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,&
        & z, ldz, work, lwork, iwork, ifail, info)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(wp), intent(out) :: z(ldz, *)
      real(wp), intent(in) :: vl
      real(wp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(wp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine ssyevx
    pure subroutine dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,&
        & z, ldz, work, lwork, iwork, ifail, info)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(wp), intent(out) :: z(ldz, *)
      real(wp), intent(in) :: vl
      real(wp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(wp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine dsyevx
  end interface la_syevx

  !> Computes selected eigenvalues, and optionally, eigenvectors
  !> of a real generalized symmetric-definite eigenproblem, of the form
  !> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
  !> and B are assumed to be symmetric and B is also positive definite.
  !> Eigenvalues and eigenvectors can be selected by specifying either a
  !> range of values or a range of indices for the desired eigenvalues.
  interface la_sygvx
    pure subroutine ssygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il,  &
        & iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(wp), intent(out) :: z(ldz, *)
      real(wp), intent(in) :: vl
      real(wp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(wp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine ssygvx
    pure subroutine dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il,  &
        & iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(wp), intent(out) :: z(ldz, *)
      real(wp), intent(in) :: vl
      real(wp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(wp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine dsygvx
  end interface la_sygvx

  !> Computes all the eigenvalues, and optionally, the eigenvectors
  !> of a real generalized symmetric-definite eigenproblem, of the form
  !> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
  !> B are assumed to be symmetric and B is also positive definite.
  !> If eigenvectors are desired, it uses a divide and conquer algorithm.
  interface la_sygvd
    pure subroutine ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork,    &
        & iwork, liwork, info)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine ssygvd
    pure subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork,    &
        & iwork, liwork, info)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine dsygvd
    module procedure :: la_sygvd_rdp
    module procedure :: la_sygvd_rsp
  end interface la_sygvd

  
   !> Computes all eigenvalues and, optionally, eigenvectors of a
   !> real symmetric matrix A with divide and conquer algorithm
   interface la_syevd
      pure subroutine ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
         import :: ik, sp
         integer, parameter :: wp = sp
         real(wp), intent(inout) :: a(lda, *)
         real(wp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(wp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
         integer(ik), intent(inout) :: iwork(*)
         integer(ik), intent(in) :: liwork
      end subroutine ssyevd
      pure subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
         import :: ik, dp
         integer, parameter :: wp = dp
         real(wp), intent(inout) :: a(lda, *)
         real(wp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(wp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
         integer(ik), intent(inout) :: iwork(*)
         integer(ik), intent(in) :: liwork
      end subroutine dsyevd
      module procedure :: la_syevd_rdp
      module procedure :: la_syevd_rsp
   end interface la_syevd
contains
   subroutine la_sygvd_rsp(a, b, w, info, itype, jobz, uplo, pr)
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(:,:)
      real(wp), intent(inout) :: b(:,:)
      real(wp), intent(out) :: w(:)
      integer(ik), intent(out) :: info
      integer(ik), intent(in), optional :: itype
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(ik), intent(in), optional :: pr

      character(len=1) :: job, upl
      integer(ik) :: n, lwork, liwork, ityp

      !> workspace
      real(wp), allocatable :: work(:)
      integer(ik), allocatable :: iwork(:)
      logical :: show
      integer(i4) :: err

      if (present(pr)) then
         show = pr > 1
      else
         show = .false.
      endif

      !  eigenvalue problem type !
      ityp = 1
      if (present(itype)) ityp =  itype

      ! if eigenvalues or eigenvalues + eigenvectors!
      job = 'V'
      if (present(jobz)) job = jobz

      ! upper or lower triangle to use !
      upl = 'U'
      if (present(uplo)) upl = uplo
      n = size(a, 2)

!      if (allocated(ctx)) then
!         if (show) &
!            write(stdout, '(3x, a)') 'Lapack: cuda_ssygvd'
!         call cuda_ssygvd(ctx, n, a, b, w, err)
!
!         if (err /= 0) &
!            error stop 'Error: cuda_ssygvd failed'
!      else
         if (show) &
            write(stdout, '(3x, a)') 'Lapack: ssygvd'

         lwork = -1
         liwork = -1
         w = 0.0_wp
         allocate(work(1))
         allocate(iwork(1))
         call la_sygvd(ityp, job, upl, n, a, n, b, n, w, work, lwork, iwork, liwork, info)
         
         lwork = int(work(1))
         liwork = iwork(1)
         deallocate(work, iwork)
         allocate(work(lwork))
         allocate(iwork(liwork))
         call la_sygvd(ityp, job, upl, n, a, n, b, n, w, work, lwork, iwork, liwork, info)
         deallocate(work, iwork)
!      endif
   
   end subroutine la_sygvd_rsp

   subroutine la_sygvd_rdp(a, b, w, info, itype, jobz, uplo, pr)
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(:,:)
      real(wp), intent(inout) :: b(:,:)
      real(wp), intent(out) :: w(:)
      integer(ik), intent(out) :: info
      integer(ik), intent(in), optional :: itype
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(ik), intent(in), optional :: pr

      character(len=1) :: job, upl
      integer :: n, lwork, liwork, ityp

      !> workspace
      real(wp), allocatable :: work(:)
      integer(ik), allocatable :: iwork(:)
      logical :: show
      integer(i4) :: err

      if (present(pr)) then
         show = pr > 1
      else
         show = .false.
      endif
      ! show = merge( pr > 1, .false., present(pr))

      !  eigenvalue problem type !
      ityp = 1
      if (present(itype)) ityp =  itype
      
      ! if eigenvalues or eigenvalues + eigenvectors!
      job = 'V'
      if (present(jobz)) job = jobz

      ! upper or lower triangle to use !
      upl = 'U'
      if (present(uplo)) upl = uplo
      n = size(a, 2)

!      if (allocated(ctx)) then
!         if (show) &
!            write(stdout, '(3x, a)') 'Lapack: cuda_dsygvd'
!         call cuda_dsygvd(ctx, n, a, b, w, err)
!
!         if (err /= 0) &
!            error stop 'Error: cuda_dsygvd failed'
!      else
         if (show) &
            write(stdout, '(3x, a)') 'Lapack: dsygvd'

         lwork = -1
         liwork = -1
         w = 0.0_wp
         allocate(work(1))
         allocate(iwork(1))
         call la_sygvd(ityp, job, upl, n, a, n, b, n, w, work, lwork, iwork, liwork, info)
         
         lwork = idint(work(1))
         liwork = iwork(1)
         deallocate(work, iwork)
         allocate(work(lwork))
         allocate(iwork(liwork))
         call la_sygvd(ityp, job, upl, n, a, n, b, n, w, work, lwork, iwork, liwork, info)
         deallocate(work, iwork)
      
!      endif

   end subroutine la_sygvd_rdp



   subroutine la_syevd_rsp(a, w, info, jobz, uplo, pr)
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(:,:)
      real(wp), intent(out) :: w(:)
      integer(ik), intent(out) :: info
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(ik), intent(in), optional :: pr

      character(len=1) :: job, upl
      integer(ik) :: n, lwork, liwork

      !> workspace
      real(wp), allocatable :: work(:)
      integer(ik), allocatable :: iwork(:)
      integer(i4) :: err
      logical :: show

      if (present(pr)) then
         show = pr > 1
      else
         show = .false.
      endif

      ! if eigenvalues or eigenvalues + eigenvectors!
      job = 'V'
      if (present(jobz)) job = jobz

      ! upper or lower triangle to use !
      upl = 'U'
      if (present(uplo)) upl = uplo
      n = size(a, 2)

!      if (allocated(ctx)) then
!         if (show) &
!            write(stdout, '(3x, a)') 'Lapack: cuda_ssyevd'
!         call cuda_ssyevd(ctx, n, a, w, err)
!
!         if (err /= 0) &
!            error stop 'Error: cuda_ssyevd failed.'
!      else
         if (show) &
            write(stdout, '(3x, a)') 'Lapack: ssyevd'

         lwork = -1
         liwork = -1
         w = 0.0_wp
         allocate(work(1))
         allocate(iwork(1))
         call la_syevd(job, upl, n, a, n, w, work, lwork, iwork, liwork, info)
         
         lwork = int(work(1))
         liwork = iwork(1)
         deallocate(work, iwork)
         allocate(work(lwork))
         allocate(iwork(liwork))
         call la_syevd(job, upl, n, a, n, w, work, lwork, iwork, liwork, info)
         deallocate(work, iwork)
         
!      endif

   end subroutine la_syevd_rsp

   subroutine la_syevd_rdp(a, w, info, jobz, uplo, pr)
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(:,:)
      real(wp), intent(out) :: w(:)
      integer(ik), intent(out) :: info
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(ik), intent(in), optional :: pr

      character(len=1) :: job, upl
      integer :: n, lwork, liwork

      !> workspace
      real(wp), allocatable :: work(:)
      integer(ik), allocatable :: iwork(:)
      integer(i4) :: err
      logical :: show

      if (present(pr)) then
         show = pr > 1
      else
         show = .false.
      endif

      ! if eigenvalues or eigenvalues + eigenvectors!
      job = 'V'
      if (present(jobz)) job = jobz

      ! upper or lower triangle to use !
      upl = 'U'
      if (present(uplo)) upl = uplo
      n = size(a, 2)

!      if (allocated(ctx)) then
!         if (show) &
!            write(stdout, '(3x, a)') 'Lapack: cuda_dsyevd'
!         call cuda_dsyevd(ctx, n, a, w, err)
!
!         if (err /= 0) &
!            error stop 'Error: cuda_dsyevd failed.'
!      else
         if (show) &
            write(stdout, '(3x, a)') 'Lapack: dsyevd'

         ! inquery !
         lwork = -1
         liwork = -1
         w = 0.0_wp
         allocate(work(1))
         allocate(iwork(1))
         call la_syevd(job, upl, n, a, n, w, work, lwork, iwork, liwork, info)
         
         lwork = idint(work(1))
         !print*,"lwork",lwork
         liwork = iwork(1)
         deallocate(work, iwork)
         allocate(work(lwork))
         allocate(iwork(liwork))
         call la_syevd(job, upl, n, a, n, w, work, lwork, iwork, liwork, info)
         deallocate(work, iwork)
!      endif

   end subroutine la_syevd_rdp

end module gtb_lapack_eig
