module gtb_lapack_eig
  use gtb_accuracy, only: ik, sp, dp
  implicit none
  private

  public :: la_syev, la_syevx, la_sygvx, la_sygvd

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
  end interface la_sygvd

end module gtb_lapack_eig
