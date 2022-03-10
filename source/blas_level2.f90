!> Interfaces to level 2 BLAS routines.
module gtb_la_level2
  use gtb_accuracy, only : ik, sp, dp
  implicit none
  private

  public :: la_gbmv, la_gemv, la_ger, la_gerc, la_geru, la_sbmv, la_spmv, &
    & la_hbmv, la_hemv, la_spr2, la_spr, la_syr2, la_syr, la_her2, la_her, &
    & la_symv, la_hpmv, la_hpr2, la_hpr, la_tbmv, la_tbsv, la_tpmv, la_tpsv, &
    & la_trmv, la_trsv

  !> Performs one of the matrix-vector operations
  !>
  !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are vectors and A is an
  !> m by n band matrix, with kl sub-diagonals and ku super-diagonals.
  interface la_gbmv
    pure subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sgbmv
    pure subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dgbmv

    pure subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine cgbmv
    pure subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: kl
      integer(ik), intent(in) :: ku
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zgbmv

    module procedure :: la_gbmv_rsp
    module procedure :: la_gbmv_csp
    module procedure :: la_gbmv_rdp
    module procedure :: la_gbmv_cdp
  end interface la_gbmv

  !> Performs one of the matrix-vector operations
  !>
  !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are vectors and A is an
  !> m by n matrix.
  interface la_gemv
    pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sgemv
    pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dgemv

    pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine cgemv
    pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: trans
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zgemv

    module procedure :: la_gemv_rsp
    module procedure :: la_gemv_csp
    module procedure :: la_gemv_rdp
    module procedure :: la_gemv_cdp
  end interface la_gemv


  !> Performs the rank 1 operation
  !>
  !>    A := alpha*x*y**T + A,
  !>
  !> where alpha is a scalar, x is an m element vector, y is an n element
  !> vector and A is an m by n matrix.
  interface la_ger
    pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(inout) :: a(lda, *)
    end subroutine sger
    pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(inout) :: a(lda, *)
    end subroutine dger

    module procedure :: la_ger_rsp
    module procedure :: la_ger_rdp
  end interface la_ger

  !> Performs the rank 1 operation
  !>
  !>    A := alpha*x*y**H + A,
  !>
  !> where alpha is a scalar, x is an m element vector, y is an n element
  !> vector and A is an m by n matrix.
  interface la_gerc
    pure subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine cgerc
    pure subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine zgerc

    module procedure :: la_gerc_csp
    module procedure :: la_gerc_cdp
  end interface la_gerc

  !> Performs the rank 1 operation
  !>
  !>    A := alpha*x*y**T + A,
  !>
  !> where alpha is a scalar, x is an m element vector, y is an n element
  !> vector and A is an m by n matrix.
  interface la_geru
    pure subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine cgeru
    pure subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine zgeru

    module procedure :: la_geru_csp
    module procedure :: la_geru_cdp
  end interface la_geru

  !> Performs the matrix-vector  operation
  !>
  !>    y := alpha*A*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are n element vectors and
  !> A is an n by n symmetric band matrix, with k super-diagonals.
  interface la_sbmv
    pure subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine ssbmv
    pure subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dsbmv

    module procedure :: la_sbmv_rsp
    module procedure :: la_sbmv_rdp
  end interface la_sbmv

  !> Performs the matrix-vector operation
  !>
  !>    y := alpha*A*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are n element vectors and
  !> A is an n by n symmetric matrix, supplied in packed form.
  interface la_spmv
    pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: ap(*)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sspmv
    pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: ap(*)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dspmv

    module procedure :: la_spmv_rsp
    module procedure :: la_spmv_rdp
  end interface la_spmv

  !> Performs the matrix-vector  operation
  !>
  !>    y := alpha*A*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are n element vectors and
  !> A is an n by n hermitian band matrix, with k super-diagonals.
  interface la_hbmv
    pure subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine chbmv
    pure subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zhbmv

    module procedure :: la_hbmv_csp
    module procedure :: la_hbmv_cdp
  end interface la_hbmv

  !> Performs the matrix-vector  operation
  !>
  !>    y := alpha*A*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are n element vectors and
  !> A is an n by n hermitian matrix.
  interface la_hemv
    pure subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine chemv
    pure subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zhemv

    module procedure :: la_hemv_csp
    module procedure :: la_hemv_cdp
  end interface la_hemv

  !> Performs the symmetric rank 2 operation
  !>
  !>    A := alpha*x*y**T + alpha*y*x**T + A,
  !>
  !> where alpha is a scalar, x and y are n element vectors and A is an
  !> n by n symmetric matrix, supplied in packed form.
  interface la_spr2
    pure subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(inout) :: ap(*)
    end subroutine sspr2
    pure subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(inout) :: ap(*)
    end subroutine dspr2

    module procedure :: la_spr2_rsp
    module procedure :: la_spr2_rdp
  end interface la_spr2

  !> Performs the symmetric rank 1 operation
  !>
  !>    A := alpha*x*x**T + A,
  !>
  !> where alpha is a real scalar, x is an n element vector and A is an
  !> n by n symmetric matrix, supplied in packed form.
  interface la_spr
    pure subroutine sspr(uplo, n, alpha, x, incx, ap)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: ap(*)
    end subroutine sspr
    pure subroutine dspr(uplo, n, alpha, x, incx, ap)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: ap(*)
    end subroutine dspr

    module procedure :: la_spr_rsp
    module procedure :: la_spr_rdp
  end interface la_spr

  !> Performs the symmetric rank 2 operation
  !>
  !>    A := alpha*x*y**T + alpha*y*x**T + A,
  !>
  !> where alpha is a scalar, x and y are n element vectors and A is an n
  !> by n symmetric matrix.
  interface la_syr2
    pure subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(inout) :: a(lda, *)
    end subroutine ssyr2
    pure subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(inout) :: a(lda, *)
    end subroutine dsyr2

    module procedure :: la_syr2_rsp
    module procedure :: la_syr2_rdp
  end interface la_syr2

  !> Performs the symmetric rank 1 operation
  !>
  !>    A := alpha*x*x**T + A,
  !>
  !> where alpha is a real scalar, x is an n element vector and A is an
  !> n by n symmetric matrix.
  interface la_syr
    pure subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: a(lda, *)
    end subroutine ssyr
    pure subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: a(lda, *)
    end subroutine dsyr

    module procedure :: la_syr_rsp
    module procedure :: la_syr_rdp
  end interface la_syr

  !> Performs the hermitian rank 2 operation
  !>
  !>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
  !>
  !> where alpha is a scalar, x and y are n element vectors and A is an n
  !> by n hermitian matrix.
  interface la_her2
    pure subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine cher2
    pure subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine zher2

    module procedure :: la_her2_csp
    module procedure :: la_her2_cdp
  end interface la_her2

  !> Performs the hermitian rank 1 operation
  !>
  !>    A := alpha*x*x**H + A,
  !>
  !> where alpha is a real scalar, x is an n element vector and A is an
  !> n by n hermitian matrix.
  interface la_her
    pure subroutine cher(uplo, n, alpha, x, incx, a, lda)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine cher
    pure subroutine zher(uplo, n, alpha, x, incx, a, lda)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: a(lda, *)
    end subroutine zher

    module procedure :: la_her_csp
    module procedure :: la_her_cdp
  end interface la_her

  !> Performs the matrix-vector  operation
  !>
  !>    y := alpha*A*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are n element vectors and
  !> A is an n by n symmetric matrix.
  interface la_symv
    pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine ssymv
    pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: beta
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dsymv

    module procedure :: la_symv_rsp
    module procedure :: la_symv_rdp
  end interface la_symv

  !> Performs the matrix-vector operation
  !>
  !>    y := alpha*A*x + beta*y,
  !>
  !> where alpha and beta are scalars, x and y are n element vectors and
  !> A is an n by n hermitian matrix, supplied in packed form.
  interface la_hpmv
    pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: ap(*)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine chpmv
    pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: ap(*)
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: beta
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zhpmv

    module procedure :: la_hpmv_csp
    module procedure :: la_hpmv_cdp
  end interface la_hpmv

  !> Performs the hermitian rank 2 operation
  !>
  !>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
  !>
  !> where alpha is a scalar, x and y are n element vectors and A is an
  !> n by n hermitian matrix, supplied in packed form.
  interface la_hpr2
    pure subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: ap(*)
    end subroutine chpr2
    pure subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp), intent(inout) :: ap(*)
    end subroutine zhpr2

    module procedure :: la_hpr2_csp
    module procedure :: la_hpr2_cdp
  end interface la_hpr2

  !> Performs the hermitian rank 1 operation
  !>
  !>    A := alpha*x*x**H + A,
  !>
  !> where alpha is a real scalar, x is an n element vector and A is an
  !> n by n hermitian matrix, supplied in packed form.
  interface la_hpr
    pure subroutine chpr(uplo, n, alpha, x, incx, ap)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: ap(*)
    end subroutine chpr
    pure subroutine zhpr(uplo, n, alpha, x, incx, ap)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: alpha
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: ap(*)
    end subroutine zhpr

    module procedure :: la_hpr_csp
    module procedure :: la_hpr_cdp
  end interface la_hpr

  !> Performs one of the matrix-vector operations
  !>
  !>    x := A*x,   or   x := A**T*x,
  !>
  !> where x is an n element vector and  A is an n by n unit, or non-unit,
  !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.
  interface la_tbmv
    pure subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stbmv
    pure subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtbmv

    pure subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctbmv
    pure subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztbmv

    module procedure :: la_tbmv_rsp
    module procedure :: la_tbmv_csp
    module procedure :: la_tbmv_rdp
    module procedure :: la_tbmv_cdp
  end interface la_tbmv

  !> Solves one of the systems of equations
  !>
  !>    A*x = b,   or   A**T*x = b,
  !>
  !> where b and x are n element vectors and A is an n by n unit, or
  !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
  !> diagonals.
  !>
  !> No test for singularity or near-singularity is included in this
  !> routine. Such tests must be performed before calling this routine.
  interface la_tbsv
    pure subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stbsv
    pure subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtbsv

    pure subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctbsv
    pure subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: k
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztbsv

    module procedure :: la_tbsv_rsp
    module procedure :: la_tbsv_csp
    module procedure :: la_tbsv_rdp
    module procedure :: la_tbsv_cdp
  end interface la_tbsv

  !> Performs one of the matrix-vector operations
  !>
  !>    x := A*x,   or   x := A**T*x,
  !>
  !> where x is an n element vector and  A is an n by n unit, or non-unit,
  !> upper or lower triangular matrix, supplied in packed form.
  interface la_tpmv
    pure subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: ap(*)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stpmv
    pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: ap(*)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtpmv

    pure subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: ap(*)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctpmv
    pure subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: ap(*)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztpmv

    module procedure :: la_tpmv_rsp
    module procedure :: la_tpmv_csp
    module procedure :: la_tpmv_rdp
    module procedure :: la_tpmv_cdp
  end interface la_tpmv

  !> Solves one of the systems of equations
  !>
  !>    A*x = b,   or   A**T*x = b,
  !>
  !> where b and x are n element vectors and A is an n by n unit, or
  !> non-unit, upper or lower triangular matrix, supplied in packed form.
  !>
  !> No test for singularity or near-singularity is included in this
  !> routine. Such tests must be performed before calling this routine.
  interface la_tpsv
    pure subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: ap(*)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine stpsv
    pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: ap(*)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtpsv

    pure subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: ap(*)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctpsv
    pure subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: ap(*)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztpsv

    module procedure :: la_tpsv_rsp
    module procedure :: la_tpsv_csp
    module procedure :: la_tpsv_rdp
    module procedure :: la_tpsv_cdp
  end interface la_tpsv

  !> Performs one of the matrix-vector operations
  !
  !>    x := A*x,   or   x := A**T*x,
  !>
  !> where x is an n element vector and  A is an n by n unit, or non-unit,
  !> upper or lower triangular matrix.
  interface la_trmv
    pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine strmv
    pure subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtrmv

    pure subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctrmv
    pure subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztrmv

    module procedure :: la_trmv_rsp
    module procedure :: la_trmv_csp
    module procedure :: la_trmv_rdp
    module procedure :: la_trmv_cdp
  end interface la_trmv

  !> Solves one of the systems of equations
  !>
  !>    A*x = b,   or   A**T*x = b,
  !>
  !> where b and x are n element vectors and A is an n by n unit, or
  !> non-unit, upper or lower triangular matrix.
  !>
  !> No test for singularity or near-singularity is included in this
  !> routine. Such tests must be performed before calling this routine.
  interface la_trsv
    pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine strsv
    pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dtrsv

    pure subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ctrsv
    pure subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: lda
      character, intent(in) :: uplo
      character, intent(in) :: trans
      character, intent(in) :: diag
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine ztrsv

    module procedure :: la_trsv_rsp
    module procedure :: la_trsv_csp
    module procedure :: la_trsv_rdp
    module procedure :: la_trsv_cdp
  end interface la_trsv


contains

  pure subroutine la_gbmv_rsp(amat, xvec, yvec, kl, m, alpha, beta, trans)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    integer, intent(in), optional :: kl
    integer, intent(in), optional :: m
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    integer(ik) :: n, ku, lda
    integer(ik) :: kl_, m_
    real(wp) :: a, b
    character :: tra

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    kl_ = (lda-1)/2
    if (present(kl)) kl_ = kl
    m_ = n
    if (present(m)) m_ = m
    ku = lda-kl_-1
    call la_gbmv(tra, m_, n, kl_, ku, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gbmv_rsp

  pure subroutine la_gbmv_csp(amat, xvec, yvec, kl, m, alpha, beta, trans)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    integer, intent(in), optional :: kl
    integer, intent(in), optional :: m
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    integer(ik) :: n, ku, lda
    integer(ik) :: kl_, m_
    complex(wp) :: a, b
    character :: tra

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    kl_ = (lda-1)/2
    if (present(kl)) kl_ = kl
    m_ = n
    if (present(m)) m_ = m
    ku = lda-kl_-1
    call la_gbmv(tra, m_, n, kl_, ku, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gbmv_csp

  pure subroutine la_gbmv_rdp(amat, xvec, yvec, kl, m, alpha, beta, trans)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    integer, intent(in), optional :: kl
    integer, intent(in), optional :: m
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    integer(ik) :: n, ku, lda
    integer(ik) :: kl_, m_
    real(wp) :: a, b
    character :: tra

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    kl_ = (lda-1)/2
    if (present(kl)) kl_ = kl
    m_ = n
    if (present(m)) m_ = m
    ku = lda-kl_-1
    call la_gbmv(tra, m_, n, kl_, ku, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gbmv_rdp

  pure subroutine la_gbmv_cdp(amat, xvec, yvec, kl, m, alpha, beta, trans)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    integer, intent(in), optional :: kl
    integer, intent(in), optional :: m
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    integer(ik) :: n, ku, lda
    integer(ik) :: kl_, m_
    complex(wp) :: a, b
    character :: tra

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    kl_ = (lda-1)/2
    if (present(kl)) kl_ = kl
    m_ = n
    if (present(m)) m_ = m
    ku = lda-kl_-1
    call la_gbmv(tra, m_, n, kl_, ku, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gbmv_cdp



  pure subroutine la_gemv_rsp(amat, xvec, yvec, alpha, beta, trans)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    real(wp) :: a, b
    character :: tra
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_gemv(tra, m, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gemv_rsp

  pure subroutine la_gemv_csp(amat, xvec, yvec, alpha, beta, trans)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    complex(wp) :: a, b
    character :: tra
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_gemv(tra, m, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gemv_csp

  pure subroutine la_gemv_rdp(amat, xvec, yvec, alpha, beta, trans)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    real(wp) :: a, b
    character :: tra
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_gemv(tra, m, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gemv_rdp

  pure subroutine la_gemv_cdp(amat, xvec, yvec, alpha, beta, trans)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta
    character, intent(in), optional :: trans

    complex(wp) :: a, b
    character :: tra
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    tra = 'n'
    if (present(trans)) tra = trans
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_gemv(tra, m, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_gemv_cdp


  pure subroutine la_ger_rsp(amat, xvec, yvec, alpha)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(inout) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(in) :: yvec(:)
    real(wp), intent(in), optional :: alpha

    real(wp) :: a
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_ger(m, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_ger_rsp

  pure subroutine la_gerc_csp(amat, xvec, yvec, alpha)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(inout) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp), intent(in), optional :: alpha

    complex(wp) :: a
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_gerc(m, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_gerc_csp

  pure subroutine la_geru_csp(amat, xvec, yvec, alpha)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(inout) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp), intent(in), optional :: alpha

    complex(wp) :: a
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_geru(m, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_geru_csp

  pure subroutine la_ger_rdp(amat, xvec, yvec, alpha)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(inout) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(in) :: yvec(:)
    real(wp), intent(in), optional :: alpha

    real(wp) :: a
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_ger(m, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_ger_rdp

  pure subroutine la_gerc_cdp(amat, xvec, yvec, alpha)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(inout) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp), intent(in), optional :: alpha

    complex(wp) :: a
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_gerc(m, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_gerc_cdp

  pure subroutine la_geru_cdp(amat, xvec, yvec, alpha)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(inout) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp), intent(in), optional :: alpha

    complex(wp) :: a
    integer(ik) :: m, n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call la_geru(m, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_geru_cdp


  pure subroutine la_sbmv_rsp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula
    real(wp) :: a, b
    integer(ik) :: n, k, lda

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(xvec)
    call la_sbmv(ula, n, k, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_sbmv_rsp

  pure subroutine la_hbmv_csp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula
    complex(wp) :: a, b
    integer(ik) :: n, k, lda

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(xvec)
    call la_hbmv(ula, n, k, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_hbmv_csp

  pure subroutine la_sbmv_rdp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula
    real(wp) :: a, b
    integer(ik) :: n, k, lda

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(xvec)
    call la_sbmv(ula, n, k, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_sbmv_rdp

  pure subroutine la_hbmv_cdp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula
    complex(wp) :: a, b
    integer(ik) :: n, k, lda

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(xvec)
    call la_hbmv(ula, n, k, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_hbmv_cdp


  pure subroutine la_spmv_rsp(apmat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: apmat(:)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula
    real(wp) :: a, b
    integer(ik) :: n

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_spmv(ula, n, a, apmat, xvec, 1, b, yvec, 1)
  end subroutine la_spmv_rsp

  pure subroutine la_hpmv_csp(apmat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: apmat(:)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula
    complex(wp) :: a, b
    integer(ik) :: n

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_hpmv(ula, n, a, apmat, xvec, 1, b, yvec, 1)
  end subroutine la_hpmv_csp

  pure subroutine la_spmv_rdp(apmat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: apmat(:)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula
    real(wp) :: a, b
    integer(ik) :: n

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_spmv(ula, n, a, apmat, xvec, 1, b, yvec, 1)
  end subroutine la_spmv_rdp

  pure subroutine la_hpmv_cdp(apmat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: apmat(:)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula
    complex(wp) :: a, b
    integer(ik) :: n

    a = 1
    if (present(alpha)) a = alpha
    b = 0
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_hpmv(ula, n, a, apmat, xvec, 1, b, yvec, 1)
  end subroutine la_hpmv_cdp


  pure subroutine la_symv_rsp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula
    real(wp) :: a, b
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_symv(ula, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_symv_rsp

  pure subroutine la_hemv_csp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula
    complex(wp) :: a, b
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_hemv(ula, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_hemv_csp

  pure subroutine la_symv_rdp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: amat(:, :)
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta

    character :: ula
    real(wp) :: a, b
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_symv(ula, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_symv_rdp

  pure subroutine la_hemv_cdp(amat, xvec, yvec, uplo, alpha, beta)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: amat(:, :)
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp), intent(in), optional :: beta

    character :: ula
    complex(wp) :: a, b
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    b = 0.0_wp
    if (present(beta)) b = beta
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_hemv(ula, n, a, amat, lda, xvec, 1, b, yvec, 1)
  end subroutine la_hemv_cdp


  pure subroutine la_spr2_rsp(apmat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = sp
    real(wp), intent(inout) :: apmat(:)
    real(wp), intent(in) :: xvec(:)
    real(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_spr2(ula, n, a, xvec, 1, yvec, 1, apmat)
  end subroutine la_spr2_rsp

  pure subroutine la_hpr2_csp(apmat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = sp
    complex(wp), intent(inout) :: apmat(:)
    complex(wp), intent(in) :: xvec(:)
    complex(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha

    character :: ula
    complex(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_hpr2(ula, n, a, xvec, 1, yvec, 1, apmat)
  end subroutine la_hpr2_csp

  pure subroutine la_spr2_rdp(apmat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = dp
    real(wp), intent(inout) :: apmat(:)
    real(wp), intent(in) :: xvec(:)
    real(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_spr2(ula, n, a, xvec, 1, yvec, 1, apmat)
  end subroutine la_spr2_rdp

  pure subroutine la_hpr2_cdp(apmat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = dp
    complex(wp), intent(inout) :: apmat(:)
    complex(wp), intent(in) :: xvec(:)
    complex(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha

    character :: ula
    complex(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_hpr2(ula, n, a, xvec, 1, yvec, 1, apmat)
  end subroutine la_hpr2_cdp


  pure subroutine la_spr_rsp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = sp
    real(wp), intent(inout) :: amat(:)
    real(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_spr(ula, n, a, xvec, 1, amat)
  end subroutine la_spr_rsp

  pure subroutine la_hpr_csp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = sp
    complex(wp), intent(inout) :: amat(:)
    complex(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_hpr(ula, n, a, xvec, 1, amat)
  end subroutine la_hpr_csp

  pure subroutine la_spr_rdp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = dp
    real(wp), intent(inout) :: amat(:)
    real(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_spr(ula, n, a, xvec, 1, amat)
  end subroutine la_spr_rdp

  pure subroutine la_hpr_cdp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = dp
    complex(wp), intent(inout) :: amat(:)
    complex(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_hpr(ula, n, a, xvec, 1, amat)
  end subroutine la_hpr_cdp


  pure subroutine la_syr_rsp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = sp
    real(wp), intent(inout) :: amat(:, :)
    real(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_syr(ula, n, a, xvec, 1, amat, lda)
  end subroutine la_syr_rsp

  pure subroutine la_her_csp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = sp
    complex(wp), intent(inout) :: amat(:, :)
    complex(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_her(ula, n, a, xvec, 1, amat, lda)
  end subroutine la_her_csp

  pure subroutine la_syr_rdp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = dp
    real(wp), intent(inout) :: amat(:, :)
    real(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_syr(ula, n, a, xvec, 1, amat, lda)
  end subroutine la_syr_rdp

  pure subroutine la_her_cdp(amat, xvec, uplo, alpha)
    integer, parameter :: wp = dp
    complex(wp), intent(inout) :: amat(:, :)
    complex(wp), intent(in) :: xvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_her(ula, n, a, xvec, 1, amat, lda)
  end subroutine la_her_cdp


  pure subroutine la_syr2_rsp(amat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = sp
    real(wp), intent(inout) :: amat(:, :)
    real(wp), intent(in) :: xvec(:)
    real(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_syr2(ula, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_syr2_rsp

  pure subroutine la_her2_csp(amat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = sp
    complex(wp), intent(inout) :: amat(:, :)
    complex(wp), intent(in) :: xvec(:)
    complex(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha

    character :: ula
    complex(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_her2(ula, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_her2_csp

  pure subroutine la_syr2_rdp(amat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = dp
    real(wp), intent(inout) :: amat(:, :)
    real(wp), intent(in) :: xvec(:)
    real(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    real(wp), intent(in), optional :: alpha

    character :: ula
    real(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_syr2(ula, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_syr2_rdp

  pure subroutine la_her2_cdp(amat, xvec, yvec, uplo, alpha)
    integer, parameter :: wp = dp
    complex(wp), intent(inout) :: amat(:, :)
    complex(wp), intent(in) :: xvec(:)
    complex(wp), intent(in) :: yvec(:)
    character, intent(in), optional :: uplo
    complex(wp), intent(in), optional :: alpha

    character :: ula
    complex(wp) :: a
    integer(ik) :: n, lda

    a = 1.0_wp
    if (present(alpha)) a = alpha
    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_her2(ula, n, a, xvec, 1, yvec, 1, amat, lda)
  end subroutine la_her2_cdp


  pure subroutine la_tbmv_rsp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbmv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbmv_rsp

  pure subroutine la_tbmv_csp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbmv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbmv_csp

  pure subroutine la_tbmv_rdp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbmv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbmv_rdp

  pure subroutine la_tbmv_cdp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbmv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbmv_cdp


  pure subroutine la_tbsv_rsp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbsv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbsv_rsp

  pure subroutine la_tbsv_csp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbsv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbsv_csp

  pure subroutine la_tbsv_rdp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbsv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbsv_rdp

  pure subroutine la_tbsv_cdp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, k, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    k = size(amat, 1)-1
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_tbsv(ula, tra, dia, n, k, amat, lda, xvec, 1)
  end subroutine la_tbsv_cdp


  pure subroutine la_tpmv_rsp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: apmat(:)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpmv(ula, tra, dia, n, apmat, xvec, 1)
  end subroutine la_tpmv_rsp

  pure subroutine la_tpmv_csp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: apmat(:)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpmv(ula, tra, dia, n, apmat, xvec,  1)
  end subroutine la_tpmv_csp

  pure subroutine la_tpmv_rdp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: apmat(:)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpmv(ula, tra, dia, n, apmat, xvec, 1)
  end subroutine la_tpmv_rdp

  pure subroutine la_tpmv_cdp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: apmat(:)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpmv(ula, tra, dia, n, apmat, xvec,  1)
  end subroutine la_tpmv_cdp


  pure subroutine la_tpsv_rsp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: apmat(:)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpsv(ula, tra, dia, n, apmat, xvec, 1)
  end subroutine la_tpsv_rsp

  pure subroutine la_tpsv_csp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: apmat(:)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpsv(ula, tra, dia, n, apmat, xvec, 1)
  end subroutine la_tpsv_csp

  pure subroutine la_tpsv_rdp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: apmat(:)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpsv(ula, tra, dia, n, apmat, xvec, 1)
  end subroutine la_tpsv_rdp

  pure subroutine la_tpsv_cdp(apmat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: apmat(:)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    n = size(xvec)
    call la_tpsv(ula, tra, dia, n, apmat, xvec, 1)
  end subroutine la_tpsv_cdp


  pure subroutine la_trmv_rsp(amat,  xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trmv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trmv_rsp

  pure subroutine la_trmv_csp(amat,  xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trmv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trmv_csp

  pure subroutine la_trmv_rdp(amat,  xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trmv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trmv_rdp

  pure subroutine la_trmv_cdp(amat,  xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trmv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trmv_cdp


  pure subroutine la_trsv_rsp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trsv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trsv_rsp

  pure subroutine la_trsv_csp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trsv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trsv_csp

  pure subroutine la_trsv_rdp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: amat(:, :)
    real(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trsv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trsv_rdp

  pure subroutine la_trsv_cdp(amat, xvec, uplo, trans, diag)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: amat(:, :)
    complex(wp), intent(inout) :: xvec(:)
    character, intent(in), optional :: uplo
    character, intent(in), optional :: trans
    character, intent(in), optional :: diag

    character :: ula, tra, dia
    integer(ik) :: n, lda

    dia = 'n'
    tra = 'n'
    ula = 'u'
    if (present(diag)) dia = diag
    if (present(trans)) tra = trans
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call la_trsv(ula, tra, dia, n, amat, lda, xvec, 1)
  end subroutine la_trsv_cdp


end module gtb_la_level2
