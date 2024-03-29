!> Interfaces to level 1 BLAS routines.
module gtb_la_level1
  use gtb_accuracy, only : ik, sp, dp
  implicit none
  private

  public :: la_axpy, la_copy, la_dot, la_dotc, la_dotu, la_dsdot, la_sdsdot, &
    & la_rotg, la_rotm, la_rotmg, la_scal, la_rot, la_swap, la_abs1, &
    & la_asum, la_nrm2, la_iamax

  !> Constant times a vector plus a vector.
  interface la_axpy
    pure subroutine saxpy(n, a, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine saxpy
    pure subroutine daxpy(n, a, x, incx, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine daxpy

    pure subroutine caxpy(n, a, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine caxpy
    pure subroutine zaxpy(n, a, x, incx, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zaxpy

    module procedure :: la_axpy_rsp
    module procedure :: la_axpy_csp
    module procedure :: la_axpy_rdp
    module procedure :: la_axpy_cdp
  end interface la_axpy

  !> Copies a vector, x, to a vector, y.
  interface la_copy
    pure subroutine scopy(n, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine scopy
    pure subroutine dcopy(n, x, incx, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dcopy

    pure subroutine ccopy(n, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine ccopy
    pure subroutine zcopy(n, x, incx, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zcopy

    module procedure :: la_copy_rsp
    module procedure :: la_copy_csp
    module procedure :: la_copy_rdp
    module procedure :: la_copy_cdp
  end interface la_copy

  !> Forms the dot product of two vectors.
  interface la_dot
    pure function sdot(n, x, incx, y, incy) result(dot)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp) :: dot
    end function sdot
    pure function ddot(n, x, incx, y, incy) result(dot)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp) :: dot
    end function ddot

    module procedure :: la_dot_rsp
    module procedure :: la_dot_rdp
  end interface la_dot

  !> Compute the inner product of two vectors with extended
  !> precision accumulation and result.
  !>
  !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
  !> DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
  !> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !> defined in a similar way using INCY.
  interface la_dsdot
    pure function dsdot(n, sx, incx, sy, incy)
      import :: ik, dp, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: sy(*)
      integer(ik), intent(in) :: incy
      real(dp) :: dsdot
    end function dsdot

    module procedure :: la_dsdot_rspdp
  end interface la_dsdot

  !> Compute the inner product of two vectors with extended
  !> precision accumulation.
  !>
  !> Returns S.P. result with dot product accumulated in D.P.
  !> SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
  !> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !> defined in a similar way using INCY.
  interface la_sdsdot
    pure function sdsdot(n, sb, sx, incx, sy, incy)
      import :: ik, sp
      integer(ik), intent(in) :: n
      real(sp), intent(in) :: sb
      real(sp), intent(in) :: sx(*)
      integer(ik), intent(in) :: incx
      real(sp), intent(in) :: sy(*)
      integer(ik), intent(in) :: incy
      real(sp) :: sdsdot
    end function sdsdot

    module procedure :: la_sdsdot_rsp
  end interface la_sdsdot

  !> Forms the dot product of two complex vectors
  !>
  !>    DOTC = X^H * Y
  interface la_dotc
    pure function cdotc(n, x, incx, y, incy) result(dotc)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp) :: dotc
    end function cdotc
    pure function zdotc(n, x, incx, y, incy) result(dotc)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp) :: dotc
    end function zdotc

    module procedure :: la_dotc_csp
    module procedure :: la_dotc_cdp
  end interface la_dotc

  !> Forms the dot product of two complex vectors
  !>
  !>    DOTU = X^T * Y
  interface la_dotu
    pure function cdotu(n, x, incx, y, incy) result(dotu)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp) :: dotu
    end function cdotu
    pure function zdotu(n, x, incx, y, incy) result(dotu)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(in) :: y(*)
      integer(ik), intent(in) :: incy
      complex(wp) :: dotu
    end function zdotu

    module procedure :: la_dotu_csp
    module procedure :: la_dotu_cdp
  end interface la_dotu

  !> The computation uses the formulas
  !>    sigma = sgn(a)    if |a| >  |b|
  !>          = sgn(b)    if |b| >= |a|
  !>    r = sigma*sqrt( a**2 + b**2 )
  !>    c = 1; s = 0      if r = 0
  !>    c = a/r; s = b/r  if r != 0
  !> The subroutine also computes
  !>    z = s    if |a| > |b|,
  !>      = 1/c  if |b| >= |a| and c != 0
  !>      = 1    if c = 0
  !> This allows c and s to be reconstructed from z as follows:
  !>    If z = 1, set c = 0, s = 1.
  !>    If |z| < 1, set c = sqrt(1 - z**2) and s = z.
  !>    If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).
  interface la_rotg
    pure subroutine srotg(a, b, c, s)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(in) :: a
      real(wp), intent(in) :: b
      real(wp), intent(out) :: c
      real(wp), intent(out) :: s
    end subroutine srotg
    pure subroutine drotg(a, b, c, s)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(in) :: a
      real(wp), intent(in) :: b
      real(wp), intent(out) :: c
      real(wp), intent(out) :: s
    end subroutine drotg

    pure subroutine crotg(a, b, c, s)
      import :: ik, sp
      integer, parameter :: wp = sp
      complex(wp), intent(in) :: a
      complex(wp), intent(in) :: b
      real(wp), intent(out) :: c
      complex(wp), intent(out) :: s
    end subroutine crotg
    pure subroutine zrotg(a, b, c, s)
      import :: ik, dp
      integer, parameter :: wp = dp
      complex(wp), intent(in) :: a
      complex(wp), intent(in) :: b
      real(wp), intent(out) :: c
      complex(wp), intent(out) :: s
    end subroutine zrotg
  end interface la_rotg

  !> Apply the modified givens transformation, h, to the 2 by n matrix
  !>
  !> (dx**t) , where **t indicates transpose. the elements of dx are in (dy**t)
  !>
  !> dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
  !> lx = (-incx)*n, and similarly for sy using ly and incy.
  !> with dparam(1)=dflag, h has one of the following forms.
  !>
  !> dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
  !>
  !>   (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
  !> h=(          )    (          )    (          )    (          )
  !>   (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
  interface la_rotm
    pure subroutine srotm(n, x, incx, y, incy, param)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(in) :: param(5)
    end subroutine srotm
    pure subroutine drotm(n, x, incx, y, incy, param)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(in) :: param(5)
    end subroutine drotm

    module procedure :: la_rotm_rsp
    module procedure :: la_rotm_rdp
  end interface la_rotm

  !> construct the modified givens transformation matrix h which zeros
  !> the second component of the 2-vector  (dsqrt(dd1)*dx1,dsqrt(dd2)
  !> dy2)**t. With dparam(1)=dflag, h has one of the following forms..
  !>
  !> dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
  !>
  !>   (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
  !> h=(          )    (          )    (          )    (          )
  !>   (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
  !> Locations 2-4 of dparam contain dh11, dh21, dh12, and dh22
  !> respectively. (values of 1.d0, -1.d0, or 0.d0 implied by the
  !> value of dparam(1) are not stored in dparam.)
  !>
  !> The values of gamsq and rgamsq set in the data statement may be
  !> inexact.  this is ok as they are only used for testing the size
  !> of dd1 and dd2.  all actual scaling of data is done using gam.
  interface la_rotmg
    pure subroutine srotmg(d1, d2, x1, y1, param)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: d1
      real(wp), intent(inout) :: d2
      real(wp), intent(inout) :: x1
      real(wp), intent(in) :: y1
      real(wp), intent(out) :: param(5)
    end subroutine srotmg
    pure subroutine drotmg(d1, d2, x1, y1, param)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: d1
      real(wp), intent(inout) :: d2
      real(wp), intent(inout) :: x1
      real(wp), intent(in) :: y1
      real(wp), intent(out) :: param(5)
    end subroutine drotmg
  end interface la_rotmg

  !> Scales a vector by a constant.
  interface la_scal
    pure subroutine sscal(n, a, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine sscal
    pure subroutine dscal(n, a, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine dscal

    pure subroutine cscal(n, a, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine cscal
    pure subroutine zscal(n, a, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: a
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine zscal

    pure subroutine csscal(n, a, x, incx)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine csscal
    pure subroutine zdscal(n, a, x, incx)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
    end subroutine zdscal

    module procedure :: la_scal_rsp
    module procedure :: la_scal_csp
    module procedure :: la_scal_rcsp
    module procedure :: la_scal_rdp
    module procedure :: la_scal_cdp
    module procedure :: la_scal_rcdp
  end interface la_scal

  !> Applies a plane rotation.
  interface la_rot
    pure subroutine srot(n, x, incx, y, incy, c, s)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(in) :: c
      real(wp), intent(in) :: s
    end subroutine srot
    pure subroutine drot(n, x, incx, y, incy, c, s)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(in) :: c
      real(wp), intent(in) :: s
    end subroutine drot

    pure subroutine csrot(n, x, incx, y, incy, c, s)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(in) :: c
      real(wp), intent(in) :: s
    end subroutine csrot
    pure subroutine zdrot(n, x, incx, y, incy, c, s)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
      real(wp), intent(in) :: c
      real(wp), intent(in) :: s
    end subroutine zdrot

    module procedure :: la_rot_rsp
    module procedure :: la_rot_csp
    module procedure :: la_rot_rdp
    module procedure :: la_rot_cdp
  end interface la_rot

  !> Interchanges two vectors.
  interface la_swap
    pure subroutine sswap(n, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine sswap
    pure subroutine dswap(n, x, incx, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine dswap

    pure subroutine cswap(n, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine cswap
    pure subroutine zswap(n, x, incx, y, incy)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(inout) :: x(*)
      integer(ik), intent(in) :: incx
      complex(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine zswap

    module procedure :: la_swap_rsp
    module procedure :: la_swap_csp
    module procedure :: la_swap_rdp
    module procedure :: la_swap_cdp
  end interface la_swap

  !> Computes |Re(.)| + |Im(.)| of a complex number
  interface la_abs1
    pure function scabs1(z) result(abs1)
      import :: ik, sp
      integer, parameter :: wp = sp
      complex(wp), intent(in) :: z
      real(wp) :: abs1
    end function scabs1
    pure function dzabs1(z) result(abs1)
      import :: ik, dp
      integer, parameter :: wp = dp
      complex(wp), intent(in) :: z
      real(wp) :: abs1
    end function dzabs1
  end interface la_abs1

  !> Takes the sum of the absolute values.
  interface la_asum
    pure function sasum(n, x, incx) result(asum)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: asum
    end function sasum
    pure function dasum(n, x, incx) result(asum)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: asum
    end function dasum

    pure function scasum(n, x, incx) result(asum)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: asum
    end function scasum
    pure function dzasum(n, x, incx) result(asum)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: asum
    end function dzasum

    module procedure :: la_asum_rsp
    module procedure :: la_asum_csp
    module procedure :: la_asum_rdp
    module procedure :: la_asum_cdp
  end interface la_asum

  !> Returns the euclidean norm of a vector via the function
  !> name, so that
  !>
  !>    NRM2 := sqrt( x'*x )
  interface la_nrm2
    pure function snrm2(n, x, incx) result(nrm2)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: nrm2
    end function snrm2
    pure function dnrm2(n, x, incx) result(nrm2)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: nrm2
    end function dnrm2

    pure function scnrm2(n, x, incx) result(nrm2)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: nrm2
    end function scnrm2
    pure function dznrm2(n, x, incx) result(nrm2)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp) :: nrm2
    end function dznrm2

    module procedure :: la_nrm2_rsp
    module procedure :: la_nrm2_csp
    module procedure :: la_nrm2_rdp
    module procedure :: la_nrm2_cdp
  end interface la_nrm2

  !> Finds the index of the first element having maximum absolute value.
  interface la_iamax
    pure function isamax(n, x, incx) result(iamax)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: iamax
    end function isamax
    pure function idamax(n, x, incx) result(iamax)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: iamax
    end function idamax

    pure function icamax(n, x, incx) result(iamax)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: iamax
    end function icamax
    pure function izamax(n, x, incx) result(iamax)
      import :: ik, dp
      integer, parameter :: wp = dp
      integer(ik), intent(in) :: n
      complex(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      integer(ik) :: iamax
    end function izamax

    module procedure :: la_iamax_rsp
    module procedure :: la_iamax_csp
    module procedure :: la_iamax_rdp
    module procedure :: la_iamax_cdp
  end interface la_iamax

contains

  pure subroutine la_axpy_rsp(xvec, yvec, alpha)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in), optional :: alpha

    integer(ik) :: n
    real(wp) :: a

    n = size(xvec)
    a = 1.0_wp
    if (present(alpha)) a = alpha
    call la_axpy(n, a, xvec, 1, yvec, 1)
  end subroutine la_axpy_rsp

  pure subroutine la_axpy_csp(xvec, yvec, alpha)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    complex(wp), intent(in), optional :: alpha

    integer(ik) :: n
    complex(wp) :: a

    n = size(xvec)
    a = 1.0_wp
    if (present(alpha)) a = alpha
    call la_axpy(n, a, xvec, 1, yvec, 1)
  end subroutine la_axpy_csp

  pure subroutine la_axpy_rdp(xvec, yvec, alpha)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in), optional :: alpha

    integer(ik) :: n
    real(wp) :: a

    n = size(xvec)
    a = 1.0_wp
    if (present(alpha)) a = alpha
    call la_axpy(n, a, xvec, 1, yvec, 1)
  end subroutine la_axpy_rdp

  pure subroutine la_axpy_cdp(xvec, yvec, alpha)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    complex(wp), intent(in), optional :: alpha

    integer(ik) :: n
    complex(wp) :: a

    n = size(xvec)
    a = 1.0_wp
    if (present(alpha)) a = alpha
    call la_axpy(n, a, xvec, 1, yvec, 1)
  end subroutine la_axpy_cdp


  pure subroutine la_copy_rsp(xvec, yvec)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_copy(n, xvec, 1, yvec, 1)
  end subroutine la_copy_rsp

  pure subroutine la_copy_csp(xvec, yvec)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_copy(n, xvec, 1, yvec, 1)
  end subroutine la_copy_csp

  pure subroutine la_copy_rdp(xvec, yvec)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_copy(n, xvec, 1, yvec, 1)
  end subroutine la_copy_rdp

  pure subroutine la_copy_cdp(xvec, yvec)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_copy(n, xvec, 1, yvec, 1)
  end subroutine la_copy_cdp


  pure function la_dot_rsp(xvec, yvec) result(dot)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(in) :: yvec(:)
    real(wp) :: dot

    integer(ik) :: n

    n = size(xvec)
    dot = la_dot(n, xvec, 1, yvec, 1)
  end function la_dot_rsp

  pure function la_dotc_csp(xvec, yvec) result(dotc)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp) :: dotc

    integer(ik) :: n

    n = size(xvec)
    dotc = la_dotc(n, xvec, 1, yvec, 1)
  end function la_dotc_csp

  pure function la_dotu_csp(xvec, yvec) result(dotu)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp) :: dotu

    integer(ik) :: n

    n = size(xvec)
    dotu = la_dotu(n, xvec, 1, yvec, 1)
  end function la_dotu_csp

  pure function la_dot_rdp(xvec, yvec) result(dot)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp), contiguous, intent(in) :: yvec(:)
    real(wp) :: dot

    integer(ik) :: n

    n = size(xvec)
    dot = la_dot(n, xvec, 1, yvec, 1)
  end function la_dot_rdp

  pure function la_dotc_cdp(xvec, yvec) result(dotc)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp) :: dotc

    integer(ik) :: n

    n = size(xvec)
    dotc = la_dotc(n, xvec, 1, yvec, 1)
  end function la_dotc_cdp

  pure function la_dotu_cdp(xvec, yvec) result(dotu)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    complex(wp), contiguous, intent(in) :: yvec(:)
    complex(wp) :: dotu

    integer(ik) :: n

    n = size(xvec)
    dotu = la_dotu(n, xvec, 1, yvec, 1)
  end function la_dotu_cdp


  pure function la_dsdot_rspdp(xvec, yvec) result(dot)
    real(sp), contiguous, intent(in) :: xvec(:)
    real(sp), contiguous, intent(in) :: yvec(:)
    real(dp) :: dot

    integer(ik) :: n

    n = size(xvec)
    dot = la_dsdot(n, xvec, 1, yvec, 1)
  end function la_dsdot_rspdp

  pure function la_sdsdot_rsp(beta, xvec, yvec) result(dot)
    real(sp), intent(in) :: beta
    real(sp), contiguous, intent(in) :: xvec(:)
    real(sp), contiguous, intent(in) :: yvec(:)
    real(sp) :: dot

    integer(ik) :: n

    n = size(xvec)
    dot = la_sdsdot(n, beta, xvec, 1, yvec, 1)
  end function la_sdsdot_rsp

  pure subroutine la_rotm_rsp(xvec, yvec, param)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(inout) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), contiguous, intent(in) :: param(:)

    integer(ik) :: n

    n = size(xvec)
    call la_rotm(n, xvec, 1, yvec, 1, param)
  end subroutine la_rotm_rsp

  pure subroutine la_rotm_rdp(xvec, yvec, param)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(inout) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), contiguous, intent(in) :: param(:)

    integer(ik) :: n

    n = size(xvec)
    call la_rotm(n, xvec, 1, yvec, 1, param)
  end subroutine la_rotm_rdp


  pure subroutine la_scal_rsp(alpha, xvec)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: alpha
    real(wp), contiguous, intent(inout) :: xvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_scal(n, alpha, xvec, 1)
  end subroutine la_scal_rsp

  pure subroutine la_scal_rcsp(alpha, xvec)
    integer, parameter :: wp = sp
    real(wp), intent(in) :: alpha
    complex(wp), contiguous, intent(inout) :: xvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_scal(n, alpha, xvec, 1)
  end subroutine la_scal_rcsp

  pure subroutine la_scal_csp(alpha, xvec)
    integer, parameter :: wp = sp
    complex(wp), intent(in) :: alpha
    complex(wp), contiguous, intent(inout) :: xvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_scal(n, alpha, xvec, 1)
  end subroutine la_scal_csp

  pure subroutine la_scal_rdp(alpha, xvec)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: alpha
    real(wp), contiguous, intent(inout) :: xvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_scal(n, alpha, xvec, 1)
  end subroutine la_scal_rdp

  pure subroutine la_scal_rcdp(alpha, xvec)
    integer, parameter :: wp = dp
    real(wp), intent(in) :: alpha
    complex(wp), contiguous, intent(inout) :: xvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_scal(n, alpha, xvec, 1)
  end subroutine la_scal_rcdp

  pure subroutine la_scal_cdp(alpha, xvec)
    integer, parameter :: wp = dp
    complex(wp), intent(in) :: alpha
    complex(wp), contiguous, intent(inout) :: xvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_scal(n, alpha, xvec, 1)
  end subroutine la_scal_cdp


  pure subroutine la_rot_rsp(xvec, yvec, c, s)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(inout) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s

    integer(ik) :: n

    n = size(xvec)
    call la_rot(n, xvec, 1, yvec, 1, c, s)
  end subroutine la_rot_rsp

  pure subroutine la_rot_csp(xvec, yvec, c, s)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(inout) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s

    integer(ik) :: n

    n = size(xvec)
    call la_rot(n, xvec, 1, yvec, 1, c, s)
  end subroutine la_rot_csp

  pure subroutine la_rot_rdp(xvec, yvec, c, s)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(inout) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s

    integer(ik) :: n

    n = size(xvec)
    call la_rot(n, xvec, 1, yvec, 1, c, s)
  end subroutine la_rot_rdp

  pure subroutine la_rot_cdp(xvec, yvec, c, s)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(inout) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s

    integer(ik) :: n

    n = size(xvec)
    call la_rot(n, xvec, 1, yvec, 1, c, s)
  end subroutine la_rot_cdp


  pure subroutine la_swap_rsp(xvec, yvec)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(inout) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_swap(n, xvec, 1, yvec, 1)
  end subroutine la_swap_rsp

  pure subroutine la_swap_csp(xvec, yvec)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(inout) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_swap(n, xvec, 1, yvec, 1)
  end subroutine la_swap_csp

  pure subroutine la_swap_rdp(xvec, yvec)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(inout) :: xvec(:)
    real(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_swap(n, xvec, 1, yvec, 1)
  end subroutine la_swap_rdp

  pure subroutine la_swap_cdp(xvec, yvec)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(inout) :: xvec(:)
    complex(wp), contiguous, intent(inout) :: yvec(:)

    integer(ik) :: n

    n = size(xvec)
    call la_swap(n, xvec, 1, yvec, 1)
  end subroutine la_swap_cdp


  pure function la_asum_rsp(xvec) result(asum)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: asum

    integer(ik) :: n

    n = size(xvec)
    asum = la_asum(n, xvec, 1)
  end function la_asum_rsp

  pure function la_asum_csp(xvec) result(asum)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: asum

    integer(ik) :: n

    n = size(xvec)
    asum = la_asum(n, xvec, 1)
  end function la_asum_csp

  pure function la_asum_rdp(xvec) result(asum)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: asum

    integer(ik) :: n

    n = size(xvec)
    asum = la_asum(n, xvec, 1)
  end function la_asum_rdp

  pure function la_asum_cdp(xvec) result(asum)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: asum

    integer(ik) :: n

    n = size(xvec)
    asum = la_asum(n, xvec, 1)
  end function la_asum_cdp


  pure function la_nrm2_rsp(xvec) result(nrm2)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: nrm2

    integer(ik) :: n

    n = size(xvec)
    nrm2 = la_nrm2(n, xvec, 1)
  end function la_nrm2_rsp

  pure function la_nrm2_csp(xvec) result(nrm2)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: nrm2

    integer(ik) :: n

    n = size(xvec)
    nrm2 = la_nrm2(n, xvec, 1)
  end function la_nrm2_csp

  pure function la_nrm2_rdp(xvec) result(nrm2)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: nrm2

    integer(ik) :: n

    n = size(xvec)
    nrm2 = la_nrm2(n, xvec, 1)
  end function la_nrm2_rdp

  pure function la_nrm2_cdp(xvec) result(nrm2)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    real(wp) :: nrm2

    integer(ik) :: n

    n = size(xvec)
    nrm2 = la_nrm2(n, xvec, 1)
  end function la_nrm2_cdp


  pure function la_iamax_rsp(xvec) result(iamax)
    integer, parameter :: wp = sp
    real(wp), contiguous, intent(in) :: xvec(:)
    integer(ik) :: iamax

    integer(ik) :: n

    n = size(xvec)
    iamax = la_iamax(n, xvec, 1)
  end function la_iamax_rsp

  pure function la_iamax_csp(xvec) result(iamax)
    integer, parameter :: wp = sp
    complex(wp), contiguous, intent(in) :: xvec(:)
    integer(ik) :: iamax

    integer(ik) :: n

    n = size(xvec)
    iamax = la_iamax(n, xvec, 1)
  end function la_iamax_csp

  pure function la_iamax_rdp(xvec) result(iamax)
    integer, parameter :: wp = dp
    real(wp), contiguous, intent(in) :: xvec(:)
    integer(ik) :: iamax

    integer(ik) :: n

    n = size(xvec)
    iamax = la_iamax(n, xvec, 1)
  end function la_iamax_rdp

  pure function la_iamax_cdp(xvec) result(iamax)
    integer, parameter :: wp = dp
    complex(wp), contiguous, intent(in) :: xvec(:)
    integer(ik) :: iamax

    integer(ik) :: n

    n = size(xvec)
    iamax = la_iamax(n, xvec, 1)
  end function la_iamax_cdp


end module gtb_la_level1
