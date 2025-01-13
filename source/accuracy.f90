!> Kind parameters for real and integer variables
module gtb_accuracy
  implicit none
  public

  !> Single precision floating point accuracy
  integer, parameter :: sp = kind(1.0e0)

  !> Double precision floating point accuracy
  integer, parameter :: dp = kind(1.0d0)

  !> Wanted precision for real variables
  integer, parameter :: wp = dp

  !> Default integer accuracy
  integer, parameter :: ik = kind(1)

  !> 4-byte Integer
  integer, parameter :: i4 = selected_int_kind(9)

  !> Long length for integers
  integer, parameter :: i8 = selected_int_kind(18)

  !> Default logical accuracy
  integer, parameter :: lk = kind(.true.)

end module gtb_accuracy
