
! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Versioning information on this library.
module dftd4_version
   implicit none
   private

   public :: dftd4_version_string, dftd4_version_compact
   public :: get_dftd4_version


   !> String representation of the dftd4 version
   character(len=*), parameter :: dftd4_version_string = "3.2.0"

   !> Numeric representation of the dftd4 version
   integer, parameter :: dftd4_version_compact(3) = [3, 2, 0]


contains


!> Getter function to retrieve dftd4 version
subroutine get_dftd4_version(major, minor, patch, string)

   !> Major version number of the dftd4 version
   integer, intent(out), optional :: major

   !> Minor version number of the dftd4 version
   integer, intent(out), optional :: minor

   !> Patch version number of the dftd4 version
   integer, intent(out), optional :: patch

   !> String representation of the dftd4 version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = dftd4_version_compact(1)
   end if
   if (present(minor)) then
      minor = dftd4_version_compact(2)
   end if
   if (present(patch)) then
      patch = dftd4_version_compact(3)
   end if
   if (present(string)) then
      string = dftd4_version_string
   end if

end subroutine get_dftd4_version


end module dftd4_version

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_version
   implicit none
   private

   public :: mctc_version_string, mctc_version_compact
   public :: get_mctc_version


   !> String representation of the mctc-lib version
   character(len=*), parameter :: mctc_version_string = "0.2.1"

   !> Numeric representation of the mctc-lib version
   integer, parameter :: mctc_version_compact(3) = [0, 2, 1]


contains


!> Getter function to retrieve mctc-lib version
subroutine get_mctc_version(major, minor, patch, string)

   !> Major version number of the mctc-lib version
   integer, intent(out), optional :: major

   !> Minor version number of the mctc-lib version
   integer, intent(out), optional :: minor

   !> Patch version number of the mctc-lib version
   integer, intent(out), optional :: patch

   !> String representation of the mctc-lib version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = mctc_version_compact(1)
   end if
   if (present(minor)) then
      minor = mctc_version_compact(2)
   end if
   if (present(patch)) then
      patch = mctc_version_compact(3)
   end if
   if (present(string)) then
      string = mctc_version_string
   end if

end subroutine get_mctc_version


end module mctc_version

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Numerical storage size parameters for real and integer values
module mctc_env_accuracy
   implicit none
   public

   !> Single precision real numbers
   integer, parameter :: sp = selected_real_kind(6)

   !> Double precision real numbers
   integer, parameter :: dp = selected_real_kind(15)

   !> Wanted precision
   integer, parameter :: wp = dp

   !> Char length for integers
   integer, parameter :: i1 = selected_int_kind(2)

   !> Short length for integers
   integer, parameter :: i2 = selected_int_kind(4)

   !> Length of default integers
   integer, parameter :: i4 = selected_int_kind(9)

   !> Long length for integers
   integer, parameter :: i8 = selected_int_kind(18)


end module mctc_env_accuracy

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Central registry for error codes
module mctc_env_error
   implicit none
   private

   public :: mctc_stat, error_type
   public :: fatal_error


   !> Possible error codes
   type :: enum_stat

      !> Successful run
      integer :: success = 0

      !> Internal error:
      integer :: fatal = 1

   end type enum_stat

   !> Actual enumerator for return states
   type(enum_stat), parameter :: mctc_stat = enum_stat()


   !> Error message
   type :: error_type

      !> Error code
      integer :: stat

      !> Payload of the error
      character(len=:), allocatable :: message

   end type error_type


contains


!> A fatal error is encountered
subroutine fatal_error(error, message, stat)

   !> Instance of the error
   type(error_type), allocatable, intent(out) :: error

   !> A detailed message describing the error and (optionally) offering advice
   character(len=*), intent(in), optional :: message

   !> Overwrite of the error code
   integer, intent(in), optional :: stat

   allocate(error)

   if (present(stat)) then
      error%stat = stat
   else
      error%stat = mctc_stat%fatal
   end if

   if (present(message)) then
      error%message = message
   else
      error%message = "Fatal error"
   end if

end subroutine fatal_error


end module mctc_env_error

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Module collecting commands to conveniently interface with system commands
module mctc_env_system
   implicit none
   private

   public :: get_argument, get_variable
   public :: is_windows, is_unix


contains


!> Obtain the command line argument at a given index
subroutine get_argument(idx, arg)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: idx

   !> Command line argument
   character(len=:), allocatable, intent(out) :: arg

   integer :: length, stat

   call get_command_argument(idx, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: arg, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_command_argument(idx, arg, status=stat)
      if (stat /= 0) then
         deallocate(arg)
         return
      end if
   end if

end subroutine get_argument


!> Obtain the value of an environment variable
subroutine get_variable(var, val)

   !> Name of variable
   character(len=*), intent(in) :: var

   !> Value of variable
   character(len=:), allocatable, intent(out) :: val

   integer :: length, stat

   call get_environment_variable(var, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: val, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_environment_variable(var, val, status=stat)
      if (stat /= 0) then
         deallocate(val)
         return
      end if
   end if

end subroutine get_variable


!> Try to determine if we run on Windows and don't have POSIX compliance around
function is_windows()

   !> Operating system seems to be Windows
   logical :: is_windows

   character(len=:), allocatable :: tmp

   is_windows = .false.
   call get_variable('OS', tmp)
   if (allocated(tmp)) then
      is_windows = index(tmp, 'Windows_NT') > 0
   end if
   if (.not.is_windows) then
      call get_variable('OSTYPE', tmp)
      if (allocated(tmp)) then
         is_windows = index(tmp, 'win') > 0 .or. index(tmp, 'msys') > 0
      end if
   end if

end function is_windows


!> Try to determine if we run on Unix and probably can rely on POSIX compliance
function is_unix()

   !> Operating system seems to be Unix
   logical :: is_unix

   character(len=:), allocatable :: tmp

   is_unix = .not. is_windows()

end function is_unix


end module mctc_env_system

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> File type support
module mctc_io_filetype
   implicit none
   private

   public :: filetype, get_filetype


   !> Possible file types
   type :: enum_filetype

      !> Unknown file type
      integer :: unknown = 0

      !> xyz-format
      integer :: xyz = 1

      !> Turbomole coordinate format
      integer :: tmol = 2

      !> mol-format
      integer :: molfile = 3

      !> Vasp coordinate input
      integer :: vasp = 4

      !> Protein database format
      integer :: pdb = 5

      !> Structure data format
      integer :: sdf = 6

      !> GenFormat of DFTB+
      integer :: gen = 7

      !> Gaussian external format
      integer :: gaussian = 8

   end type enum_filetype

   !> File type enumerator
   type(enum_filetype), parameter :: filetype = enum_filetype()


contains


elemental function get_filetype(file) result(ftype)

   !> Name of the file
   character(len=*), intent(in) :: file

   !> File type from extension
   integer :: ftype

   integer :: iext, isep

   ftype = filetype%unknown
   iext = index(file, '.', back=.true.)
   isep = scan(file, '\/', back=.true.)

   if (iext > isep .and. iext > 0) then
      select case(to_lower(file(iext+1:)))
      case('coord', 'tmol')
         ftype = filetype%tmol
      case('xyz', 'log')
         ftype = filetype%xyz
      case('mol')
         ftype = filetype%molfile
      case('sdf')
         ftype = filetype%sdf
      case('poscar', 'contcar', 'vasp')
         ftype = filetype%vasp
      case('pdb')
         ftype = filetype%pdb
      case('gen')
         ftype = filetype%gen
      case('ein')
         ftype = filetype%gaussian
      end select
      if (ftype /= filetype%unknown) return
   else
      iext = len(file) + 1
   end if

   if (iext > isep) then
      select case(to_lower(file(isep+1:iext-1)))
      case('coord')
         ftype = filetype%tmol
      case('poscar', 'contcar')
         ftype = filetype%vasp
      end select
   end if

end function get_filetype


!> Convert input string to lowercase
elemental function to_lower(str) result(lcstr)

   !> Input string
   character(len=*), intent(in) :: str

   !> Lowercase version of string
   character(len=len(str)):: lcstr

   integer :: ilen, iquote, i, iav, iqc
   integer, parameter :: offset = iachar('A') - iachar('a')

   ilen = len(str)
   iquote = 0
   lcstr = str

   do i = 1, ilen
      iav = iachar(str(i:i))
      if (iquote == 0 .and. (iav == 34 .or.iav == 39)) then
         iquote = 1
         iqc = iav
         cycle
      end if
      if (iquote == 1 .and. iav==iqc) then
         iquote=0
         cycle
      end if
      if (iquote == 1) cycle
      if (iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr(i:i) = achar(iav - offset)
      else
         lcstr(i:i) = str(i:i)
      end if
   end do

end function to_lower


end module mctc_io_filetype

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_utils
   implicit none
   private

   public :: getline


contains


subroutine getline(unit, line, iostat, iomsg)

   !> Formatted IO unit
   integer, intent(in) :: unit

   !> Line to read
   character(len=:), allocatable, intent(out) :: line

   !> Status of operation
   integer, intent(out) :: iostat

   !> Error message
   character(len=:), allocatable, optional :: iomsg

   integer, parameter :: bufsize = 512
   character(len=bufsize) :: buffer
   character(len=bufsize) :: msg
   integer :: size
   integer :: stat

   allocate(character(len=0) :: line)
   do
      read(unit, '(a)', advance='no', iostat=stat, iomsg=msg, size=size) &
         & buffer
      if (stat > 0) exit
      line = line // buffer(:size)
      if (stat < 0) then
         if (is_iostat_eor(stat)) then
            stat = 0
         end if
         exit
      end if
   end do

   if (stat /= 0) then
      if (present(iomsg)) iomsg = trim(msg)
   end if
   iostat = stat

end subroutine getline


end module mctc_io_utils

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Versioning information on this library.
module multicharge_version
   implicit none
   private

   public :: multicharge_version_string, multicharge_version_compact
   public :: get_multicharge_version


   !> String representation of the multicharge version
   character(len=*), parameter :: multicharge_version_string = "0.1.0"

   !> Numeric representation of the multicharge version
   integer, parameter :: multicharge_version_compact(3) = [0, 1, 0]


contains


!> Getter function to retrieve multicharge version
subroutine get_multicharge_version(major, minor, patch, string)

   !> Major version number of the multicharge version
   integer, intent(out), optional :: major

   !> Minor version number of the multicharge version
   integer, intent(out), optional :: minor

   !> Patch version number of the multicharge version
   integer, intent(out), optional :: patch

   !> String representation of the multicharge version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = multicharge_version_compact(1)
   end if
   if (present(minor)) then
      minor = multicharge_version_compact(2)
   end if
   if (present(patch)) then
      patch = multicharge_version_compact(3)
   end if
   if (present(string)) then
      string = multicharge_version_string
   end if

end subroutine get_multicharge_version


end module multicharge_version

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Public API reexport of environment library
module mctc_env
   use mctc_env_accuracy, only : sp, dp, wp, i1, i2, i4, i8
   use mctc_env_error, only : error_type, fatal_error, mctc_stat
   use mctc_env_system, only : get_argument, get_variable, &
      & is_unix, is_windows
   implicit none
   public

end module mctc_env

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Provides a light-weight testing framework for usage in projects depending on
!> the tool chain library.
!>
!> Testsuites are defined by a [[collect_interface]] returning a set of
!> [[unittest_type]] objects. To create a new test use the [[new_unittest]]
!> constructor, which requires a test identifier and a procedure with a
!> [[test_interface]] compatible signature. The error status is communicated
!> by the allocation status of an [[error_type]].
!>
!> The necessary boilerplate code to setup the test entry point is just
!>
!>```fortran
!>program tester
!>   use, intrinsic :: iso_fortran_env, only : error_unit
!>   use mctc_env_testing, only : run_testsuite, new_testsuite, testsuite_type
!>   use test_suite1, only : collect_suite1
!>   use test_suite2, only : collect_suite2
!>   implicit none
!>   integer :: stat, ii
!>   type(testsuite_type), allocatable :: testsuites(:)
!>   character(len=*), parameter :: fmt = '("#", *(1x, a))'
!>
!>   stat = 0
!>
!>   testsuites = [ &
!>      & new_testsuite("suite1", collect_suite1), &
!>      & new_testsuite("suite2", collect_suite2) &
!>      & ]
!>
!>   do ii = 1, size(testsuites)
!>      write(error_unit, fmt) "Testing:", testsuites(ii)%name
!>      call run_testsuite(testsuites(ii)%collect, error_unit, stat)
!>   end do
!>
!>   if (stat > 0) then
!>      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
!>      error stop
!>   end if
!>
!>end program tester
!>```
!>
!> Every test is defined in a separate module using a ``collect`` function, which
!> is exported and added to the ``testsuites`` array in the test runner.
!> All test have a simple interface with just an allocatable [[error_type]] as
!> output to provide the test results.
!>
!>```fortran
!>module test_suite1
!>   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
!>   implicit none
!>   private
!>
!>   public :: collect_suite1
!>
!>contains
!>
!>!> Collect all exported unit tests
!>subroutine collect_suite1(testsuite)
!>   !> Collection of tests
!>   type(unittest_type), allocatable, intent(out) :: testsuite(:)
!>
!>   testsuite = [ &
!>      & new_unittest("valid", test_valid), &
!>      & new_unittest("invalid", test_invalid, should_fail=.true.) &
!>      & ]
!>
!>end subroutine collect_suite1
!>
!>subroutine test_valid(error)
!>   type(error_type), allocatable, intent(out) :: error
!>   ! ...
!>end subroutine test_valid
!>
!>subroutine test_invalid(error)
!>   type(error_type), allocatable, intent(out) :: error
!>   ! ...
!>end subroutine test_invalid
!>
!>end module test_suite1
!>```
!>
!> For an example setup checkout the ``test/`` directory in this project.
module mctc_env_testing
   use mctc_env_error, only : error_type, mctc_stat
   use mctc_env_accuracy, only : sp, dp, i1, i2, i4, i8
   implicit none
   private

   public :: run_testsuite, run_selected, new_unittest, new_testsuite
   public :: select_test, select_suite
   public :: unittest_type, testsuite_type, error_type
   public :: check, test_failed
   public :: test_interface, collect_interface


   interface check
      module procedure :: check_stat
      module procedure :: check_logical
      module procedure :: check_float_sp
      module procedure :: check_float_dp
      module procedure :: check_int_i1
      module procedure :: check_int_i2
      module procedure :: check_int_i4
      module procedure :: check_int_i8
      module procedure :: check_bool
      module procedure :: check_string
   end interface check


   abstract interface
      !> Entry point for tests
      subroutine test_interface(error)
         import :: error_type

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

      end subroutine test_interface
   end interface


   !> Declaration of a unit test
   type :: unittest_type

      !> Name of the test
      character(len=:), allocatable :: name

      !> Entry point of the test
      procedure(test_interface), pointer, nopass :: test => null()

      !> Whether test is supposed to fail
      logical :: should_fail = .false.

   end type unittest_type


   abstract interface
      !> Collect all tests
      subroutine collect_interface(testsuite)
         import :: unittest_type

         !> Collection of tests
         type(unittest_type), allocatable, intent(out) :: testsuite(:)

      end subroutine collect_interface
   end interface


   !> Collection of unit tests
   type :: testsuite_type

      !> Name of the testsuite
      character(len=:), allocatable :: name

      !> Entry point of the test
      procedure(collect_interface), pointer, nopass :: collect => null()

   end type testsuite_type


   character(len=*), parameter :: fmt = '(1x, *(1x, a))'
   character(len=*), parameter :: indent = repeat(" ", 5) // repeat(".", 3)


contains


!> Driver for testsuite
subroutine run_testsuite(collect, unit, stat)

   !> Collect tests
   procedure(collect_interface) :: collect

   !> Unit for IO
   integer, intent(in) :: unit

   !> Number of failed tests
   integer, intent(inout) :: stat

   type(unittest_type), allocatable :: testsuite(:)
   integer :: ii

   call collect(testsuite)

   !$omp parallel do shared(testsuite, unit) reduction(+:stat)
   do ii = 1, size(testsuite)
      !$omp critical(mctc_env_testsuite)
      write(unit, '(1x, 3(1x, a), 1x, "(", i0, "/", i0, ")")') &
         & "Starting", testsuite(ii)%name, "...", ii, size(testsuite)
      !$omp end critical(mctc_env_testsuite)
      call run_unittest(testsuite(ii), unit, stat)
   end do

end subroutine run_testsuite


!> Driver for selective testing
subroutine run_selected(collect, name, unit, stat)

   !> Collect tests
   procedure(collect_interface) :: collect

   !> Name of the selected test
   character(len=*), intent(in) :: name

   !> Unit for IO
   integer, intent(in) :: unit

   !> Number of failed tests
   integer, intent(inout) :: stat

   type(unittest_type), allocatable :: testsuite(:)
   integer :: ii

   call collect(testsuite)

   ii = select_test(testsuite, name)

   if (ii > 0 .and. ii <= size(testsuite)) then
      call run_unittest(testsuite(ii), unit, stat)
   else
      write(unit, fmt) "Available tests:"
      do ii = 1, size(testsuite)
         write(unit, fmt) "-", testsuite(ii)%name
      end do
      stat = -huge(ii)
   end if

end subroutine run_selected


!> Run a selected unit test
subroutine run_unittest(test, unit, stat)

   !> Unit test
   type(unittest_type), intent(in) :: test

   !> Unit for IO
   integer, intent(in) :: unit

   !> Number of failed tests
   integer, intent(inout) :: stat

   type(error_type), allocatable :: error

   call test%test(error)
   !$omp critical(mctc_env_testsuite)
   if (allocated(error) .neqv. test%should_fail) then
      if (test%should_fail) then
         write(unit, fmt) indent, test%name, "[UNEXPECTED PASS]"
      else
         write(unit, fmt) indent, test%name, "[FAILED]"
      end if
      stat = stat + 1
   else
      if (test%should_fail) then
         write(unit, fmt) indent, test%name, "[EXPECTED FAIL]"
      else
         write(unit, fmt) indent, test%name, "[PASSED]"
      end if
   end if
   if (allocated(error)) then
      write(unit, fmt) "Message:", error%message
   end if
   !$omp end critical(mctc_env_testsuite)

end subroutine run_unittest


!> Select a unit test from all available tests
function select_test(tests, name) result(pos)

   !> Name identifying the test suite
   character(len=*), intent(in) :: name

   !> Available unit tests
   type(unittest_type) :: tests(:)

   !> Selected test suite
   integer :: pos

   integer :: it

   pos = 0
   do it = 1, size(tests)
      if (name == tests(it)%name) then
         pos = it
         exit
      end if
   end do

end function select_test


!> Select a test suite from all available suites
function select_suite(suites, name) result(pos)

   !> Name identifying the test suite
   character(len=*), intent(in) :: name

   !> Available test suites
   type(testsuite_type) :: suites(:)

   !> Selected test suite
   integer :: pos

   integer :: it

   pos = 0
   do it = 1, size(suites)
      if (name == suites(it)%name) then
         pos = it
         exit
      end if
   end do

end function select_suite


!> Register a new unit test
function new_unittest(name, test, should_fail) result(self)

   !> Name of the test
   character(len=*), intent(in) :: name

   !> Entry point for the test
   procedure(test_interface) :: test

   !> Whether test is supposed to error or not
   logical, intent(in), optional :: should_fail

   !> Newly registered test
   type(unittest_type) :: self

   self%name = name
   self%test => test
   if (present(should_fail)) self%should_fail = should_fail

end function new_unittest


!> Register a new testsuite
function new_testsuite(name, collect) result(self)

   !> Name of the testsuite
   character(len=*), intent(in) :: name

   !> Entry point to collect tests
   procedure(collect_interface) :: collect

   !> Newly registered testsuite
   type(testsuite_type) :: self

   self%name = name
   self%collect => collect

end function new_testsuite


subroutine check_stat(error, stat, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Status of operation
   integer, intent(in) :: stat

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (stat /= mctc_stat%success) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Non-zero exit code encountered", more)
      end if
   end if

end subroutine check_stat


subroutine check_logical(error, expression, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Result of logical operator
   logical, intent(in) :: expression

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (.not.expression) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Condition not fullfilled", more)
      end if
   end if

end subroutine check_logical


subroutine check_float_dp(error, actual, expected, message, more, thr, rel)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found floating point value
   real(dp), intent(in) :: actual

   !> Expected floating point value
   real(dp), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   !> Allowed threshold for matching floating point values
   real(dp), intent(in), optional :: thr

   !> Check for relative errors instead
   logical, intent(in), optional :: rel

   logical :: relative
   real(dp) :: diff, threshold

   if (present(thr)) then
      threshold = thr
   else
      threshold = epsilon(expected)
   end if

   if (present(rel)) then
      relative = rel
   else
      relative = .false.
   end if

   if (relative) then
      diff = abs(actual - expected) / expected
   else
      diff = abs(actual - expected)
   end if

   if (diff > threshold) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Floating point value missmatch", more)
      end if
   end if

end subroutine check_float_dp


subroutine check_float_sp(error, actual, expected, message, more, thr, rel)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found floating point value
   real(sp), intent(in) :: actual

   !> Expected floating point value
   real(sp), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   !> Allowed threshold for matching floating point values
   real(sp), intent(in), optional :: thr

   !> Check for relative errors instead
   logical, intent(in), optional :: rel

   logical :: relative
   real(sp) :: diff, threshold

   if (present(thr)) then
      threshold = thr
   else
      threshold = epsilon(expected)
   end if

   if (present(rel)) then
      relative = rel
   else
      relative = .false.
   end if

   if (relative) then
      diff = abs(actual - expected) / expected
   else
      diff = abs(actual - expected)
   end if

   if (diff > threshold) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Floating point value missmatch", more)
      end if
   end if

end subroutine check_float_sp


subroutine check_int_i1(error, actual, expected, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found integer value
   integer(i1), intent(in) :: actual

   !> Expected integer value
   integer(i1), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (expected /= actual) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Integer value missmatch", more)
      end if
   end if

end subroutine check_int_i1


subroutine check_int_i2(error, actual, expected, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found integer value
   integer(i2), intent(in) :: actual

   !> Expected integer value
   integer(i2), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (expected /= actual) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Integer value missmatch", more)
      end if
   end if

end subroutine check_int_i2


subroutine check_int_i4(error, actual, expected, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found integer value
   integer(i4), intent(in) :: actual

   !> Expected integer value
   integer(i4), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (expected /= actual) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Integer value missmatch", more)
      end if
   end if

end subroutine check_int_i4


subroutine check_int_i8(error, actual, expected, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found integer value
   integer(i8), intent(in) :: actual

   !> Expected integer value
   integer(i8), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (expected /= actual) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Integer value missmatch", more)
      end if
   end if

end subroutine check_int_i8


subroutine check_bool(error, actual, expected, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found boolean value
   logical, intent(in) :: actual

   !> Expected boolean value
   logical, intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (expected .neqv. actual) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Logical value missmatch", more)
      end if
   end if

end subroutine check_bool


subroutine check_string(error, actual, expected, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Found boolean value
   character(len=*), intent(in) :: actual

   !> Expected boolean value
   character(len=*), intent(in) :: expected

   !> A detailed message describing the error
   character(len=*), intent(in), optional :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   if (expected /= actual) then
      if (present(message)) then
         call test_failed(error, message, more)
      else
         call test_failed(error, "Character value missmatch", more)
      end if
   end if

end subroutine check_string


subroutine test_failed(error, message, more)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> A detailed message describing the error
   character(len=*), intent(in) :: message

   !> Another line of error message
   character(len=*), intent(in), optional :: more

   allocate(error)
   error%stat = mctc_stat%fatal

   if (present(more)) then
      error%message = message // new_line('a') // more
   else
      error%message = message
   end if

end subroutine test_failed


end module mctc_env_testing

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Numerical constants
module mctc_io_constants
   use mctc_env_accuracy, only : wp
   implicit none
   private

   public :: pi, codata


   !> Ratio between a circles diameter and its circumfence
   real(wp), parameter :: pi = 3.1415926535897932384626433832795029_wp


   !> Natural constants defining the SI unit base
   type :: enum_codata

      !> Planck's constant
      real(wp) :: h = 6.6260715e-34_wp ! J·s = kg·m²·s⁻¹

      !> Speed of light in vacuum
      real(wp) :: c = 299792458.0_wp ! m·s⁻¹

      !> Boltzmann's constant
      real(wp) :: kb = 1.380649e-23_wp ! J·K⁻¹ = kg·m²·s⁻²·K⁻¹

      !> Avogadro's number
      real(wp) :: NA = 6.02214076e23_wp ! mol⁻¹

      !> Elementary charge
      real(wp) :: e = 1.602176634e-19_wp ! C

      !> fine structure constant (CODATA2018)
      real(wp) :: alpha = 1.0_wp/137.035999046_wp ! dimensionless

      !> electron rest mass
      real(wp) :: me = 9.10938356e-31_wp ! kg

   end type enum_codata

   !> Actual collection of natural constants
   type(enum_codata), parameter :: codata = enum_codata()


end module mctc_io_constants

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Reallocation implementation for resizing arrays
module mctc_io_resize
   use mctc_env_accuracy, only : wp
   implicit none
   private

   public :: resize


   !> Overloaded resize interface
   interface resize
      module procedure :: resize_char
      module procedure :: resize_int
      module procedure :: resize_real
      module procedure :: resize_real_2d
   end interface resize


   !> Initial size for dynamic sized arrays
   integer, parameter :: initial_size = 64


contains


!> Reallocate list of integers
pure subroutine resize_int(var, n)

   !> Instance of the array to be resized
   integer, allocatable, intent(inout) :: var(:)

   !> Dimension of the final array size
   integer, intent(in), optional :: n

   integer, allocatable :: tmp(:)
   integer :: this_size, new_size

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if

end subroutine resize_int


!> Reallocate list of characters
pure subroutine resize_char(var, n)

   !> Instance of the array to be resized
   character(len=*), allocatable, intent(inout) :: var(:)

   !> Dimension of the final array size
   integer, intent(in), optional :: n

   character(len=:), allocatable :: tmp(:)
   integer :: this_size, new_size

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if

end subroutine resize_char


!> Reallocate list of reals
pure subroutine resize_real(var, n)

   !> Instance of the array to be resized
   real(wp), allocatable, intent(inout) :: var(:)

   !> Dimension of the final array size
   integer, intent(in), optional :: n

   real(wp), allocatable :: tmp(:)
   integer :: this_size, new_size

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if

end subroutine resize_real


!> Reallocate list of reals
pure subroutine resize_real_2d(var, n)

   !> Instance of the array to be resized
   real(wp), allocatable, intent(inout) :: var(:,:)

   !> Dimension of the final array size
   integer, intent(in), optional :: n

   real(wp), allocatable :: tmp(:,:)
   integer :: order, this_size, new_size

   if (allocated(var)) then
      order = size(var, 1)
      this_size = size(var, 2)
      call move_alloc(var, tmp)
   else
      order = 3
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(order, new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 2), size(var, 2))
      var(:, :this_size) = tmp(:, :this_size)
      deallocate(tmp)
   end if

end subroutine resize_real_2d


end module mctc_io_resize

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_structure_info
   use mctc_env_accuracy, only : wp
   implicit none
   private

   public :: pdb_data, sdf_data, structure_info
   public :: resize


   !> Atomic pdb data type.
   !>
   !> keeps information from PDB input that is currently not used by the
   !> caller program (like residues or chains) but is needed to write
   !> the PDB output eventually
   !>
   !>     ATOM   2461  HA3 GLY A 153     -10.977  -7.661   2.011  1.00  0.00           H
   !>     TER    2462      GLY A 153
   !>     a6----i5---xa4--aa3-xai4--axxxf8.3----f8.3----f8.3----f6.2--f6.2--xxxxxxa4--a2a2
   !>     HETATM 2463  CHA HEM A 154       9.596 -13.100  10.368  1.00  0.00           C
   type :: pdb_data
      logical :: het = .false.
      integer :: charge = 0
      integer :: residue_number = 0
      character(len=4) :: name = ' '
      character(len=1) :: loc = ' '
      character(len=3) :: residue = ' '
      character(len=1) :: chains = ' '
      character(len=1) :: code = ' '
      character(len=4) :: segid = ' '
   end type pdb_data


   !> SDF atomic data.
   !>
   !> We only support some entries, the rest is simply dropped.
   !> the format is: ddcccssshhhbbbvvvHHHrrriiimmmnnneee
   type :: sdf_data
      integer :: isotope = 0   !< d field
      integer :: charge = 0    !< c field
      integer :: hydrogens = 0 !< h field
      integer :: valence = 0   !< v field
   end type sdf_data


   !> structure input info
   !>
   !> contains informations from different input file formats
   type :: structure_info

      !> Vasp coordinate scaling information
      real(wp) :: scale = 1.0_wp

      !> Vasp selective dynamics keyword is present
      logical :: selective = .false.

      !> SDF 2D structure present
      logical :: two_dimensional = .false.

      !> SDF hydrogen query present or PDB without hydrogen atoms found
      logical :: missing_hydrogen = .false.

      !> Periodic coordinates should use preferrably cartesian coordinates
      logical :: cartesian = .true.

      !> Lattice information should use preferrably lattice vectors
      logical :: lattice = .true.

      !> Unit of the lattice vectors should be in Angstrom if possible
      logical :: angs_lattice = .false.

      !> Unit of the atomic coordinates should be in Angstrom if possible
      logical :: angs_coord = .false.

   end type structure_info


   interface resize
      module procedure resize_pdb_data
   end interface


contains


subroutine resize_pdb_data(var, n)
   type(pdb_data), allocatable, intent(inout) :: var(:)
   integer, intent(in), optional :: n
   type(pdb_data), allocatable :: tmp(:)
   integer :: length, current_length
   current_length = size(var)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(tmp(length), source=pdb_data())
      tmp(:current_length) = var(:current_length)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(var(length), source=pdb_data())
   endif
end subroutine resize_pdb_data


end module mctc_io_structure_info

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Interface to BLAS library
module dftd4_blas
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: d4_gemv, blas_gemv


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface d4_gemv
      module procedure :: d4_sgemv
      module procedure :: d4_dgemv
      module procedure :: d4_sgemv312
      module procedure :: d4_sgemv321
      module procedure :: d4_dgemv312
      module procedure :: d4_dgemv321
   end interface d4_gemv


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface blas_gemv
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgemv
   end interface blas_gemv


contains


subroutine d4_sgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call d4_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine d4_sgemv312


subroutine d4_sgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call d4_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine d4_sgemv321


subroutine d4_dgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call d4_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine d4_dgemv312


subroutine d4_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call d4_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine d4_dgemv321


pure subroutine d4_sgemv(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine d4_sgemv


pure subroutine d4_dgemv(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine d4_dgemv


end module dftd4_blas

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Realspace cutoff and lattice point generator utilities
module dftd4_cutoff
   use mctc_env, only : wp
   implicit none
   private

   public :: realspace_cutoff, get_lattice_points


   !> Coordination number cutoff
   real(wp), parameter :: cn_default = 30.0_wp

   !> Two-body interaction cutoff
   real(wp), parameter :: disp2_default = 60.0_wp

   !> Three-body interaction cutoff
   real(wp), parameter :: disp3_default = 40.0_wp


   !> Collection of real space cutoffs
   type :: realspace_cutoff
      sequence

      !> Coordination number cutoff
      real(wp) :: cn = cn_default

      !> Two-body interaction cutoff
      real(wp) :: disp2 = disp2_default

      !> Three-body interaction cutoff
      real(wp) :: disp3 = disp3_default

   end type realspace_cutoff


   interface get_lattice_points
      module procedure :: get_lattice_points_cutoff
      module procedure :: get_lattice_points_rep_3d
   end interface get_lattice_points


contains


!> Generate lattice points from repeatitions
subroutine get_lattice_points_rep_3d(lat, rep, origin, trans)

   !> Lattice vectors
   real(wp), intent(in) :: lat(:, :)

   !> Repeatitions of lattice points to generate
   integer, intent(in) :: rep(:)

   !> Include the origin in the generated lattice points
   logical, intent(in) :: origin

   !> Generated lattice points
   real(wp), allocatable, intent(out) :: trans(:, :)

   integer :: itr, ix, iy, iz, jx, jy, jz

   itr = 0
   if (origin) then
      allocate(trans(3, product(2*rep+1)))
      do ix = 0, rep(1)
         do iy = 0, rep(2)
            do iz = 0, rep(3)
               do jx = 1, merge(-1, 1, ix > 0), -2
                  do jy = 1, merge(-1, 1, iy > 0), -2
                     do jz = 1, merge(-1, 1, iz > 0), -2
                        itr = itr + 1
                        trans(:, itr) = lat(:, 1)*ix*jx &
                           & + lat(:, 2)*iy*jy + lat(:, 3)*iz*jz
                     end do
                  end do
               end do
            end do
         end do
      end do
   else
      allocate(trans(3, product(2*rep+1)-1))
      do ix = 0, rep(1)
         do iy = 0, rep(2)
            do iz = 0, rep(3)
               if (ix == 0 .and. iy == 0 .and. iz == 0) cycle
               do jx = 1, merge(-1, 1, ix > 0), -2
                  do jy = 1, merge(-1, 1, iy > 0), -2
                     do jz = 1, merge(-1, 1, iz > 0), -2
                        itr = itr + 1
                        trans(:, itr) = lat(:, 1)*ix*jx &
                           & + lat(:, 2)*iy*jy + lat(:, 3)*iz*jz
                     end do
                  end do
               end do
            end do
         end do
      end do
   end if

end subroutine get_lattice_points_rep_3d


!> Create lattice points within a given cutoff
subroutine get_lattice_points_cutoff(periodic, lat, rthr, trans)

   !> Periodic dimensions
   logical, intent(in) :: periodic(:)

   !> Real space cutoff
   real(wp), intent(in) :: rthr

   !> Lattice parameters
   real(wp), intent(in) :: lat(:, :)

   !> Generated lattice points
   real(wp), allocatable, intent(out) :: trans(:, :)

   integer :: rep(3)

   if (.not.any(periodic)) then
      allocate(trans(3, 1))
      trans(:, :) = 0.0_wp
   else
      call get_translations(lat, rthr, rep)
      call get_lattice_points(lat, rep, .true., trans)
   end if

end subroutine get_lattice_points_cutoff


!> Generate a supercell based on a realspace cutoff, this subroutine
!> doesn't know anything about the convergence behaviour of the
!> associated property.
pure subroutine get_translations(lat, rthr, rep)
   real(wp), intent(in) :: rthr
   real(wp), intent(in) :: lat(3, 3)
   integer, intent(out) :: rep(3)
   real(wp) :: normx(3), normy(3), normz(3)
   real(wp) :: cos10, cos21, cos32

   ! find normal to the plane...
   call crossproduct(lat(:, 2), lat(:, 3), normx)
   call crossproduct(lat(:, 3), lat(:, 1), normy)
   call crossproduct(lat(:, 1), lat(:, 2), normz)
   ! ...normalize it...
   normx = normx/norm2(normx)
   normy = normy/norm2(normy)
   normz = normz/norm2(normz)
   ! cos angles between normals and lattice vectors
   cos10 = sum(normx*lat(:, 1))
   cos21 = sum(normy*lat(:, 2))
   cos32 = sum(normz*lat(:, 3))
   rep(1) = ceiling(abs(rthr/cos10))
   rep(2) = ceiling(abs(rthr/cos21))
   rep(3) = ceiling(abs(rthr/cos32))

contains

   pure subroutine crossproduct(a, b, c)
      real(wp), intent(in) :: a(3)
      real(wp), intent(in) :: b(3)
      real(wp), intent(out) :: c(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
   end subroutine crossproduct

end subroutine get_translations


end module dftd4_cutoff

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Conversion factors
module mctc_io_convert
   use mctc_env_accuracy, only : wp
   use mctc_io_constants, only : pi, codata
   implicit none
   private

   !> Reducted Planck's constant
   real(wp), parameter :: hbar = codata%h/(2.0_wp*pi) ! J·s = kg·m²·s⁻¹

   !> Bohr radius
   real(wp), parameter :: bohr = hbar/(codata%me*codata%c*codata%alpha) ! m

   !> Hartree energy
   real(wp), parameter :: hartree = codata%me*codata%c**2*codata%alpha**2 ! J = kg·m²·s⁻²

   !> Conversion factor from bohr to Ångström
   real(wp), public, parameter :: autoaa = bohr * 1e10_wp

   !> Conversion factor from Ångström to bohr
   real(wp), public, parameter :: aatoau = 1.0_wp/autoaa

   !> Conversion factor from hartree to electron volts
   real(wp), public, parameter :: autoeV = hartree/codata%e

   !> Conversion factor from electron volts to hartree
   real(wp), public, parameter :: evtoau = 1.0_wp/autoev

   !> Coversion factor between calorine and joule
   real(wp), public, parameter :: caltoj = 4.184_wp

   !> Coversion factor between joule and calorine
   real(wp), public, parameter :: jtocal = 1.0_wp/caltoj

   !> Conversion from hartree to kJ/mol
   real(wp), public, parameter :: autokj = hartree*codata%na*1e-3_wp

   !> Conversion from kJ/mol to hartree
   real(wp), public, parameter :: kjtoau = 1.0_wp/autokj

   !> Conversion from hartree to kcal/mol
   real(wp), public, parameter :: autokcal = autokJ*Jtocal

   !> Conversion from kcal/mol to hartree
   real(wp), public, parameter :: kcaltoau = 1.0_wp/autokcal

   !> Conversion from hartree to reciprocal centimeters
   real(wp), public, parameter :: autorcm = hartree/(codata%h*codata%c)*1e-2_wp

   !> Conversion from reciprocal centimeters to hartree
   real(wp), public, parameter :: rcmtoau = 1.0_wp/autorcm

   !> Conversion from hartree to nanometers (wavelength)
   real(wp), public, parameter :: autonm = codata%h*codata%c/hartree * 1e+9_wp

   !> Conversion from nanometers (wavelength) to hartree
   real(wp), public, parameter :: nmtoau = 1.0_wp/autonm

   !> Conversion from electron mass (a.u.) to kg
   real(wp), public, parameter :: autokg = codata%me

   !> Conversion from kg to electron mass (a.u.)
   real(wp), public, parameter :: kgtoau = 1.0_wp/autokg

   !> Molecular mass per mole (g/mol) to electron mass (a.u.)
   real(wp), public, parameter :: autogmol = codata%me*codata%na*1e+3_wp

   !> Electron mass (a.u.) to molecular mass per mole (g/mol)
   real(wp), public, parameter :: gmoltoau = 1.0_wp/autogmol

   !> Molecular mass per mole (g/mol) to kg
   real(wp), public, parameter :: gmoltokg = gmoltoau*autokg

   !> kg to molecular mass per mole (g/mol)
   real(wp), public, parameter :: kgtogmol = 1.0_wp/gmoltokg

   !> Coulomb to atomic charge units
   real(wp), public, parameter :: autoc = codata%e

   !> Atomic charge units to Coulomb
   real(wp), public, parameter :: ctoau = 1.0_wp/autoc


end module mctc_io_convert

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Simple algebraic functions
module mctc_io_math
   use mctc_env_accuracy, only : wp
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: matdet_3x3, matinv_3x3, crossprod, eigval_3x3, eigvec_3x3


   real(wp), parameter :: twothirdpi = 2.0_wp * pi / 3.0_wp


contains


!> Performs a direct calculation of the inverse of a 3×3 matrix.
!
!  reference: http://fortranwiki.org/fortran/show/Matrix+inversion
pure function matinv_3x3(a) result(b)

   !> Matrix
   real(wp), intent(in) :: a(3, 3)

   !> Inverse matrix
   real(wp) :: b(3, 3)

   real(wp) :: detinv

   ! Calculate the inverse determinant of the matrix
   detinv = 1.0_wp/matdet_3x3(a)

   ! Calculate the inverse of the matrix
   b(1, 1) = +detinv * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
   b(2, 1) = -detinv * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
   b(3, 1) = +detinv * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
   b(1, 2) = -detinv * (a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
   b(2, 2) = +detinv * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
   b(3, 2) = -detinv * (a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
   b(1, 3) = +detinv * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
   b(2, 3) = -detinv * (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
   b(3, 3) = +detinv * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))

end function matinv_3x3


!> Determinat of 3×3 matrix
pure function matdet_3x3(a) result (det)

   !> Matrix
   real(wp), intent(in) :: a(3, 3)

   !> Determinant
   real(wp) :: det

   det =  a(1, 1) * a(2, 2) * a(3, 3)  &
      & - a(1, 1) * a(2, 3) * a(3, 2)  &
      & - a(1, 2) * a(2, 1) * a(3, 3)  &
      & + a(1, 2) * a(2, 3) * a(3, 1)  &
      & + a(1, 3) * a(2, 1) * a(3, 2)  &
      & - a(1, 3) * a(2, 2) * a(3, 1)

end function matdet_3x3


!> Implements the cross/vector product between two 3D vectors
pure function crossprod(a,b) result(c)

   !> First vector
   real(wp), intent(in) :: a(3)

   !> Second vector
   real(wp), intent(in) :: b(3)

   !> Orthogonal vector
   real(wp) :: c(3)

   c(1) = a(2) * b(3) - b(2) * a(3)
   c(2) = a(3) * b(1) - b(3) * a(1)
   c(3) = a(1) * b(2) - b(1) * a(2)

end function crossprod


!> Calculates eigenvalues based on the trigonometric solution of A = pB + qI
pure subroutine eigval_3x3(a, w)

   !> The symmetric input matrix
   real(wp), intent(in) :: a(3, 3)

   !> Contains eigenvalues on exit
   real(wp), intent(out) :: w(3)

   real(wp) :: q, p, r

   r = a(1, 2) * a(1, 2) + a(1, 3) * a(1, 3) + a(2, 3) * a(2, 3)
   q = (a(1, 1) + a(2, 2) + a(3, 3)) / 3.0_wp
   w(1) = a(1, 1) - q
   w(2) = a(2, 2) - q
   w(3) = a(3, 3) - q
   p = sqrt((w(1) * w(1) + w(2) * w(2) + w(3) * w(3) + 2*r) / 6.0_wp)
   r = (w(1) * (w(2) * w(3) - a(2, 3) * a(2, 3)) &
      & - a(1, 2) * (a(1, 2) * w(3) - a(2, 3) * a(1, 3)) &
      & + a(1, 3) * (a(1, 2) * a(2, 3) - w(2) * a(1, 3))) / (p*p*p) * 0.5_wp

   if (r <= -1.0_wp) then
      r = 0.5_wp * twothirdpi
   else if (r >= 1.0_wp) then
      r = 0.0_wp
   else
      r = acos(r) / 3.0_wp
   end if

   w(3) = q + 2 * p * cos(r)
   w(1) = q + 2 * p * cos(r + twothirdpi)
   w(2) = 3 * q - w(1) - w(3)

end subroutine eigval_3x3


!> Calculates eigenvector using an analytical method based on vector cross
!  products.
pure subroutine eigvec_3x3(a, w, q)
   real(wp), intent(inout) :: a(3,3)
   real(wp), intent(out) :: w(3)
   real(wp), intent(out) :: q(3,3)

   real(wp), parameter :: eps = epsilon(1.0_wp)
   real(wp) norm, n1, n2, n3, precon
   integer :: i

   w(1) = max(abs(a(1, 1)), abs(a(1, 2)))
   w(2) = max(abs(a(1, 3)), abs(a(2, 2)))
   w(3) = max(abs(a(2, 3)), abs(a(3, 3)))
   precon = max(w(1), max(w(2), w(3)))

   ! null matrix
   if (precon < eps) then
      w(1) = 0.0_wp
      w(2) = 0.0_wp
      w(3) = 0.0_wp
      q(1, 1) = 1.0_wp
      q(2, 2) = 1.0_wp
      q(3, 3) = 1.0_wp
      q(1, 2) = 0.0_wp
      q(1, 3) = 0.0_wp
      q(2, 3) = 0.0_wp
      q(2, 1) = 0.0_wp
      q(3, 1) = 0.0_wp
      q(3, 2) = 0.0_wp
      return
   end if

   norm = 1.0_wp / precon

   a(1, 1) = a(1, 1) * norm
   a(1, 2) = a(1, 2) * norm
   a(2, 2) = a(2, 2) * norm
   a(1, 3) = a(1, 3) * norm
   a(2, 3) = a(2, 3) * norm
   a(3, 3) = a(3, 3) * norm

   ! Calculate eigenvalues
   call eigval_3x3(a, w)

   ! Compute first eigenvector
   a(1, 1) = a(1, 1) - w(1)
   a(2, 2) = a(2, 2) - w(1)
   a(3, 3) = a(3, 3) - w(1)

   q(1, 1) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
   q(2, 1) = a(1, 3) * a(1, 2) - a(1, 1) * a(2, 3)
   q(3, 1) = a(1, 1) * a(2, 2) - a(1, 2) * a(1, 2)
   q(1, 2) = a(1, 2) * a(3, 3) - a(1, 3) * a(2, 3)
   q(2, 2) = a(1, 3) * a(1, 3) - a(1, 1) * a(3, 3)
   q(3, 2) = a(1, 1) * a(2, 3) - a(1, 2) * a(1, 3)
   q(1, 3) = a(2, 2) * a(3, 3) - a(2, 3) * a(2, 3)
   q(2, 3) = a(2, 3) * a(1, 3) - a(1, 2) * a(3, 3)
   q(3, 3) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
   n1 = q(1, 1) * q(1, 1) + q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)
   n2 = q(1, 2) * q(1, 2) + q(2, 2) * q(2, 2) + q(3, 2) * q(3, 2)
   n3 = q(1, 3) * q(1, 3) + q(2, 3) * q(2, 3) + q(3, 3) * q(3, 3)

   norm = n1
   i = 1
   if (n2 > norm) then
      i = 2
      norm = n1
   end if
   if (n3 > norm) then
      i = 3
   end if

   if (i == 1) then
      norm = sqrt(1.0_wp / n1)
      q(1, 1) = q(1, 1) * norm
      q(2, 1) = q(2, 1) * norm
      q(3, 1) = q(3, 1) * norm
   else if (i == 2) then
      norm = sqrt(1.0_wp / n2)
      q(1, 1) = q(1, 2) * norm
      q(2, 1) = q(2, 2) * norm
      q(3, 1) = q(3, 2) * norm
   else
      norm = sqrt(1.0_wp / n3)
      q(1, 1) = q(1, 3) * norm
      q(2, 1) = q(2, 3) * norm
      q(3, 1) = q(3, 3) * norm
   end if

   ! Robustly compute a right-hand orthonormal set (ev1, u, v)
   if (abs(q(1, 1)) > abs(q(2, 1))) then
      norm = sqrt(1.0_wp / (q(1, 1) * q(1, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = -q(3, 1) * norm
      q(2, 2) = 0.0_wp
      q(3, 2) = +q(1, 1) * norm
   else
      norm = sqrt(1.0_wp / (q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = 0.0_wp
      q(2, 2) = +q(3, 1) * norm
      q(3, 2) = -q(2, 1) * norm
   end if
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   ! Reset A
   a(1, 1) = a(1, 1) + w(1)
   a(2, 2) = a(2, 2) + w(1)
   a(3, 3) = a(3, 3) + w(1)

   ! A*U
   n1 = a(1, 1) * q(1, 2) + a(1, 2) * q(2, 2) + a(1, 3) * q(3, 2)
   n2 = a(1, 2) * q(1, 2) + a(2, 2) * q(2, 2) + a(2, 3) * q(3, 2)
   n3 = a(1, 3) * q(1, 2) + a(2, 3) * q(2, 2) + a(3, 3) * q(3, 2)

   ! A*V, note out of order computation
   a(3, 3) = a(1, 3) * q(1, 3) + a(2, 3) * q(2, 3) + a(3, 3) * q(3, 3)
   a(1, 3) = a(1, 1) * q(1, 3) + a(1, 2) * q(2, 3) + a(1, 3) * q(3, 3)
   a(2, 3) = a(1, 2) * q(1, 3) + a(2, 2) * q(2, 3) + a(2, 3) * q(3, 3)

   ! UT*(A*U) - l2*E
   n1 = q(1, 2) * n1 + q(2, 2) * n2 + q(3, 2) * n3 - w(2)
   ! UT*(A*V)
   n2 = q(1, 2) * a(1, 3) + q(2, 2) * a(2, 3) + q(3, 2) * a(3, 3)
   ! VT*(A*V) - l2*E
   n3 = q(1, 3) * a(1, 3) + q(2, 3) * a(2, 3) + q(3, 3) * a(3, 3) - w(2)

   if (abs(n1) >= abs(n3)) then
      norm = max(abs(n1), abs(n2))
      if (norm > eps) then
         if (abs(n1) >= abs(n2)) then
            n2 = n2 / n1
            n1 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n1
         else
            n1 = n1 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n1 * n1))
            n1 = n1 * n2
         end if
         q(1, 2) = n2 * q(1, 2) - n1 * q(1, 3)
         q(2, 2) = n2 * q(2, 2) - n1 * q(2, 3)
         q(3, 2) = n2 * q(3, 2) - n1 * q(3, 3)
      end if
   else
      norm = max(abs(n3), abs(n2))
      if (norm > eps) then
         if (abs(n3) >= abs(n2)) then
            n2 = n2 / n3
            n3 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n3
         else
            n3 = n3 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n3 * n3))
            n3 = n3 * n2
         end if
         q(1, 2) = n3 * q(1, 2) - n2 * q(1, 3)
         q(2, 2) = n3 * q(2, 2) - n2 * q(2, 3)
         q(3, 2) = n3 * q(3, 2) - n2 * q(3, 3)
      end if
   end if

   ! Calculate third eigenvector from cross product
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   w(1) = w(1) * precon
   w(2) = w(2) * precon
   w(3) = w(3) * precon

end subroutine eigvec_3x3


end module mctc_io_math

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Handle conversion between element symbols and atomic numbers
module mctc_io_symbols
   use mctc_io_resize, only : resize
   implicit none
   private

   public :: symbol_length
   public :: symbol_to_number, number_to_symbol, number_to_lcsymbol
   public :: to_number, to_symbol, to_lcsymbol
   public :: get_identity, collect_identical


   !> Get chemical identity
   interface get_identity
      module procedure :: get_identity_number
      module procedure :: get_identity_symbol
   end interface get_identity


   !> Maximum allowed length of element symbols
   integer, parameter :: symbol_length = 4


   !> Periodic system of elements
   character(len=2), parameter :: pse(118) = [ &
      & 'H ','He', &
      & 'Li','Be','B ','C ','N ','O ','F ','Ne', &
      & 'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
      & 'K ','Ca', &
      & 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
      &           'Ga','Ge','As','Se','Br','Kr', &
      & 'Rb','Sr', &
      & 'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', &
      &           'In','Sn','Sb','Te','I ','Xe', &
      & 'Cs','Ba', &
      & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
      & 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
      &           'Tl','Pb','Bi','Po','At','Rn', &
      & 'Fr','Ra', &
      & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No', &
      & 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn', &
      &           'Nh','Fl','Mc','Lv','Ts','Og' ]


   !> Lower case version of the periodic system of elements
   character(len=2), parameter :: lcpse(118) = [ &
      & 'h ','he', &
      & 'li','be','b ','c ','n ','o ','f ','ne', &
      & 'na','mg','al','si','p ','s ','cl','ar', &
      & 'k ','ca', &
      & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
      &           'ga','ge','as','se','br','kr', &
      & 'rb','sr', &
      & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd', &
      &           'in','sn','sb','te','i ','xe', &
      & 'cs','ba','la', &
      & 'ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
      & 'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
      &           'tl','pb','bi','po','at','rn', &
      & 'fr','ra','ac', &
      & 'th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no', &
      & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn', &
      &           'nh','fl','mc','lv','ts','og' ]


   !> ASCII offset between lowercase and uppercase letters
   integer, parameter :: offset = iachar('a') - iachar('A')


contains


!> Convert element symbol to atomic number
elemental subroutine symbol_to_number(number, symbol)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Atomic number
   integer, intent(out) :: number

   character(len=2) :: lcsymbol
   integer :: i, j, k, l

   number = 0
   lcsymbol = '  '

   k = 0
   do j = 1, len_trim(symbol)
      if (k > 2) exit
      l = iachar(symbol(j:j))
      if (k >= 1 .and. l == iachar(' ')) exit
      if (k >= 1 .and. l == 9) exit
      if (l >= iachar('A') .and. l <= iachar('Z')) l = l + offset
      if (l >= iachar('a') .and. l <= iachar('z')) then
         k = k+1
         if (k > 2) exit
         lcsymbol(k:k) = achar(l)
      endif
   enddo

   do i = 1, size(lcpse)
      if (lcsymbol == lcpse(i)) then
         number = i
         exit
      endif
   enddo

   if (number == 0) then
      select case(lcsymbol)
      case('d ', 't ')
         number = 1
      end select
   end if

end subroutine symbol_to_number


!> Convert atomic number to element symbol
elemental subroutine number_to_symbol(symbol, number)

   !> Atomic number
   integer, intent(in) :: number

   !> Element symbol
   character(len=2), intent(out) :: symbol

   if (number <= 0 .or. number > size(pse)) then
      symbol = '--'
   else
      symbol = pse(number)
   endif

end subroutine number_to_symbol


!> Convert atomic number to element symbol
elemental subroutine number_to_lcsymbol(symbol, number)

   !> Atomic number
   integer, intent(in) :: number

   !> Element symbol
   character(len=2), intent(out) :: symbol

   if (number <= 0 .or. number > size(lcpse)) then
      symbol = '--'
   else
      symbol = lcpse(number)
   endif

end subroutine number_to_lcsymbol


!> Convert element symbol to atomic number
elemental function to_number(symbol) result(number)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Atomic number
   integer :: number

   call symbol_to_number(number, symbol)

end function to_number


!> Convert atomic number to element symbol
elemental function to_symbol(number) result(symbol)

   !> Atomic number
   integer,intent(in) :: number

   !> Element symbol
   character(len=2) :: symbol

   call number_to_symbol(symbol, number)

end function to_symbol


!> Convert atomic number to element symbol
elemental function to_lcsymbol(number) result(symbol)

   !> Atomic number
   integer,intent(in) :: number

   !> Element symbol
   character(len=2) :: symbol

   call number_to_lcsymbol(symbol, number)

end function to_lcsymbol


!> Get chemical identity from a list of atomic numbers
pure subroutine get_identity_number(nid, identity, number)

   !> Number of unique species
   integer, intent(out) :: nid

   !> Ordinal numbers
   integer, intent(in) :: number(:)

   !> Chemical identity
   integer, intent(out) :: identity(:)

   integer, allocatable :: itmp(:)
   integer :: nat, iat, iid

   nat = size(identity)
   allocate(itmp(nat))
   nid = 0
   do iat = 1, nat
      iid = find_number(itmp(:nid), number(iat))
      if (iid == 0) then
         call append_number(itmp, nid, number(iat))
         iid = nid
      end if
      identity(iat) = iid
   end do

end subroutine get_identity_number


!> Get chemical identity from a list of element symbols
pure subroutine get_identity_symbol(nid, identity, symbol)

   !> Number of unique species
   integer, intent(out) :: nid

   !> Element symbols
   character(len=symbol_length), intent(in) :: symbol(:)

   !> Chemical identity
   integer, intent(out) :: identity(:)

   character(len=symbol_length), allocatable :: stmp(:)
   integer :: nat, iat, iid

   nat = size(identity)
   allocate(stmp(nat))
   nid = 0
   do iat = 1, nat
      iid = find_symbol(stmp(:nid), symbol(iat))
      if (iid == 0) then
         call append_symbol(stmp, nid, symbol(iat))
         iid = nid
      end if
      identity(iat) = iid
   end do

end subroutine get_identity_symbol


!> Establish a mapping between unique atom types and species
pure subroutine collect_identical(identity, mapping)

   !> Chemical identity
   integer, intent(in) :: identity(:)

   !> Mapping from unique atoms
   integer, intent(out) :: mapping(:)

   integer :: iid, iat

   do iid = 1, size(mapping)
      do iat = 1, size(identity)
         if (identity(iat) == iid) then
            mapping(iid) = iat
            exit
         end if
      end do
   end do

end subroutine collect_identical


!> Find element symbol in an unordered list, all entries are required to be unique
pure function find_symbol(list, symbol) result(position)

   !> List of element symbols
   character(len=*), intent(in) :: list(:)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Position of the symbol in list if found, otherwise zero
   integer :: position
   integer :: isym

   position = 0
   do isym = 1, size(list)
      if (symbol == list(isym)) then
         position = isym
         exit
      end if
   end do

end function find_symbol


!> Find atomic number in an unordered list, all entries are required to be unique
pure function find_number(list, number) result(position)

   !> List of atomic numbers
   integer, intent(in) :: list(:)

   !> Atomic number
   integer, intent(in) :: number

   !> Position of the number in list if found, otherwise zero
   integer :: position
   integer :: inum

   position = 0
   do inum = 1, size(list)
      if (number == list(inum)) then
         position = inum
         exit
      end if
   end do

end function find_number


!> Append an element symbol to an unsorted list, to ensure no dublicates search
!> for the element symbol first
pure subroutine append_symbol(list, nlist, symbol)

   !> List of element symbols
   character(len=*), allocatable, intent(inout) :: list(:)

   !> Current occupied size of list
   integer, intent(inout) :: nlist

   !> Elements symbol
   character(len=*), intent(in) :: symbol

   if (nlist >= size(list)) then
      call resize(list)
   end if

   nlist = nlist + 1
   list(nlist) = symbol

end subroutine append_symbol


!> Append an atomic number to an unsorted list, to ensure no dublicates search
!> for the atomic number first
pure subroutine append_number(list, nlist, number)

   !> List of atomic number
   integer, allocatable, intent(inout) :: list(:)

   !> Current occupied size of list
   integer, intent(inout) :: nlist

   !> Atomic number
   integer, intent(in) :: number

   if (nlist >= size(list)) then
      call resize(list)
   end if

   nlist = nlist + 1
   list(nlist) = number

end subroutine append_number


end module mctc_io_symbols

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Interface to BLAS library for matrix-vector and matrix-matrix operations
module multicharge_blas
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: symv, gemv, gemm


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface gemv
      module procedure :: mchrg_sgemv
      module procedure :: mchrg_dgemv
      module procedure :: mchrg_sgemv312
      module procedure :: mchrg_sgemv321
      module procedure :: mchrg_dgemv312
      module procedure :: mchrg_dgemv321
   end interface gemv

   !> Performs the matrix-vector  operation
   !>
   !>    y := alpha*A*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are n element vectors and
   !> A is an n by n symmetric matrix.
   interface symv
      module procedure :: mchrg_ssymv
      module procedure :: mchrg_dsymv
   end interface symv

   !> Performs one of the matrix-matrix operations
   !>
   !>    C := alpha*A*B + beta*C,
   !>
   !> or
   !>
   !>    C := alpha*B*A + beta*C,
   !>
   !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
   !> C are  m by n matrices.
   interface gemm
      module procedure :: mchrg_sgemm
      module procedure :: mchrg_dgemm
      module procedure :: mchrg_sgemm323
      module procedure :: mchrg_sgemm233
      module procedure :: mchrg_sgemm332
      module procedure :: mchrg_dgemm323
      module procedure :: mchrg_dgemm233
      module procedure :: mchrg_dgemm332
   end interface gemm


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface blas_gemv
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgemv
   end interface blas_gemv

   !> Performs the matrix-vector  operation
   !>
   !>    y := alpha*A*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are n element vectors and
   !> A is an n by n symmetric matrix.
   interface blas_symv
      pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ssymv
      pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dsymv
   end interface blas_symv

   !> Performs one of the matrix-matrix operations
   !>
   !>    C := alpha*op( A )*op( B ) + beta*C,
   !>
   !> where  op( X ) is one of
   !>
   !>    op( X ) = X   or   op( X ) = X**T,
   !>
   !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
   !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   interface blas_gemm
      pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine sgemm
      pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine dgemm
   end interface blas_gemm


contains


subroutine mchrg_sgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mchrg_sgemv312


subroutine mchrg_sgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mchrg_sgemv321


subroutine mchrg_dgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mchrg_dgemv312


subroutine mchrg_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mchrg_dgemv321


pure subroutine mchrg_sgemv(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_sgemv


pure subroutine mchrg_dgemv(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_dgemv


pure subroutine mchrg_ssymv(amat, xvec, yvec, uplo, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: ula
   real(sp) :: a, b
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_ssymv


pure subroutine mchrg_dsymv(amat, xvec, yvec, uplo, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: ula
   real(dp) :: a, b
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_dsymv


pure subroutine mchrg_sgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(sp) :: a, b
   integer :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine mchrg_sgemm


pure subroutine mchrg_dgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(dp) :: a, b
   integer :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine mchrg_dgemm


subroutine mchrg_sgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mchrg_sgemm323


subroutine mchrg_sgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mchrg_sgemm233


subroutine mchrg_sgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mchrg_sgemm332


subroutine mchrg_dgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mchrg_dgemm323


subroutine mchrg_dgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mchrg_dgemm233


subroutine mchrg_dgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mchrg_dgemm332


end module multicharge_blas

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_cutoff
   use mctc_env, only : wp
   implicit none
   private

   public :: get_lattice_points

   interface get_lattice_points
      module procedure :: get_lattice_points_cutoff
      module procedure :: get_lattice_points_rep_3d
   end interface get_lattice_points


contains


subroutine get_lattice_points_rep_3d(lat, rep, origin, trans)
   real(wp), intent(in) :: lat(:, :)
   integer, intent(in) :: rep(:)
   logical, intent(in) :: origin
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer :: itr, ix, iy, iz, jx, jy, jz

   itr = 0
   if (origin) then
      allocate(trans(3, product(2*rep+1)))
      do ix = 0, rep(1)
         do iy = 0, rep(2)
            do iz = 0, rep(3)
               do jx = 1, merge(-1, 1, ix > 0), -2
                  do jy = 1, merge(-1, 1, iy > 0), -2
                     do jz = 1, merge(-1, 1, iz > 0), -2
                        itr = itr + 1
                        trans(:, itr) = lat(:, 1)*ix*jx &
                           & + lat(:, 2)*iy*jy + lat(:, 3)*iz*jz
                     end do
                  end do
               end do
            end do
         end do
      end do
   else
      allocate(trans(3, product(2*rep+1)-1))
      do ix = 0, rep(1)
         do iy = 0, rep(2)
            do iz = 0, rep(3)
               if (ix == 0 .and. iy == 0 .and. iz == 0) cycle
               do jx = 1, merge(-1, 1, ix > 0), -2
                  do jy = 1, merge(-1, 1, iy > 0), -2
                     do jz = 1, merge(-1, 1, iz > 0), -2
                        itr = itr + 1
                        trans(:, itr) = lat(:, 1)*ix*jx &
                           & + lat(:, 2)*iy*jy + lat(:, 3)*iz*jz
                     end do
                  end do
               end do
            end do
         end do
      end do
   end if
end subroutine get_lattice_points_rep_3d


subroutine get_lattice_points_cutoff(periodic, lat, rthr, trans)
   logical, intent(in) :: periodic(:)
   real(wp), intent(in) :: rthr
   real(wp), intent(in) :: lat(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer :: rep(3)

   if (.not.any(periodic)) then
      allocate(trans(3, 1))
      trans(:, :) = 0.0_wp
   else
      call get_translations(lat, rthr, rep)
      call get_lattice_points(lat, rep, .true., trans)
   end if

end subroutine get_lattice_points_cutoff


!> generate a supercell based on a realspace cutoff, this subroutine
!  doesn't know anything about the convergence behaviour of the
!  associated property.
pure subroutine get_translations(lat, rthr, rep)
   real(wp), intent(in)  :: rthr
   real(wp), intent(in)  :: lat(3, 3)
   integer, intent(out) :: rep(3)
   real(wp) :: normx(3), normy(3), normz(3)
   real(wp) :: cos10, cos21, cos32

   ! find normal to the plane...
   call crossproduct(lat(:, 2), lat(:, 3), normx)
   call crossproduct(lat(:, 3), lat(:, 1), normy)
   call crossproduct(lat(:, 1), lat(:, 2), normz)
   ! ...normalize it...
   normx = normx/norm2(normx)
   normy = normy/norm2(normy)
   normz = normz/norm2(normz)
   ! cos angles between normals and lattice vectors
   cos10 = sum(normx*lat(:, 1))
   cos21 = sum(normy*lat(:, 2))
   cos32 = sum(normz*lat(:, 3))
   rep(1) = ceiling(abs(rthr/cos10))
   rep(2) = ceiling(abs(rthr/cos21))
   rep(3) = ceiling(abs(rthr/cos32))

contains

   pure subroutine crossproduct(a, b, c)
      real(wp), intent(in)  :: a(3), b(3)
      real(wp), intent(out) :: c(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
   end subroutine crossproduct

end subroutine get_translations


end module multicharge_cutoff

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_lapack
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: sytrf, sytrs, sytri

   interface sytrf
      module procedure :: mchrg_ssytrf
      module procedure :: mchrg_dsytrf
   end interface sytrf

   interface sytrs
      module procedure :: mchrg_ssytrs
      module procedure :: mchrg_ssytrs1
      module procedure :: mchrg_ssytrs3
      module procedure :: mchrg_dsytrs
      module procedure :: mchrg_dsytrs1
      module procedure :: mchrg_dsytrs3
   end interface sytrs

   interface sytri
      module procedure :: mchrg_ssytri
      module procedure :: mchrg_dsytri
   end interface sytri


   interface lapack_sytrf
      pure subroutine ssytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssytrf
      pure subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsytrf
   end interface lapack_sytrf

   interface lapack_sytrs
      pure subroutine ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ssytrs
      pure subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dsytrs
   end interface lapack_sytrs

   interface lapack_sytri
      pure subroutine ssytri(uplo, n, a, lda, ipiv, work, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(in) :: work(*)
      end subroutine ssytri
      pure subroutine dsytri(uplo, n, a, lda, ipiv, work, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(in) :: work(*)
      end subroutine dsytri
   end interface lapack_sytri


contains


subroutine mchrg_ssytrf(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, lwork, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, stat)
   if (stat == 0) then
      lwork = nint(test(1))
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
      else
         stat = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytrf failed"
   end if
end subroutine mchrg_ssytrf


subroutine mchrg_dsytrf(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, lwork, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, stat)
   if (stat == 0) then
      lwork = nint(test(1))
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
      else
         stat = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytrf failed"
   end if
end subroutine mchrg_dsytrf


subroutine mchrg_ssytrs(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, nrhs, lda, ldb
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytrs failed"
   end if
end subroutine mchrg_ssytrs


subroutine mchrg_dsytrs(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, nrhs, lda, ldb
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytrs failed"
   end if
end subroutine mchrg_dsytrs


subroutine mchrg_ssytrs1(amat, bvec, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_ssytrs1


subroutine mchrg_ssytrs3(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_ssytrs3


subroutine mchrg_dsytrs1(amat, bvec, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_dsytrs1


subroutine mchrg_dsytrs3(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_dsytrs3


subroutine mchrg_ssytri(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytri failed"
   end if
end subroutine mchrg_ssytri


subroutine mchrg_dsytri(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytri failed"
   end if
end subroutine mchrg_dsytri


end module multicharge_lapack

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_utils
   use mctc_env, only : wp
   use mctc_io_math, only : matinv_3x3
   implicit none

   public :: wrap_to_central_cell


contains


subroutine wrap_to_central_cell(xyz, lattice, periodic)
   real(wp), intent(inout) :: xyz(:, :)
   real(wp), intent(in) :: lattice(:, :)
   logical, intent(in) :: periodic(:)
   real(wp) :: invlat(3, 3), vec(3)
   integer :: iat, idir

   if (.not.any(periodic)) return

   invlat = matinv_3x3(lattice)
   do iat = 1, size(xyz, 2)
      vec(:) = matmul(invlat, xyz(:, iat))
      vec(:) = shift_back_abc(vec)
      xyz(:, iat) = matmul(lattice, vec)
   end do

end subroutine wrap_to_central_cell


elemental function shift_back_abc(in) result(out)
   !> fractional coordinate in (-∞,+∞)
   real(wp),intent(in) :: in
   !> fractional coordinate in [0,1)
   real(wp) :: out
   real(wp),parameter :: p_pbc_eps = 1.0e-14_wp
   out = in
   if(in < (0.0_wp - p_pbc_eps)) &
      out = in + real(ceiling(-in),wp)
   if(in > (1.0_wp + p_pbc_eps)) &
      out = in - real(floor  ( in),wp)
   if (abs(in - 1.0_wp) < p_pbc_eps) &
      out = in - 1.0_wp
end function shift_back_abc


end module dftd4_utils

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_data_covrad
   use mctc_env, only : wp
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_covalent_rad


   !> Covalent radii for DFT-D3 coordination number
   interface get_covalent_rad
      module procedure :: get_covalent_rad_num
      module procedure :: get_covalent_rad_sym
   end interface get_covalent_rad


   integer, parameter :: max_elem = 118

   !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
   !  188-197), values for metals decreased by 10 %
   real(wp), parameter :: covalent_rad_2009(max_elem) = aatoau * [ &
      & 0.32_wp,0.46_wp, & ! H,He
      & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
      & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
      & 1.76_wp,1.54_wp, & ! K,Ca
      &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
      &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
      &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
      & 1.89_wp,1.67_wp, & ! Rb,Sr
      &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
      &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
      &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
      & 2.09_wp,1.76_wp, & ! Cs,Ba
      &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
      &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
      &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
      &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
      &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
      & 2.01_wp,1.81_wp, & ! Fr,Ra
      &      1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
      &      1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
      &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
      &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
      &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

   !> D3 covalent radii used to construct the coordination number
   real(wp), parameter :: covalent_rad_d3(max_elem) = &
      & 4.0_wp / 3.0_wp * covalent_rad_2009

contains


!> Get covalent radius for a given element symbol
elemental function get_covalent_rad_sym(sym) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Covalent radius
   real(wp) :: rad

   rad = get_covalent_rad(to_number(sym))

end function get_covalent_rad_sym


!> Get covalent radius for a given atomic number
elemental function get_covalent_rad_num(num) result(rad)

   !> Atomic number
   integer, intent(in) :: num

   !> Covalent radius
   real(wp) :: rad

   if (num > 0 .and. num <= size(covalent_rad_d3)) then
      rad = covalent_rad_d3(num)
   else
      rad = 0.0_wp
   end if

end function get_covalent_rad_num


end module dftd4_data_covrad

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_data_en
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_electronegativity


   interface get_electronegativity
      module procedure :: get_electronegativity_num
      module procedure :: get_electronegativity_sym
   end interface get_electronegativity


   integer, parameter :: max_elem = 118

   !> Pauling electronegativities, used for the covalent coordination number.
   real(wp), parameter :: pauling_en(max_elem) = [ &
      & 2.20_wp,3.00_wp, & ! H,He
      & 0.98_wp,1.57_wp,2.04_wp,2.55_wp,3.04_wp,3.44_wp,3.98_wp,4.50_wp, & ! Li-Ne
      & 0.93_wp,1.31_wp,1.61_wp,1.90_wp,2.19_wp,2.58_wp,3.16_wp,3.50_wp, & ! Na-Ar
      & 0.82_wp,1.00_wp, & ! K,Ca
      &                 1.36_wp,1.54_wp,1.63_wp,1.66_wp,1.55_wp, & ! Sc-
      &                 1.83_wp,1.88_wp,1.91_wp,1.90_wp,1.65_wp, & ! -Zn
      &                 1.81_wp,2.01_wp,2.18_wp,2.55_wp,2.96_wp,3.00_wp, & ! Ga-Kr
      & 0.82_wp,0.95_wp, & ! Rb,Sr
      &                 1.22_wp,1.33_wp,1.60_wp,2.16_wp,1.90_wp, & ! Y-
      &                 2.20_wp,2.28_wp,2.20_wp,1.93_wp,1.69_wp, & ! -Cd
      &                 1.78_wp,1.96_wp,2.05_wp,2.10_wp,2.66_wp,2.60_wp, & ! In-Xe
      & 0.79_wp,0.89_wp, & ! Cs,Ba
      &         1.10_wp,1.12_wp,1.13_wp,1.14_wp,1.15_wp,1.17_wp,1.18_wp, & ! La-Eu
      &         1.20_wp,1.21_wp,1.22_wp,1.23_wp,1.24_wp,1.25_wp,1.26_wp, & ! Gd-Yb
      &                 1.27_wp,1.30_wp,1.50_wp,2.36_wp,1.90_wp, & ! Lu-
      &                 2.20_wp,2.20_wp,2.28_wp,2.54_wp,2.00_wp, & ! -Hg
      &                 1.62_wp,2.33_wp,2.02_wp,2.00_wp,2.20_wp,2.20_wp, & ! Tl-Rn
      ! only dummies below
      & 1.50_wp,1.50_wp, & ! Fr,Ra
      &         1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Ac-Am
      &         1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Cm-No
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-Cn
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp ] ! Nh-Og

contains


!> Get electronegativity for a given element symbol
elemental function get_electronegativity_sym(sym) result(en)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Electronegativity
   real(wp) :: en

   en = get_electronegativity(to_number(sym))

end function get_electronegativity_sym


!> Get electronegativity for a given atomic number
elemental function get_electronegativity_num(num) result(en)

   !> Atomic number
   integer, intent(in) :: num

   !> Electronegativity
   real(wp) :: en

   if (num > 0 .and. num <= size(pauling_en)) then
      en = pauling_en(num)
   else
      en = 0.0_wp
   end if

end function get_electronegativity_num


end module dftd4_data_en

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_data_hardness
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_hardness


   interface get_hardness
      module procedure :: get_hardness_num
      module procedure :: get_hardness_sym
   end interface get_hardness


   integer, parameter :: max_elem = 118

  !> Element-specific chemical hardnesses for the charge scaling function used
  !> to extrapolate the C6 coefficients in DFT-D4.
  real(wp), parameter :: chemical_hardness(max_elem) = [ &
    & 0.47259288_wp, 0.92203391_wp, 0.17452888_wp, 0.25700733_wp, 0.33949086_wp, &
    & 0.42195412_wp, 0.50438193_wp, 0.58691863_wp, 0.66931351_wp, 0.75191607_wp, &
    & 0.17964105_wp, 0.22157276_wp, 0.26348578_wp, 0.30539645_wp, 0.34734014_wp, &
    & 0.38924725_wp, 0.43115670_wp, 0.47308269_wp, 0.17105469_wp, 0.20276244_wp, &
    & 0.21007322_wp, 0.21739647_wp, 0.22471039_wp, 0.23201501_wp, 0.23933969_wp, &
    & 0.24665638_wp, 0.25398255_wp, 0.26128863_wp, 0.26859476_wp, 0.27592565_wp, &
    & 0.30762999_wp, 0.33931580_wp, 0.37235985_wp, 0.40273549_wp, 0.43445776_wp, &
    & 0.46611708_wp, 0.15585079_wp, 0.18649324_wp, 0.19356210_wp, 0.20063311_wp, &
    & 0.20770522_wp, 0.21477254_wp, 0.22184614_wp, 0.22891872_wp, 0.23598621_wp, &
    & 0.24305612_wp, 0.25013018_wp, 0.25719937_wp, 0.28784780_wp, 0.31848673_wp, &
    & 0.34912431_wp, 0.37976593_wp, 0.41040808_wp, 0.44105777_wp, 0.05019332_wp, &
    & 0.06762570_wp, 0.08504445_wp, 0.10247736_wp, 0.11991105_wp, 0.13732772_wp, &
    & 0.15476297_wp, 0.17218265_wp, 0.18961288_wp, 0.20704760_wp, 0.22446752_wp, &
    & 0.24189645_wp, 0.25932503_wp, 0.27676094_wp, 0.29418231_wp, 0.31159587_wp, &
    & 0.32902274_wp, 0.34592298_wp, 0.36388048_wp, 0.38130586_wp, 0.39877476_wp, &
    & 0.41614298_wp, 0.43364510_wp, 0.45104014_wp, 0.46848986_wp, 0.48584550_wp, &
    & 0.12526730_wp, 0.14268677_wp, 0.16011615_wp, 0.17755889_wp, 0.19497557_wp, &
    & 0.21240778_wp, 0.07263525_wp, 0.09422158_wp, 0.09920295_wp, 0.10418621_wp, &
    & 0.14235633_wp, 0.16394294_wp, 0.18551941_wp, 0.22370139_wp, 0.00000000_wp, &
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp]


contains


!> Get chemical hardness for a given element symbol
elemental function get_hardness_sym(sym) result(eta)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Chemical hardness
   real(wp) :: eta

   eta = get_hardness(to_number(sym))

end function get_hardness_sym


!> Get chemical hardness for a given atomic number
elemental function get_hardness_num(num) result(eta)

   !> Atomic number
   integer, intent(in) :: num

   !> Chemical hardness
   real(wp) :: eta

   if (num > 0 .and. num <= size(chemical_hardness)) then
      eta = chemical_hardness(num)
   else
      eta = 0.0_wp
   end if

end function get_hardness_num


end module dftd4_data_hardness

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_data_r4r2
   use mctc_env, only : wp
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_r4r2_val


   !> Get r4/r2 expectation value
   interface get_r4r2_val
      module procedure :: get_r4r2_val_num
      module procedure :: get_r4r2_val_sym
   end interface get_r4r2_val


   integer, parameter :: max_elem = 118

   !  PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
   !  rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
   !  He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638,
   !  Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
   !  not replaced but recalculated (PBE0/cc-pVQZ) were
   !   H: 8.0589 ->10.9359, Li:29.0974 ->39.7226, Be:14.8517 ->17.7460
   !  also new super heavies Cn,Nh,Fl,Lv,Og
   real(wp), parameter :: r4_over_r2(max_elem) = [  &
      &  8.0589_wp, 3.4698_wp, & ! H,He
      & 29.0974_wp,14.8517_wp,11.8799_wp, 7.8715_wp, 5.5588_wp, 4.7566_wp, 3.8025_wp, 3.1036_wp, & ! Li-Ne
      & 26.1552_wp,17.2304_wp,17.7210_wp,12.7442_wp, 9.5361_wp, 8.1652_wp, 6.7463_wp, 5.6004_wp, & ! Na-Ar
      & 29.2012_wp,22.3934_wp, & ! K,Ca
      &           19.0598_wp,16.8590_wp,15.4023_wp,12.5589_wp,13.4788_wp, & ! Sc-
      &           12.2309_wp,11.2809_wp,10.5569_wp,10.1428_wp, 9.4907_wp, & ! -Zn
      &                      13.4606_wp,10.8544_wp, 8.9386_wp, 8.1350_wp, 7.1251_wp, 6.1971_wp, & ! Ga-Kr
      & 30.0162_wp,24.4103_wp, & ! Rb,Sr
      &            20.3537_wp,17.4780_wp,13.5528_wp,11.8451_wp,11.0355_wp, & ! Y-
      &            10.1997_wp, 9.5414_wp, 9.0061_wp, 8.6417_wp, 8.9975_wp, & ! -Cd
      &                       14.0834_wp,11.8333_wp,10.0179_wp, 9.3844_wp, 8.4110_wp, 7.5152_wp, & ! In-Xe
      & 32.7622_wp,27.5708_wp, & ! Cs,Ba
      &            23.1671_wp,21.6003_wp,20.9615_wp,20.4562_wp,20.1010_wp,19.7475_wp,19.4828_wp, & ! La-Eu
      &            15.6013_wp,19.2362_wp,17.4717_wp,17.8321_wp,17.4237_wp,17.1954_wp,17.1631_wp, & ! Gd-Yb
      &            14.5716_wp,15.8758_wp,13.8989_wp,12.4834_wp,11.4421_wp, & ! Lu-
      &            10.2671_wp, 8.3549_wp, 7.8496_wp, 7.3278_wp, 7.4820_wp, & ! -Hg
      &                       13.5124_wp,11.6554_wp,10.0959_wp, 9.7340_wp, 8.8584_wp, 8.0125_wp, & ! Tl-Rn
      & 29.8135_wp,26.3157_wp, & ! Fr,Ra
      &            19.1885_wp,15.8542_wp,16.1305_wp,15.6161_wp,15.1226_wp,16.1576_wp, 0.0000_wp, & ! Ac-Am
      &             0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, & ! Cm-No
      &             0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, & ! Lr-
      &             0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 5.4929_wp, & ! -Cn
      &                        6.7286_wp, 6.5144_wp,10.9169_wp,10.3600_wp, 9.4723_wp, 8.6641_wp ] ! Nh-Og

   integer :: idum
   real(wp), parameter :: sqrt_z_r4_over_r2(max_elem) = &
      &  sqrt(0.5_wp*(r4_over_r2*[(sqrt(real(idum,wp)),idum=1,max_elem)]))


contains


!> Get r4/r2 expectation value for a given element symbol
elemental function get_r4r2_val_sym(sym) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> r4/r2 expectation value
   real(wp) :: rad

   rad = get_r4r2_val(to_number(sym))

end function get_r4r2_val_sym


!> Get r4/r2 expectation value for a given atomic number
elemental function get_r4r2_val_num(num) result(rad)

   !> Atomic number
   integer, intent(in) :: num

   !> r4/r2 expectation value
   real(wp) :: rad

   if (num > 0 .and. num <= size(sqrt_z_r4_over_r2)) then
      rad = sqrt_z_r4_over_r2(num)
   else
      rad = 0.0_wp
   end if

end function get_r4r2_val_num


end module dftd4_data_r4r2

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_data_zeff
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_effective_charge


   interface get_effective_charge
      module procedure :: get_effective_charge_num
      module procedure :: get_effective_charge_sym
   end interface get_effective_charge


   integer, parameter :: max_elem = 118


  !> Effective nuclear charges from the def2-ECPs used for calculating the reference
  !> polarizibilities for DFT-D4.
  real(wp), parameter :: effective_nuclear_charge(max_elem) = [ &
    &   1,                                                 2,  & ! H-He
    &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
    &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
    &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
    &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
    !  just copy & paste from above
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26] ! Rf-Og

contains


!> Get effective nuclear charge for a given element symbol
elemental function get_effective_charge_sym(sym) result(zeff)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Effective nuclear charge
   real(wp) :: zeff

   zeff = get_effective_charge(to_number(sym))

end function get_effective_charge_sym


!> Get effective nuclear charge for a given atomic number
elemental function get_effective_charge_num(num) result(zeff)

   !> Atomic number
   integer, intent(in) :: num

   !> Effective nuclear charge
   real(wp) :: zeff

   if (num > 0 .and. num <= size(effective_nuclear_charge)) then
      zeff = effective_nuclear_charge(num)
   else
      zeff = 0.0_wp
   end if

end function get_effective_charge_num


end module dftd4_data_zeff

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Basic structure representation of the system of interest
module mctc_io_structure
   use mctc_env_accuracy, only : wp
   use mctc_io_symbols, only : to_number, to_symbol, symbol_length, get_identity, &
      & collect_identical
   use mctc_io_structure_info, only : structure_info, pdb_data, sdf_data
   implicit none
   private

   public :: structure_type, new_structure, new


   !> Structure representation
   type :: structure_type

      !> Number of atoms
      integer :: nat = 0

      !> Number of unique species
      integer :: nid = 0

      !> Number of bonds
      integer :: nbd = 0

      !> Species identifier
      integer, allocatable :: id(:)

      !> Atomic number for each species
      integer, allocatable :: num(:)

      !> Element symbol for each species
      character(len=symbol_length), allocatable :: sym(:)

      !> Cartesian coordinates, in Bohr
      real(wp), allocatable :: xyz(:, :)

      !> Number of unpaired electrons
      integer :: uhf = 0

      !> Total charge
      real(wp) :: charge = 0.0_wp

      !> Lattice parameters
      real(wp), allocatable :: lattice(:, :)

      !> Periodic directions
      logical, allocatable :: periodic(:)

      !> Bond indices
      integer, allocatable :: bond(:, :)

      !> Comment, name or identifier for this structure
      character(len=:), allocatable :: comment

      !> Vendor specific structure annotations
      type(structure_info) :: info = structure_info()

      !> SDF atomic data annotations
      type(sdf_data), allocatable :: sdf(:)

      !> PDB atomic data annotations
      type(pdb_data), allocatable :: pdb(:)

   end type structure_type


   interface new
      module procedure :: new_structure
      module procedure :: new_structure_num
      module procedure :: new_structure_sym
   end interface


contains


!> Constructor for structure representations
subroutine new_structure(self, num, sym, xyz, charge, uhf, lattice, periodic, &
      & info, bond)

   !> Instance of the structure representation
   type(structure_type), intent(out) :: self

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Element symbols
   character(len=*), intent(in) :: sym(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Total charge
   real(wp), intent(in), optional :: charge

   !> Number of unpaired electrons
   integer, intent(in), optional :: uhf

   !> Lattice parameters
   real(wp), intent(in), optional :: lattice(:, :)

   !> Periodic directions
   logical, intent(in), optional :: periodic(:)

   !> Vendor specific structure information
   type(structure_info), intent(in), optional :: info

   !> Bond topology of the system
   integer, intent(in), optional :: bond(:, :)

   integer :: ndim, iid
   integer, allocatable :: map(:)

   ndim = min(size(num, 1), size(xyz, 2), size(sym, 1))

   self%nat = ndim
   allocate(self%id(ndim))
   allocate(self%xyz(3, ndim))

   if (present(lattice)) then
      self%lattice = lattice
   else
      allocate(self%lattice(0, 0))
   end if

   if (present(periodic)) then
      self%periodic = periodic
   else
      if (present(lattice)) then
         allocate(self%periodic(3))
         self%periodic(:) = .true.
      else
         allocate(self%periodic(1))
         self%periodic(:) = .false.
      end if
   end if

   call get_identity(self%nid, self%id, sym)
   allocate(map(self%nid))
   call collect_identical(self%id, map)

   allocate(self%num(self%nid))
   allocate(self%sym(self%nid))
   do iid = 1, self%nid
      self%num(iid) = num(map(iid))
      self%sym(iid) = sym(map(iid))
   end do
   self%xyz(:, :) = xyz(:, :ndim)

   if (present(charge)) then
      self%charge = charge
   else
      self%charge = 0.0_wp
   end if

   if (present(uhf)) then
      self%uhf = uhf
   else
      self%uhf = 0
   end if

   if (present(info)) then
      self%info = info
   else
      self%info = structure_info()
   end if

   if (present(bond)) then
      self%nbd = size(bond, 2)
      self%bond = bond
   end if

end subroutine new_structure


!> Simplified constructor for structure representations
subroutine new_structure_num(self, num, xyz, charge, uhf, lattice, periodic, &
      & info, bond)

   !> Instance of the structure representation
   type(structure_type), intent(out) :: self

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Total charge
   real(wp), intent(in), optional :: charge

   !> Number of unpaired electrons
   integer, intent(in), optional :: uhf

   !> Lattice parameters
   real(wp), intent(in), optional :: lattice(:, :)

   !> Periodic directions
   logical, intent(in), optional :: periodic(:)

   !> Vendor specific structure information
   type(structure_info), intent(in), optional :: info

   !> Bond topology of the system
   integer, intent(in), optional :: bond(:, :)

   integer :: ndim, iat
   character(len=symbol_length), allocatable :: sym(:)

   ndim = min(size(num, 1), size(xyz, 2))
   allocate(sym(ndim))
   do iat = 1, ndim
      sym(iat) = to_symbol(num(iat))
   end do

   call new_structure(self, num, sym, xyz, charge, uhf, lattice, periodic, &
      & info, bond)

end subroutine new_structure_num


!> Simplified constructor for structure representations
subroutine new_structure_sym(self, sym, xyz, charge, uhf, lattice, periodic, &
      & info, bond)

   !> Instance of the structure representation
   type(structure_type), intent(out) :: self

   !> Element symbols
   character(len=*), intent(in) :: sym(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Total charge
   real(wp), intent(in), optional :: charge

   !> Number of unpaired electrons
   integer, intent(in), optional :: uhf

   !> Lattice parameters
   real(wp), intent(in), optional :: lattice(:, :)

   !> Periodic directions
   logical, intent(in), optional :: periodic(:)

   !> Vendor specific structure information
   type(structure_info), intent(in), optional :: info

   !> Bond topology of the system
   integer, intent(in), optional :: bond(:, :)

   integer :: ndim, iat
   integer, allocatable :: num(:)

   ndim = min(size(sym, 1), size(xyz, 2))
   allocate(num(ndim))
   do iat = 1, ndim
      num(iat) = to_number(sym(iat))
   end do

   call new_structure(self, num, sym, xyz, charge, uhf, lattice, periodic, &
      & info, bond)

end subroutine new_structure_sym


end module mctc_io_structure

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_ewald
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_alpha

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   abstract interface
      !> Returns the max. value of a term in the reciprocal space part of the Ewald
      !> summation for a given vector length.
      pure function get_rec_term_gen(gg, alpha, vol) result(gTerm)
         import :: wp

         !> Length of the reciprocal space vector
         real(wp), intent(in) :: gg

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol

         !> Reciprocal term
         real(wp) :: gTerm

      end function get_rec_term_gen
   end interface

contains


subroutine get_alpha(lattice, alpha)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), intent(out) :: alpha
   real(wp) :: vol, rec_lat(3, 3)

   vol = abs(matdet_3x3(lattice))
   rec_lat = twopi*transpose(matinv_3x3(lattice))

   call search_alpha(lattice, rec_lat, vol, eps, alpha)

end subroutine get_alpha


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, alpha)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:,:)
   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:,:)
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance
   !> Optimal alpha
   real(wp), intent(out) :: alpha

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = 1.0e-8_wp
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:,:)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:,:)**2, dim=1)))

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   elseif (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, get_rec_term_3d, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha

!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume) result(diff)

   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha

   !> Procedure pointer to reciprocal routine
   procedure(get_rec_term_gen) :: get_rec_term

   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen

   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((get_rec_term(4*rlen, alpha, volume) &
      & - get_rec_term(5*rlen, alpha, volume))) &
      & - (get_dir_term(2*dlen, alpha) - get_dir_term(3*dlen, alpha))

end function rec_dir_diff

!> Returns the max. value of a term in the real space part of the Ewald summation
!> for a given vector length.
pure function get_dir_term(rr, alpha) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha*rr)/rr

end function get_dir_term


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function get_rec_term_3d(gg, alpha, vol) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*gg*gg/(alpha**2))/(vol*gg*gg))

end function get_rec_term_3d

end module multicharge_ewald

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_data_covrad
   use mctc_env, only : wp
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_covalent_rad


   !> Covalent radii for DFT-D3 coordination number
   interface get_covalent_rad
      module procedure :: get_covalent_rad_num
      module procedure :: get_covalent_rad_sym
   end interface get_covalent_rad


   integer, parameter :: max_elem = 118

   !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
   !  188-197), values for metals decreased by 10 %
   real(wp), parameter :: covalent_rad_2009(max_elem) = aatoau * [ &
      & 0.32_wp,0.46_wp, & ! H,He
      & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
      & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
      & 1.76_wp,1.54_wp, & ! K,Ca
      &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
      &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
      &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
      & 1.89_wp,1.67_wp, & ! Rb,Sr
      &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
      &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
      &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
      & 2.09_wp,1.76_wp, & ! Cs,Ba
      &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
      &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
      &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
      &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
      &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
      & 2.01_wp,1.81_wp, & ! Fr,Ra
      &      1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
      &      1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
      &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
      &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
      &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

   !> D3 covalent radii used to construct the coordination number
   real(wp), parameter :: covalent_rad_d3(max_elem) = &
      & 4.0_wp / 3.0_wp * covalent_rad_2009

contains


!> Get covalent radius for a given element symbol
elemental function get_covalent_rad_sym(sym) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Covalent radius
   real(wp) :: rad

   rad = get_covalent_rad(to_number(sym))

end function get_covalent_rad_sym


!> Get covalent radius for a given atomic number
elemental function get_covalent_rad_num(num) result(rad)

   !> Atomic number
   integer, intent(in) :: num

   !> Covalent radius
   real(wp) :: rad

   if (num > 0 .and. num <= size(covalent_rad_d3)) then
      rad = covalent_rad_d3(num)
   else
      rad = 0.0_wp
   end if

end function get_covalent_rad_num


end module multicharge_data_covrad

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Electronegativity equilibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
module multicharge_param_eeq2019
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_eeq_chi, get_eeq_eta, get_eeq_rad, get_eeq_kcn


   !> Element-specific electronegativity for the electronegativity equilibration charges.
   interface get_eeq_chi
      module procedure get_eeq_chi_sym
      module procedure get_eeq_chi_num
   end interface get_eeq_chi

   !> Element-specific chemical hardnesses for the electronegativity equilibration charges
   interface get_eeq_eta
      module procedure :: get_eeq_eta_sym
      module procedure :: get_eeq_eta_num
   end interface get_eeq_eta

   !> Element-specific CN scaling constant for the electronegativity equilibration charges.
   interface get_eeq_kcn
      module procedure :: get_eeq_kcn_sym
      module procedure :: get_eeq_kcn_num
   end interface get_eeq_kcn

   !> Element-specific charge widths for the electronegativity equilibration charges.
   interface get_eeq_rad
      module procedure :: get_eeq_rad_sym
      module procedure :: get_eeq_rad_num
   end interface get_eeq_rad


   !> Maximum atomic number allowed in EEQ calculations
   integer, parameter :: max_elem = 86


   !> Element-specific electronegativity for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_chi(max_elem) = [&
      & 1.23695041_wp, 1.26590957_wp, 0.54341808_wp, 0.99666991_wp, 1.26691604_wp, &
      & 1.40028282_wp, 1.55819364_wp, 1.56866440_wp, 1.57540015_wp, 1.15056627_wp, &
      & 0.55936220_wp, 0.72373742_wp, 1.12910844_wp, 1.12306840_wp, 1.52672442_wp, &
      & 1.40768172_wp, 1.48154584_wp, 1.31062963_wp, 0.40374140_wp, 0.75442607_wp, &
      & 0.76482096_wp, 0.98457281_wp, 0.96702598_wp, 1.05266584_wp, 0.93274875_wp, &
      & 1.04025281_wp, 0.92738624_wp, 1.07419210_wp, 1.07900668_wp, 1.04712861_wp, &
      & 1.15018618_wp, 1.15388455_wp, 1.36313743_wp, 1.36485106_wp, 1.39801837_wp, &
      & 1.18695346_wp, 0.36273870_wp, 0.58797255_wp, 0.71961946_wp, 0.96158233_wp, &
      & 0.89585296_wp, 0.81360499_wp, 1.00794665_wp, 0.92613682_wp, 1.09152285_wp, &
      & 1.14907070_wp, 1.13508911_wp, 1.08853785_wp, 1.11005982_wp, 1.12452195_wp, &
      & 1.21642129_wp, 1.36507125_wp, 1.40340000_wp, 1.16653482_wp, 0.34125098_wp, &
      & 0.58884173_wp, 0.68441115_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
      & 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
      & 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
      & 0.56999999_wp, 0.87936784_wp, 1.02761808_wp, 0.93297476_wp, 1.10172128_wp, &
      & 0.97350071_wp, 1.16695666_wp, 1.23997927_wp, 1.18464453_wp, 1.14191734_wp, &
      & 1.12334192_wp, 1.01485321_wp, 1.12950808_wp, 1.30804834_wp, 1.33689961_wp, &
      & 1.27465977_wp]

   !> Element-specific chemical hardnesses for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_eta(max_elem) = [&
      &-0.35015861_wp, 1.04121227_wp, 0.09281243_wp, 0.09412380_wp, 0.26629137_wp, &
      & 0.19408787_wp, 0.05317918_wp, 0.03151644_wp, 0.32275132_wp, 1.30996037_wp, &
      & 0.24206510_wp, 0.04147733_wp, 0.11634126_wp, 0.13155266_wp, 0.15350650_wp, &
      & 0.15250997_wp, 0.17523529_wp, 0.28774450_wp, 0.42937314_wp, 0.01896455_wp, &
      & 0.07179178_wp,-0.01121381_wp,-0.03093370_wp, 0.02716319_wp,-0.01843812_wp, &
      &-0.15270393_wp,-0.09192645_wp,-0.13418723_wp,-0.09861139_wp, 0.18338109_wp, &
      & 0.08299615_wp, 0.11370033_wp, 0.19005278_wp, 0.10980677_wp, 0.12327841_wp, &
      & 0.25345554_wp, 0.58615231_wp, 0.16093861_wp, 0.04548530_wp,-0.02478645_wp, &
      & 0.01909943_wp, 0.01402541_wp,-0.03595279_wp, 0.01137752_wp,-0.03697213_wp, &
      & 0.08009416_wp, 0.02274892_wp, 0.12801822_wp,-0.02078702_wp, 0.05284319_wp, &
      & 0.07581190_wp, 0.09663758_wp, 0.09547417_wp, 0.07803344_wp, 0.64913257_wp, &
      & 0.15348654_wp, 0.05054344_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
      & 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
      & 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
      & 0.11000000_wp,-0.02786741_wp, 0.01057858_wp,-0.03892226_wp,-0.04574364_wp, &
      &-0.03874080_wp,-0.03782372_wp,-0.07046855_wp, 0.09546597_wp, 0.21953269_wp, &
      & 0.02522348_wp, 0.15263050_wp, 0.08042611_wp, 0.01878626_wp, 0.08715453_wp, &
      & 0.10500484_wp]

   !> Element-specific CN scaling constant for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_kcn(max_elem) = [&
      & 0.04916110_wp, 0.10937243_wp,-0.12349591_wp,-0.02665108_wp,-0.02631658_wp, &
      & 0.06005196_wp, 0.09279548_wp, 0.11689703_wp, 0.15704746_wp, 0.07987901_wp, &
      &-0.10002962_wp,-0.07712863_wp,-0.02170561_wp,-0.04964052_wp, 0.14250599_wp, &
      & 0.07126660_wp, 0.13682750_wp, 0.14877121_wp,-0.10219289_wp,-0.08979338_wp, &
      &-0.08273597_wp,-0.01754829_wp,-0.02765460_wp,-0.02558926_wp,-0.08010286_wp, &
      &-0.04163215_wp,-0.09369631_wp,-0.03774117_wp,-0.05759708_wp, 0.02431998_wp, &
      &-0.01056270_wp,-0.02692862_wp, 0.07657769_wp, 0.06561608_wp, 0.08006749_wp, &
      & 0.14139200_wp,-0.05351029_wp,-0.06701705_wp,-0.07377246_wp,-0.02927768_wp, &
      &-0.03867291_wp,-0.06929825_wp,-0.04485293_wp,-0.04800824_wp,-0.01484022_wp, &
      & 0.07917502_wp, 0.06619243_wp, 0.02434095_wp,-0.01505548_wp,-0.03030768_wp, &
      & 0.01418235_wp, 0.08953411_wp, 0.08967527_wp, 0.07277771_wp,-0.02129476_wp, &
      &-0.06188828_wp,-0.06568203_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
      &-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
      &-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
      &-0.11000000_wp,-0.03585873_wp,-0.03132400_wp,-0.05902379_wp,-0.02827592_wp, &
      &-0.07606260_wp,-0.02123839_wp, 0.03814822_wp, 0.02146834_wp, 0.01580538_wp, &
      &-0.00894298_wp,-0.05864876_wp,-0.01817842_wp, 0.07721851_wp, 0.07936083_wp, &
      & 0.05849285_wp]

   !> Element-specific charge widths for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_rad(max_elem) = [&
      & 0.55159092_wp, 0.66205886_wp, 0.90529132_wp, 1.51710827_wp, 2.86070364_wp, &
      & 1.88862966_wp, 1.32250290_wp, 1.23166285_wp, 1.77503721_wp, 1.11955204_wp, &
      & 1.28263182_wp, 1.22344336_wp, 1.70936266_wp, 1.54075036_wp, 1.38200579_wp, &
      & 2.18849322_wp, 1.36779065_wp, 1.27039703_wp, 1.64466502_wp, 1.58859404_wp, &
      & 1.65357953_wp, 1.50021521_wp, 1.30104175_wp, 1.46301827_wp, 1.32928147_wp, &
      & 1.02766713_wp, 1.02291377_wp, 0.94343886_wp, 1.14881311_wp, 1.47080755_wp, &
      & 1.76901636_wp, 1.98724061_wp, 2.41244711_wp, 2.26739524_wp, 2.95378999_wp, &
      & 1.20807752_wp, 1.65941046_wp, 1.62733880_wp, 1.61344972_wp, 1.63220728_wp, &
      & 1.60899928_wp, 1.43501286_wp, 1.54559205_wp, 1.32663678_wp, 1.37644152_wp, &
      & 1.36051851_wp, 1.23395526_wp, 1.65734544_wp, 1.53895240_wp, 1.97542736_wp, &
      & 1.97636542_wp, 2.05432381_wp, 3.80138135_wp, 1.43893803_wp, 1.75505957_wp, &
      & 1.59815118_wp, 1.76401732_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
      & 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
      & 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
      & 1.63999999_wp, 1.47055223_wp, 1.81127084_wp, 1.40189963_wp, 1.54015481_wp, &
      & 1.33721475_wp, 1.57165422_wp, 1.04815857_wp, 1.78342098_wp, 2.79106396_wp, &
      & 1.78160840_wp, 2.47588882_wp, 2.37670734_wp, 1.76613217_wp, 2.66172302_wp, &
      & 2.82773085_wp]

contains



!> Get electronegativity for species with a given symbol
elemental function get_eeq_chi_sym(symbol) result(chi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> electronegativity
   real(wp) :: chi

   chi = get_eeq_chi(to_number(symbol))

end function get_eeq_chi_sym


!> Get electronegativity for species with a given atomic number
elemental function get_eeq_chi_num(number) result(chi)

   !> Atomic number
   integer, intent(in) :: number

   !> electronegativity
   real(wp) :: chi

   if (number > 0 .and. number <= size(eeq_chi, dim=1)) then
      chi = eeq_chi(number)
   else
      chi = -1.0_wp
   end if

end function get_eeq_chi_num


!> Get hardness for species with a given symbol
elemental function get_eeq_eta_sym(symbol) result(eta)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> hardness
   real(wp) :: eta

   eta = get_eeq_eta(to_number(symbol))

end function get_eeq_eta_sym


!> Get hardness for species with a given atomic number
elemental function get_eeq_eta_num(number) result(eta)

   !> Atomic number
   integer, intent(in) :: number

   !> hardness
   real(wp) :: eta

   if (number > 0 .and. number <= size(eeq_eta, dim=1)) then
      eta = eeq_eta(number)
   else
      eta = -1.0_wp
   end if

end function get_eeq_eta_num


!> Get CN scaling for species with a given symbol
elemental function get_eeq_kcn_sym(symbol) result(kcn)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> CN scaling
   real(wp) :: kcn

   kcn = get_eeq_kcn(to_number(symbol))

end function get_eeq_kcn_sym


!> Get CN scaling for species with a given atomic number
elemental function get_eeq_kcn_num(number) result(kcn)

   !> Atomic number
   integer, intent(in) :: number

   !> CN scaling
   real(wp) :: kcn

   if (number > 0 .and. number <= size(eeq_kcn, dim=1)) then
      kcn = eeq_kcn(number)
   else
      kcn = -1.0_wp
   end if

end function get_eeq_kcn_num


!> Get charge width for species with a given symbol
elemental function get_eeq_rad_sym(symbol) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> charge width
   real(wp) :: rad

   rad = get_eeq_rad(to_number(symbol))

end function get_eeq_rad_sym


!> Get charge width for species with a given atomic number
elemental function get_eeq_rad_num(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> Charge width
   real(wp) :: rad

   if (number > 0 .and. number <= size(eeq_rad, dim=1)) then
      rad = eeq_rad(number)
   else
      rad = -1.0_wp
   end if

end function get_eeq_rad_num


end module multicharge_param_eeq2019

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Element specific data needed for the DFT-D4 model
module dftd4_data
   use dftd4_data_covrad, only : get_covalent_rad
   use dftd4_data_en, only : get_electronegativity
   use dftd4_data_hardness, only : get_hardness
   use dftd4_data_r4r2, only : get_r4r2_val
   use dftd4_data_zeff, only : get_effective_charge
   implicit none
   public

end module dftd4_data

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_ctfile
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_convert, only : aatoau
   use mctc_io_structure, only : structure_type, new
   use mctc_io_structure_info, only : sdf_data, structure_info
   use mctc_io_symbols, only : to_number, symbol_length
   use mctc_io_utils, only : getline
   implicit none
   private

   public :: read_sdf, read_molfile


contains


subroutine read_sdf(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: line
   integer :: stat

   call read_molfile(self, unit, error)
   if (allocated(error)) return

   stat = 0
   do while(stat == 0)
      call getline(unit, line, stat)
      if (index(line, '$$$$') == 1) exit
   end do
   if (stat /= 0) then
      call fatal_error(error, "Failed while reading SDF key-value pairs")
      return
   end if

end subroutine read_sdf


subroutine read_molfile(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: line
   character(len=:), allocatable :: comment
   integer :: i, iatom, jatom, ibond, btype, atomtype
   integer :: stat, length, charge(2, 15)
   integer :: number_of_atoms, number_of_bonds, number_of_atom_lists, &
      &       chiral_flag, number_of_stext_entries, i999
   integer :: list4(4), list12(12)
   real(wp) :: x, y, z
   character(len=2) :: sdf_dim
   character(len=3) :: symbol
   character(len=5) :: v2000
   integer, parameter :: ccc_to_charge(0:7) = [0, +3, +2, +1, 0, -1, -2, -3]
   logical :: two_dim
   character(len=symbol_length), allocatable :: sym(:)
   type(sdf_data), allocatable :: sdf(:)
   type(structure_info) :: info
   real(wp), allocatable :: xyz(:, :)
   integer, allocatable :: bond(:, :)

   two_dim = .false.

   call getline(unit, comment, stat)
   call getline(unit, line, stat)
   read(line, '(20x, a2)', iostat=stat) sdf_dim
   if (stat == 0) then
      two_dim = sdf_dim == '2D' .or. sdf_dim == '2d'
   end if
   call getline(unit, line, stat)
   call getline(unit, line, stat)
   read(line, '(3i3, 3x, 2i3, 12x, i3, 1x, a5)', iostat=stat) &
      & number_of_atoms, number_of_bonds, number_of_atom_lists, &
      & chiral_flag, number_of_stext_entries, i999, v2000
   if (stat /= 0) then
      call fatal_error(error, "Cannot read header of molfile")
      return
   end if

   allocate(sdf(number_of_atoms))
   allocate(xyz(3, number_of_atoms))
   allocate(sym(number_of_atoms))

   do iatom = 1, number_of_atoms
      call getline(unit, line, stat)
      read(line, '(3f10.4, 1x, a3, i2, 11i3)', iostat=stat) &
         & x, y, z, symbol, list12
      if (stat /= 0) then
         call fatal_error(error, "Cannot read coordinates from connection table")
         return
      end if
      atomtype = to_number(symbol)
      if (atomtype == 0) then
         call fatal_error(error, "Unknown atom type '"//trim(symbol)//"' in connection table")
         return
      end if
      xyz(:, iatom) = [x, y, z] * aatoau
      sym(iatom) = trim(symbol)
      sdf(iatom)%isotope = list12(1)
      sdf(iatom)%charge = ccc_to_charge(list12(2)) ! drop doublet radical
      sdf(iatom)%hydrogens = list12(4)
      sdf(iatom)%valence = list12(6)
   end do

   allocate(bond(3, number_of_bonds))
   do ibond = 1, number_of_bonds
      call getline(unit, line, stat)
      read(line, '(7i3)', iostat=stat) &
         & iatom, jatom, btype, list4
      if (stat /= 0) then
         call fatal_error(error, "Cannot read topology from connection table")
         return
      end if
      bond(:, ibond) = [ibond, jatom, btype]
   end do

   do while(stat == 0)
      call getline(unit, line, stat)
      if (index(line, 'M  END') == 1) exit
      if (index(line, 'M  CHG') == 1) then
         read(line(7:9), *) length
         read(line(10:), '(*(1x, i3, 1x, i3))') (charge(:, i), i=1, length)
         do i = 1, length
            sdf(charge(1, i))%charge = charge(2, i)
         end do
      end if
   end do
   if (stat /= 0) then
      call fatal_error(error, "Cannot read connection table")
      return
   end if

   info = structure_info(two_dimensional=two_dim, &
      & missing_hydrogen=any(sdf%hydrogens > 1))
   call new(self, sym, xyz, charge=real(sum(sdf%charge), wp), info=info, bond=bond)
   call move_alloc(sdf, self%sdf)
   if (len(comment) > 0) self%comment = comment

end subroutine read_molfile


end module mctc_io_read_ctfile

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_gaussian
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_structure, only : structure_type, new
   implicit none
   private

   public :: read_gaussian_external


contains


subroutine read_gaussian_external(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, n, mode, chrg, spin, iat, ii
   integer, allocatable :: at(:)
   real(wp), allocatable :: xyz(:,:)
   real(wp) :: coord(3), q

   read(unit, '(4i10)', iostat=stat) n, mode, chrg, spin
   if (stat.ne.0) then
      call fatal_error(error, "Could not read number of atoms, check format!")
      return
   end if

   if (n <= 0) then
      call fatal_error(error, "Found no atoms, cannot work without atoms!")
      return
   end if

   allocate(xyz(3, n))
   allocate(at(n))

   ii = 0
   do while (ii < n)
      read(unit, '(i10, 4f20.12)', iostat=stat) iat, coord, q
      if (is_iostat_end(stat)) exit
      if (stat.ne.0) then
         call fatal_error(error, "Could not read geometry from Gaussian file")
         return
      end if
      if (iat > 0) then
         ii = ii+1
         at(ii) = iat
         xyz(:, ii) = coord
      else
         call fatal_error(error, "Invalid atomic number")
         return
      end if
   end do

   call new(self, at, xyz, charge=real(chrg, wp), uhf=spin)

   if (ii /= n) then
      call fatal_error(error, "Atom number missmatch in Gaussian file")
      return
   end if

end subroutine read_gaussian_external


end module mctc_io_read_gaussian

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_genformat
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_convert, only : aatoau
   use mctc_io_structure, only : structure_type, new
   use mctc_io_structure_info, only : structure_info
   use mctc_io_symbols, only : to_number, symbol_length
   use mctc_io_utils, only : getline
   implicit none
   private

   public :: read_genformat


contains


subroutine read_genformat(mol, unit, error)

   !> Instance of the molecular structure data
   type(structure_type),intent(out) :: mol

   !> File handle
   integer,intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: line
   integer :: natoms, nspecies, iatom, dummy, isp, ilat, stat, istart, iend
   logical :: cartesian, periodic
   real(wp) :: coord(3), lattice(3, 3)
   character(len=1) :: variant
   character(len=symbol_length), allocatable :: species(:), sym(:)
   real(wp), allocatable :: xyz(:, :), abc(:, :)
   type(structure_info) :: info

   call next_line(unit, line, stat)
   read(line, *, iostat=stat) natoms, variant
   if (stat /= 0 .or. natoms < 1) then
      call fatal_error(error, 'could not read number of atoms')
      return
   end if

   allocate(species(natoms))
   allocate(sym(natoms))
   allocate(xyz(3, natoms))
   allocate(abc(3, natoms))

   select case(variant)
   case('c', 'C')
      cartesian = .true.
      periodic = .false.
   case('s', 'S')
      cartesian = .true.
      periodic = .true.
   case('f', 'F')
      cartesian = .false.
      periodic = .true.
   case default
      call fatal_error(error, 'invalid input version')
      return
   endselect

   call next_line(unit, line, stat)
   istart = 1
   iend = 1
   isp = 0
   do while(iend < len_trim(line))
      istart = verify(line(iend:), ' ') - 1 + iend
      iend = scan(line(istart:), ' ') - 1 + istart
      if (iend < istart) iend = len_trim(line)
      isp = isp + 1
      species(isp) = trim(line(istart:iend))
   end do
   nspecies = isp
   if (any(to_number(species(:nspecies)) == 0)) then
      call fatal_error(error, 'unknown atom type present')
      return
   end if

   do iatom = 1, natoms
      call next_line(unit, line, stat)
      read(line, *, iostat=stat) dummy, isp, coord
      if (stat /= 0) then
         call fatal_error(error, 'could not read coordinates from file')
         return
      end if
      sym(iatom) = species(isp)
      if (cartesian) then
         xyz(:, iatom) = coord * aatoau
      else
         abc(:, iatom) = coord
      end if
   end do

   if (periodic) then
      call next_line(unit, line, stat)
      if (stat /= 0) then
         call fatal_error(error, 'missing lattice information')
         return
      end if
      do ilat = 1, 3
         call next_line(unit, line, stat)
         read(line, *, iostat=stat) coord
         if (stat /= 0) then
            call fatal_error(error, 'could not read lattice from file')
            return
         end if
         lattice(:, ilat) = coord * aatoau
      end do
      if (.not.cartesian) then
         xyz = matmul(lattice, abc)
      end if
      info = structure_info(cartesian=cartesian)
      call new(mol, sym, xyz, lattice=lattice, info=info)
   else
      call new(mol, sym, xyz)
   end if


contains

subroutine next_line(unit, line, stat)
   integer,intent(in) :: unit
   character(len=:), allocatable, intent(out) :: line
   integer, intent(out) :: stat
   integer :: ihash

   stat = 0
   do while(stat == 0)
      call getline(unit, line, stat)
      ihash = index(line, '#')
      if (ihash > 0) line = line(:ihash-1)
      if (len_trim(line) > 0) exit
   end do
   line = trim(adjustl(line))
end subroutine next_line

end subroutine read_genformat


end module mctc_io_read_genformat

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_pdb
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_convert, only : aatoau
   use mctc_io_resize, only : resize
   use mctc_io_symbols, only : to_number, symbol_length
   use mctc_io_structure, only : structure_type, new
   use mctc_io_structure_info, only : pdb_data, resize
   use mctc_io_utils, only : getline
   implicit none
   private

   public :: read_pdb


contains


subroutine read_pdb(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type),intent(out) :: self

   !> File handle
   integer,intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: pdb_format = &
      &  '(6x,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)'
   integer, parameter :: p_initial_size = 1000 ! this is going to be a proteine

   integer :: iatom, jatom, iresidue, try, stat, atom_type
   real(wp) :: occ, temp, coords(3)
   real(wp), allocatable :: xyz(:,:)
   character(len=4) :: a_charge
   character(len=:), allocatable :: line
   character(len=symbol_length), allocatable :: sym(:)
   type(pdb_data), allocatable :: pdb(:)

   allocate(sym(p_initial_size), source=repeat(' ', symbol_length))
   allocate(xyz(3, p_initial_size), source=0.0_wp)
   allocate(pdb(p_initial_size), source=pdb_data())

   iatom = 0
   iresidue = 0

   stat = 0
   do while(stat == 0)
      call getline(unit, line, stat)
      if (index(line, 'END') == 1) exit
      if (index(line, 'ATOM') == 1 .or. index(line, 'HETATM') == 1) then
         if (iatom >= size(xyz, 2)) call resize(xyz)
         if (iatom >= size(sym)) call resize(sym)
         if (iatom >= size(pdb)) call resize(pdb)
         iatom = iatom + 1
         pdb(iatom)%het = index(line, 'HETATM') == 1
         read(line, pdb_format) &
            & jatom, pdb(iatom)%name, pdb(iatom)%loc, pdb(iatom)%residue, &
            & pdb(iatom)%chains, pdb(iatom)%residue_number, pdb(iatom)%code, &
            & coords, occ, temp, pdb(iatom)%segid, sym(iatom), a_charge
         xyz(:,iatom) = coords * aatoau
         atom_type = to_number(sym(iatom))
         if (atom_type == 0) then
            try = scan(pdb(iatom)%name, 'HCNOSPF')
            if (try > 0) sym(iatom) = pdb(iatom)%name(try:try)//' '
            pdb(iatom)%charge = 0
         else
            read(a_charge(1:1), *, iostat=stat) pdb(iatom)%charge
            if (stat /= 0) then
               stat = 0
               pdb(iatom)%charge = 0
            else
               if (a_charge(2:2) == '-') pdb(iatom)%charge = -pdb(iatom)%charge
            end if
         end if         
      end if
   end do
   if (stat /= 0) then
      call fatal_error(error, "could not read in coordinates, last line was: '"//line//"'")
      return
   end if

   call new(self, sym(:iatom), xyz(:, :iatom))
   self%pdb = pdb(:iatom)
   self%charge = sum(pdb(:iatom)%charge)

   if (.not.all(self%num > 0)) then
      call fatal_error(error, "invalid atom type found")
      return
   end if

   ! since PDB is used for biomolecules, this is a sensible check (prevents GIGO)
   if (.not.any(self%num == 1)) then
      call fatal_error(error, "You get no calculation today, please add hydrogen atoms first")
      return
   end if

end subroutine read_pdb


end module mctc_io_read_pdb

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_turbomole
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use mctc_io_resize, only : resize
   use mctc_io_structure, only : structure_type, new
   use mctc_io_structure_info, only : structure_info
   use mctc_io_symbols, only : to_number, symbol_length
   use mctc_io_utils, only : getline
   implicit none
   private

   public :: read_coord


   logical, parameter :: debug = .false.


contains


subroutine read_coord(mol, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: mol

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character, parameter :: flag = '$'
   integer, parameter :: p_initial_size = 100
   integer, parameter :: p_nlv(3) = [1, 4, 9], p_ncp(3) = [1, 3, 6]

   logical :: has_coord, has_periodic, has_lattice, has_cell
   logical :: cartesian, coord_in_bohr, lattice_in_bohr, pbc(3)
   integer :: stat, iatom, i, j, natoms, periodic, cell_vectors
   real(wp) :: latvec(9), conv, cellpar(6), lattice(3, 3)
   real(wp), allocatable :: coord(:, :), xyz(:, :)
   character(len=:), allocatable :: line, cell_string, lattice_string
   character(len=symbol_length), allocatable :: sym(:)
   type(structure_info) :: info

   allocate(sym(p_initial_size), source=repeat(' ', symbol_length))
   allocate(coord(3, p_initial_size), source=0.0_wp)

   iatom = 0
   periodic = 0
   cell_vectors = 0
   has_coord = .false.
   has_periodic = .false.
   has_lattice = .false.
   has_cell = .false.
   cartesian = .true.
   coord_in_bohr = .true.
   lattice_in_bohr = .true.
   lattice = 0.0_wp
   pbc = .false.

   stat = 0
   do while(stat == 0)
      call getline(unit, line, stat)
      if (index(line, flag) == 1) then
         if (index(line, 'end') == 2) then
            exit

         else if (.not.has_coord .and. index(line, 'coord') == 2) then
            has_coord = .true.
            ! $coord angs / $coord bohr / $coord frac
            call select_unit(line, coord_in_bohr, cartesian)
            coord_group: do while(stat == 0)
               call getline(unit, line, stat)
               if (index(line, flag) == 1) then
                  backspace(unit)
                  exit coord_group
               end if
               if (iatom >= size(coord, 2)) call resize(coord)
               if (iatom >= size(sym)) call resize(sym)
               iatom = iatom + 1
               read(line, *, iostat=stat) coord(:, iatom), sym(iatom)
            end do coord_group

         else if (.not.has_periodic .and. index(line, 'periodic') == 2) then
            has_periodic = .true.
            ! $periodic 0/1/2/3
            read(line(10:), *, iostat=stat) periodic

         else if (.not.has_lattice .and. index(line, 'lattice') == 2) then
            has_lattice = .true.
            ! $lattice bohr / $lattice angs
            call select_unit(line, lattice_in_bohr)
            cell_vectors = 0
            lattice_string = ''
            lattice_group: do while(stat == 0)
               call getline(unit, line, stat)
               if (index(line, flag) == 1) then
                  backspace(unit)
                  exit lattice_group
               end if
               cell_vectors = cell_vectors + 1
               lattice_string = lattice_string // ' ' // line
            end do lattice_group

         else if (.not.has_cell .and. index(line, 'cell') == 2) then
            has_cell = .true.
            ! $cell bohr / $cell angs
            call select_unit(line, lattice_in_bohr)
            call getline(unit, cell_string, stat)
            if (debug) print*, cell_string

         end if
      end if
   end do

   if (.not.has_coord .or. iatom == 0) then
      call fatal_error(error, "coordinates not present, cannot work without coordinates")
      return
   end if

   if (has_cell .and. has_lattice) then
      call fatal_error(error, "both lattice and cell group are present")
      return
   end if

   if (.not.has_periodic .and. (has_cell .or. has_lattice)) then
      call fatal_error(error, "cell and lattice definition is present, but periodicity is not given")
      return
   end if

   if (periodic > 0 .and. .not.(has_cell .or. has_lattice)) then
      call fatal_error(error, "system is periodic but definition of lattice/cell is missing")
      return
   end if

   if (.not.cartesian .and. periodic == 0) then
      call fatal_error(error, "fractional coordinates do not work for molecular systems")
      return
   end if

   natoms = iatom
   allocate(xyz(3, natoms))

   if (any(to_number(sym(:natoms)) == 0)) then
      call fatal_error(error, "unknown element present")
      return
   end if

   if (periodic > 0) pbc(:periodic) = .true.

   if (has_cell) then
      read(cell_string, *, iostat=stat) latvec(:p_ncp(periodic))
      if (debug) print*, latvec(:p_ncp(periodic))
      if (lattice_in_bohr) then
         conv = 1.0_wp
      else
         conv = aatoau
      end if
      select case(periodic)
      case(1)
         cellpar = [latvec(1)*conv, 1.0_wp, 1.0_wp, &
            &       pi/2, pi/2, pi/2]
      case(2)
         cellpar = [latvec(1)*conv, latvec(2)*conv, 1.0_wp, &
            &       pi/2, pi/2, latvec(3)*pi/180.0_wp]
      case(3)
         cellpar = [latvec(1:3)*conv, latvec(4:6)*pi/180.0_wp]
      end select
      call cell_to_dlat(cellpar, lattice)
   end if

   if (has_lattice) then
      if (cell_vectors /= periodic) then
         call fatal_error(error, "number of cell vectors does not match periodicity")
         return
      end if
      read(lattice_string, *, iostat=stat) latvec(:p_nlv(periodic))
      if (lattice_in_bohr) then
         conv = 1.0_wp
      else
         conv = aatoau
      end if
      j = 0
      do i = 1, periodic
         lattice(:periodic,  i) = latvec(j+1:j+periodic) * conv
         j = j + periodic
      end do
   end if

   if (cartesian) then
      if (coord_in_bohr) then
         conv = 1.0_wp
      else
         conv = aatoau
      end if
      xyz(:, :) = coord(:, :natoms) * conv
   else
      xyz = matmul(lattice, coord)
   end if

   ! save data on input format
   info = structure_info(cartesian=cartesian, lattice=has_lattice, &
      & angs_lattice=.not.lattice_in_bohr, angs_coord=.not.coord_in_bohr)
   call new(mol, sym(:natoms), xyz, lattice=lattice, periodic=pbc, info=info)

contains

   subroutine select_unit(line, in_bohr, cartesian)
      character(len=*), intent(in) :: line
      logical, intent(out) :: in_bohr
      logical, intent(out), optional :: cartesian
      in_bohr = index(line, ' angs') == 0
      if (present(cartesian)) cartesian = index(line, ' frac') == 0
   end subroutine select_unit

end subroutine read_coord


!> Calculate the lattice vectors from a set of cell parameters
pure subroutine cell_to_dlat(cellpar, lattice)

   !> Cell parameters
   real(wp), intent(in)  :: cellpar(6)

   !> Direct lattice
   real(wp), intent(out) :: lattice(3, 3)

   real(wp) :: dvol

   dvol = cell_to_dvol(cellpar)

   associate(alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
         &   alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6))

      lattice(1, 1) = alen
      lattice(2, 1) = 0.0_wp
      lattice(3, 1) = 0.0_wp
      lattice(3, 2) = 0.0_wp
      lattice(1, 2) = blen*cos(gam)
      lattice(2, 2) = blen*sin(gam)
      lattice(1, 3) = clen*cos(bet)
      lattice(2, 3) = clen*(cos(alp) - cos(bet)*cos(gam))/sin(gam);
      lattice(3, 3) = dvol/(alen*blen*sin(gam))

   end associate

end subroutine cell_to_dlat


!> Calculate the cell volume from a set of cell parameters
pure function cell_to_dvol(cellpar) result(dvol)

   !> Cell parameters
   real(wp), intent(in) :: cellpar(6)

   !> Cell volume
   real(wp) :: dvol

   real(wp) :: vol2

   associate(alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
         &   alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

      vol2 = 1.0_wp - cos(alp)**2 - cos(bet)**2 - cos(gam)**2 &
         & + 2.0_wp*cos(alp)*cos(bet)*cos(gam)

      dvol = sqrt(abs(vol2))*alen*blen*clen
      ! return negative volume instead of imaginary one (means bad cell parameters)
      if (vol2 < 0.0_wp) dvol = -dvol ! this should not happen, but who knows...

   end associate
end function cell_to_dvol


end module mctc_io_read_turbomole

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_vasp
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_convert, only : aatoau
   use mctc_io_resize, only : resize
   use mctc_io_structure, only : structure_type, new
   use mctc_io_structure_info, only : structure_info
   use mctc_io_symbols, only : to_number, symbol_length
   use mctc_io_utils, only : getline
   implicit none
   private

   public :: read_vasp


   logical, parameter :: debug = .false.


contains


subroutine read_vasp(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   logical :: selective, cartesian
   integer :: i, j, k, nn, ntype, natoms, izp, stat
   integer, allocatable :: ncount(:)
   real(wp) :: ddum, latvec(3), scalar, coord(3), lattice(3, 3)
   real(wp), allocatable :: xyz(:, :)
   character(len=:), allocatable :: line, comment
   character(len=2*symbol_length), allocatable :: args(:), args2(:)
   character(len=symbol_length), allocatable :: sym(:)
   type(structure_info) :: info

   selective = .false. ! Selective dynamics
   cartesian = .true.  ! Cartesian or direct
   lattice = 0
   stat = 0

   ntype = 0
   ! first line contains the symbols of different atom types
   call getline(unit, line, stat)
   if (stat /= 0) then
      call fatal_error(error, "Unexpected end of input encountered")
      return
   end if
   if (debug) print'(">", a)', line

   call parse_line(line, args, ntype)
   call move_alloc(line, comment)

   ! this line contains the global scaling factor,
   call getline(unit, line, stat)
   if (stat /= 0) then
      call fatal_error(error, "Unexpected end of input encountered")
      return
   end if
   if (debug) print'(">", a)', line
   read(line, *, iostat=stat) ddum
   if (stat /= 0) then
      call fatal_error(error, "Cannot read scaling factor from input")
      return
   end if
   ! the Ang->au conversion is included in the scaling factor
   if (debug) print'("->", g0)', ddum
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i = 1, 3
      call getline(unit, line, stat)
      if (stat /= 0) then
         call fatal_error(error, "Unexpected end of lattice vectors encountered")
         return
      end if
      if (debug) print'("->", a)', line
      read(line, *, iostat=stat) latvec
      if (stat /= 0) then
         call fatal_error(error, "Cannot read lattice vectors from input")
         return
      end if
      lattice(:, i) = latvec * scalar
   end do
   ! Either here are the numbers of each element,
   ! or (>vasp.5.1) here are the element symbols
   call getline(unit, line, stat)
   if (stat /= 0) then
      call fatal_error(error, "Unexpected end of input encountered")
      return
   end if
   if (debug) print'(">", a)', line

   ! try to verify that first element is actually a number
   i = max(verify(line, ' '), 1)
   j = scan(line(i:), ' ') - 2 + i
   if (j < i) j = len_trim(line)

   ! CONTCAR files have additional Element line here since vasp.5.1
   if (verify(line(i:j), '1234567890') /= 0) then
      call parse_line(line, args, ntype)
      call getline(unit, line, stat)
      if (debug) print'("->", a)', line
      if (stat /= 0) then
         call fatal_error(error, "Unexpected end of input encountered")
         return
      end if
   else
      deallocate(comment)
   end if
   call parse_line(line, args2, nn)
   if (nn /= ntype) then
      call fatal_error(error, 'Number of atom types mismatches the number of counts')
      return
   end if

   allocate(ncount(nn), source = 0)
   do i = 1, nn
      read(args2(i), *, iostat=stat) ncount(i)
      izp = to_number(args(i))
      if (izp < 1 .or. ncount(i) < 1) then
         call fatal_error(error, "Unknown element '"//trim(args(i))//"' encountered")
         return
      end if
   end do

   natoms = sum(ncount)
   allocate(sym(natoms))
   allocate(xyz(3, natoms))

   k = 0
   do i = 1, nn
      do j = 1, ncount(i)
         k = k+1
         sym(k) = trim(args(i))
      end do
   end do

   call getline(unit, line, stat)
   if (stat /= 0) then
      call fatal_error(error, "Could not read POSCAR")
      return
   end if
   if (debug) print'(">", a)', line
   line = adjustl(line)
   if (line(:1).eq.'s' .or. line(:1).eq.'S') then
      selective = .true.
      call getline(unit, line, stat)
      if (debug) print'("->", a)', line
      if (stat /= 0) then
         call fatal_error(error, "Unexpected end of input encountered")
         return
      end if
      line = adjustl(line)
   end if

   cartesian = (line(:1).eq.'c' .or. line(:1).eq.'C' .or. &
      &         line(:1).eq.'k' .or. line(:1).eq.'K')
   do i = 1, natoms
      call getline(unit, line, stat)
      if (stat /= 0) then
         call fatal_error(error, "Unexpected end of geometry encountered")
         return
      end if
      if (debug) print'("-->", a)', line
      read(line, *, iostat=stat) coord
      if (stat /= 0) then
         call fatal_error(error, "Cannot read geometry from input")
         return
      end if

      if (cartesian) then
         xyz(:, i) = coord*scalar
      else
         xyz(:, i) = matmul(lattice, coord)
      end if

   end do

   ! save information about this POSCAR for later
   info = structure_info(scale=ddum, selective=selective, cartesian=cartesian)
   call new(self, sym, xyz, lattice=lattice, info=info)
   if (allocated(comment)) self%comment = comment

end subroutine read_vasp


subroutine parse_line(line, args, nargs)
   character(len=*), intent(in) :: line
   character(len=2*symbol_length), allocatable, intent(out) :: args(:)
   integer, intent(out) :: nargs
   integer, parameter :: p_initial_size = 50
   integer :: istart, iend
   allocate(args(p_initial_size), source=repeat(' ', 2*symbol_length))
   istart = 1
   iend = 1
   nargs = 0
   do while(iend < len_trim(line))
      istart = verify(line(iend:), ' ') - 1 + iend
      iend = scan(line(istart:), ' ') - 1 + istart
      if (iend < istart) iend = len_trim(line)
      if (nargs >= size(args)) then
         call resize(args)
      end if
      nargs = nargs + 1
      args(nargs) = trim(line(istart:iend))
   end do
end subroutine parse_line


end module mctc_io_read_vasp

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read_xyz
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_convert, only : aatoau
   use mctc_io_structure, only : structure_type, new
   use mctc_io_symbols, only : to_number, to_symbol, symbol_length
   use mctc_io_utils, only : getline
   implicit none
   private

   public :: read_xyz


contains


subroutine read_xyz(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, n, iat, stat
   real(wp) :: x, y, z, conv
   real(wp), allocatable :: xyz(:, :)
   character(len=symbol_length) :: chdum
   character(len=symbol_length), allocatable :: sym(:)
   character(len=:), allocatable :: line, comment

   conv = aatoau

   read(unit, *, iostat=stat) n
   if (stat /= 0) then
      call fatal_error(error, "Could not read number of atoms, check format!")
      return
   end if

   if (n.lt.1) then
      call fatal_error(error, "Found no atoms, cannot work without atoms!")
      return
   end if

   allocate(sym(n))
   allocate(xyz(3, n))

   ! next record is a comment
   call getline(unit, comment, stat)
   if (stat /= 0) then
      call fatal_error(error, "Unexpected end of file")
      return
   end if

   ii = 0
   do while (ii < n)
      call getline(unit, line, stat)
      if (is_iostat_end(stat)) exit
      if (stat /= 0) then
         call fatal_error(error, "Could not read geometry from xyz file")
         return
      end if
      read(line, *, iostat=stat) chdum, x, y, z
      if (stat /= 0) then
         call fatal_error(error, "Could not parse coordinates from xyz file")
         return
      end if

      iat = to_number(chdum)
      if (iat <= 0) then
         read(chdum, *, iostat=stat) iat
         if (stat == 0) then
            chdum = to_symbol(iat)
         else
            iat = 0
         end if
      end if
      if (iat > 0) then
         ii = ii+1
         sym(ii) = trim(chdum)
         xyz(:, ii) = [x, y, z]*conv
      else
         call fatal_error(error, "Unknown element symbol: '"//trim(chdum)//"'")
         return
      end if
   end do

   if (ii /= n) then
      call fatal_error(error, "Atom number missmatch in xyz file")
      return
   end if

   call new(self, sym, xyz)
   if (len(comment) > 0) self%comment = comment

end subroutine read_xyz


end module mctc_io_read_xyz

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_ctfile
   use mctc_env_accuracy, only : wp
   use mctc_io_convert, only : autoaa
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_molfile, write_sdf


contains


subroutine write_sdf(self, unit, energy, gnorm)
   class(structure_type), intent(in) :: self
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gnorm
   !type(tb_buffer) :: sd_values
   character(len=:), allocatable :: line
   character(len=*), parameter :: sd_format = &
      & '("> <", a, ">", /, f20.12, /)'

   call write_molfile(self, unit)

!   sd_values = self%info
!   call sd_values%reset
!   do while(sd_values%next())
!      call sd_values%getline(line)
!      write(unit, '(a)') line
!   enddo

   if (present(energy)) then
      write(unit, sd_format) "total energy / Eh", energy
   endif

   if (present(gnorm)) then
      write(unit, sd_format) "gradient norm / Eh/a0", gnorm
   endif

   write(unit, '("$$$$")')

end subroutine write_sdf


subroutine write_molfile(self, unit, comment_line)
   class(structure_type), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: comment_line
   integer, parameter :: list4(4) = 0
   integer :: iatom, ibond, iatoms(3), list12(12)
   logical :: has_sdf_data
   integer, parameter :: charge_to_ccc(-3:3) = [7, 6, 5, 0, 3, 2, 1]
   character(len=8)  :: date
   character(len=10) :: time

   call date_and_time(date, time)

   if (present(comment_line)) then
      write(unit, '(a)') comment_line
   else
      if (allocated(self%comment)) then
         write(unit, '(a)') self%comment
      else
         write(unit, '(a)')
      end if
   end if
   write(unit, '(2x, 3x, 5x, 3a2, a4, "3D")') &
      &  date(5:6), date(7:8), date(3:4), time(:4)
   write(unit, '(a)')
   write(unit, '(3i3, 3x, 2i3, 12x, i3, 1x, a5)') &
      &  self%nat, self%nbd, 0, 0, 0, 999, 'V2000'

   has_sdf_data = allocated(self%sdf)

   do iatom = 1, self%nat
      if (has_sdf_data) then
         list12 = [self%sdf(iatom)%isotope, 0, 0, 0, 0, self%sdf(iatom)%valence, &
            & 0, 0, 0, 0, 0, 0]
      else
         list12 = 0
      endif
      write(unit, '(3f10.4, 1x, a3, i2, 11i3)') &
         & self%xyz(:, iatom)*autoaa, self%sym(self%id(iatom)), list12
   enddo

   if (self%nbd > 0) then
      if (size(self%bond, 1) > 2) then
         do ibond = 1, self%nbd
            write(unit, '(7i3)') self%bond(:3, ibond), list4
         end do
      else
         do ibond = 1, self%nbd
            write(unit, '(7i3)') self%bond(:2, ibond), 1, list4
         end do
      end if
   end if

   if (has_sdf_data) then
      if (sum(self%sdf%charge) /= nint(self%charge)) then
         write(unit, '(a, *(i3, 1x, i3, 1x, i3))') "M  CHG", 1, 1, nint(self%charge)
      else
         do iatom = 1, self%nat
            if (self%sdf(iatom)%charge /= 0) then
               write(unit, '(a, *(i3, 1x, i3, 1x, i3))') &
                  & "M  CHG", 1, iatom, self%sdf(iatom)%charge
            end if
         end do
      end if
   else
      if (nint(self%charge) /= 0) then
         write(unit, '(a, *(i3, 1x, i3, 1x, i3))') "M  CHG", 1, 1, nint(self%charge)
      end if
   end if

   write(unit, '(a)') "M  END"

end subroutine write_molfile


end module mctc_io_write_ctfile

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_gaussian
   use mctc_env_accuracy, only : wp
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_gaussian_external


contains


subroutine write_gaussian_external(mol, unit)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: iat

   write(unit, '(4i10)') mol%nat, 1, nint(mol%charge), mol%uhf
   do iat = 1, mol%nat
      write(unit, '(i10,4f20.12)') mol%num(mol%id(iat)), mol%xyz(:, iat), 0.0_wp
   end do

end subroutine write_gaussian_external


end module mctc_io_write_gaussian

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_genformat
   use mctc_env_accuracy, only : wp
   use mctc_io_convert, only : autoaa
   use mctc_io_math, only : matinv_3x3
   use mctc_io_symbols, only : to_symbol
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_genformat


contains


subroutine write_genformat(mol, unit)
   class(structure_type), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: iat, izp
   real(wp), parameter :: zero3(3) = 0.0_wp
   real(wp), allocatable :: inv_lat(:, :)
   real(wp), allocatable :: abc(:, :)

   write(unit, '(i0, 1x)', advance='no') mol%nat
   if (.not.any(mol%periodic)) then
      write(unit, '("C")') ! cluster
   else
      if (mol%info%cartesian) then
         write(unit, '("S")') ! supercell
      else
         write(unit, '("F")') ! fractional
      endif
   endif

   do izp = 1, mol%nid
      write(unit, '(1x, a)', advance='no') trim(mol%sym(izp))
   enddo
   write(unit, '(a)')

   if (.not.any(mol%periodic) .or. mol%info%cartesian) then
      ! now write the cartesian coordinates
      do iat = 1, mol%nat
         write(unit, '(2i5, 3es24.14)') iat, mol%id(iat), mol%xyz(:, iat)*autoaa
      enddo
   else
      inv_lat = matinv_3x3(mol%lattice)
      abc = matmul(inv_lat, mol%xyz)
      ! now write the fractional coordinates
      do iat = 1, mol%nat
         write(unit, '(2i5, 3es24.15)') iat, mol%id(iat), abc(:, iat)
      enddo
   endif

   if (any(mol%periodic)) then
      ! scaling factor for lattice parameters is always one
      write(unit, '(3f20.14)') zero3
      ! write the lattice parameters
      write(unit, '(3f20.14)') mol%lattice(:, :)*autoaa
   endif

end subroutine write_genformat


end module mctc_io_write_genformat

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_pdb
   use mctc_env_accuracy, only : wp
   use mctc_io_convert, only : autoaa
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_pdb


contains


subroutine write_pdb(mol, unit, number)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: unit
   integer, intent(in), optional :: number
   character(len=6) :: w1
   character(len=4) :: sym
   character(len=2) :: a_charge
   character(len=1) :: last_chain
   logical :: last_het
   integer :: offset, iat, jat
   real(wp) :: xyz(3)
   character(len=*), parameter :: pdb_format = &
      &  '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)'


   if (present(number)) write(unit, '("MODEL ",4x,i4)') number
   if (allocated(mol%pdb)) then
      offset = 0
      last_chain = mol%pdb(1)%chains
      last_het = mol%pdb(1)%het
      do iat = 1, mol%nat

         ! handle the terminator
         if (mol%pdb(iat)%het .neqv. last_het) then
            write(unit, '("TER   ",i5,6x,a3,1x,a1,i4)') iat + offset, &
               &  mol%pdb(iat-1)%residue, last_chain, mol%pdb(iat)%residue_number
            last_het = .not.last_het
            offset = offset+1
         else if (mol%pdb(iat)%chains /= last_chain) then
            write(unit, '("TER   ",i5,6x,a3,1x,a1,i4)') iat + offset, &
               &  mol%pdb(iat-1)%residue, last_chain, mol%pdb(iat)%residue_number
            offset = offset+1
         endif

         jat = iat + offset
         if (mol%pdb(iat)%het) then
            w1 = "HETATM"
         else
            w1 = "ATOM  "
         endif


         sym = adjustr(mol%sym(mol%id(iat))(1:2))
         xyz = mol%xyz(:,iat) * autoaa
         if (mol%pdb(iat)%charge < 0) then
            write(a_charge, '(i1,"-")') abs(mol%pdb(iat)%charge)
         else if (mol%pdb(iat)%charge > 0) then
            write(a_charge, '(i1,"+")') abs(mol%pdb(iat)%charge)
         else
            a_charge = '  '
         endif

         write(unit, pdb_format) &
            &  w1, jat, mol%pdb(iat)%name, mol%pdb(iat)%loc, &
            &  mol%pdb(iat)%residue, mol%pdb(iat)%chains, mol%pdb(iat)%residue_number, &
            &  mol%pdb(iat)%code, xyz, 1.0_wp, 0.0_wp, mol%pdb(iat)%segid, &
            &  sym, a_charge
      enddo
   else
      do iat = 1, mol%nat
         w1 = "HETATM"
         sym = adjustr(mol%sym(mol%id(iat))(1:2))
         xyz = mol%xyz(:,iat) * autoaa
         a_charge = '  '

         write(unit, pdb_format) &
            &  w1, iat, sym, " ", &
            &  "UNK", "A", 1, " ", xyz, 1.0_wp, 0.0_wp, "    ", &
            &  sym, "  "
      enddo
   end if

   if (present(number)) then
      write(unit, '("ENDMDL")')
   else
      write(unit, '("END")')
   endif

end subroutine write_pdb


end module mctc_io_write_pdb

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_turbomole
   use mctc_env_accuracy, only : wp
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_coord


contains


subroutine write_coord(mol, unit)
   class(structure_type), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: iat

   write(unit, '(a)') "$coord"
   do iat = 1, mol%nat
      write(unit, '(3es24.14, 6x, a)') mol%xyz(:, iat), trim(mol%sym(mol%id(iat)))
   enddo
   write(unit, '(a, 1x, i0)') "$periodic", count(mol%periodic)
   if (any(mol%periodic)) then
      write(unit, '(a)') "$lattice bohr"
      write(unit, '(3f20.14)') mol%lattice
   endif
   write(unit, '(a)') "$end"

end subroutine write_coord


end module mctc_io_write_turbomole

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_vasp
   use mctc_env_accuracy, only : wp
   use mctc_io_convert, only : autoaa
   use mctc_io_math, only : matinv_3x3
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_vasp


contains


subroutine write_vasp(self, unit, comment_line)
   class(structure_type), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: comment_line
   integer :: i, j, izp
   integer, allocatable :: kinds(:), species(:)
   real(wp), allocatable :: inv_lat(:, :)
   real(wp), allocatable :: abc(:, :)

   allocate(species(self%nat))
   allocate(kinds(self%nat), source=1)

   j = 0
   izp = 0
   do i = 1, self%nat
      if (izp.eq.self%id(i)) then
         kinds(j) = kinds(j)+1
      else
         j = j+1
         izp = self%id(i)
         species(j) = self%id(i)
      endif
   enddo

   ! use vasp 5.x format
   if (present(comment_line)) then
      write(unit, '(a)') comment_line
   else
      if (allocated(self%comment)) then
         write(unit, '(a)') self%comment
      else
         write(unit, '(a)')
      end if
   end if

   ! scaling factor for lattice parameters is always one
   write(unit, '(f20.14)') self%info%scale
   ! write the lattice parameters
   if (allocated(self%lattice)) then
      do i = 1, 3
         write(unit, '(3f20.14)') self%lattice(:, i)*autoaa/self%info%scale
      enddo
   else
      write(unit, '(3f20.14)') spread(0.0_wp, 1, 9)
   end if

   do i = 1, j
      write(unit, '(1x, a)', advance='no') self%sym(species(i))
   enddo
   write(unit, '(a)')

   ! write the count of the consequtive atom types
   do i = 1, j
      write(unit, '(1x, i0)', advance='no') kinds(i)
   enddo
   write(unit, '(a)')
   deallocate(kinds, species)

   if (self%info%selective) write(unit, '("Selective")')

   ! we write cartesian coordinates
   if (.not.allocated(self%lattice) .or. self%info%cartesian) then
      write(unit, '("Cartesian")')

      ! now write the cartesian coordinates
      do i = 1, self%nat
         write(unit, '(3f20.14)') self%xyz(:, i)*autoaa/self%info%scale
      enddo
   else
      write(unit, '("Direct")')
      inv_lat = matinv_3x3(self%lattice)
      abc = matmul(inv_lat, self%xyz)

      ! now write the fractional coordinates
      do i = 1, self%nat
         write(unit, '(3f20.14)') abc(:, i)
      enddo
   endif

end subroutine write_vasp


end module mctc_io_write_vasp

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write_xyz
   use mctc_io_convert, only : autoaa
   use mctc_io_structure, only : structure_type
   implicit none
   private

   public :: write_xyz


contains


subroutine write_xyz(self, unit, comment_line)
   class(structure_type), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: comment_line
   integer :: iat
   logical :: expo

   write(unit, '(i0)') self%nat
   if (present(comment_line)) then
      write(unit, '(a)') comment_line
   else
      if (allocated(self%comment)) then
         write(unit, '(a)') self%comment
      else
         write(unit, '(a)')
      end if
   end if
   expo = maxval(self%xyz) > 1.0e+5 .or. minval(self%xyz) < -1.0e+5
   if (expo) then
      do iat = 1, self%nat
         write(unit, '(a4, 1x, 3es24.14)') &
            & self%sym(self%id(iat)), self%xyz(:, iat)*autoaa
      enddo
   else
      do iat = 1, self%nat
         write(unit, '(a4, 1x, 3f24.14)') &
            & self%sym(self%id(iat)), self%xyz(:, iat)*autoaa
      enddo
   end if

end subroutine write_xyz


end module mctc_io_write_xyz

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_data
   use multicharge_data_covrad, only : get_covalent_rad
   implicit none
   public

end module multicharge_data

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_reference
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use dftd4_data, only : get_hardness, get_effective_charge
   implicit none
   private

   public :: get_nref, set_refcn, set_refgw, set_refq, set_refalpha

   interface get_nref
      module procedure :: get_nref_sym
      module procedure :: get_nref_num
   end interface get_nref

   interface set_refcn
      module procedure :: set_refcn_sym
      module procedure :: set_refcn_num
   end interface set_refcn

   interface set_refgw
      module procedure :: set_refgw_sym
      module procedure :: set_refgw_num
   end interface set_refgw

   interface set_refq
      module procedure :: set_refq_sym
      module procedure :: set_refq_num
   end interface set_refq

   interface set_refalpha
      module procedure :: set_refalpha_sym
      module procedure :: set_refalpha_num
   end interface set_refalpha

   integer, parameter :: max_elem = 118

   integer, dimension(max_elem)      :: refn ! for D4
   real(wp),dimension(7,max_elem)    :: refq
   real(wp),dimension(7,max_elem)    :: refh
   real(wp),dimension(7,max_elem)    :: dftq,pbcq,gffq,solq,clsq
   real(wp),dimension(7,max_elem)    :: dfth,pbch,gffh,solh,clsh
   real(wp),dimension(7,max_elem)    :: hcount
   real(wp),dimension(7,max_elem)    :: ascale
   real(wp),dimension(7,max_elem)    :: refcovcn
   real(wp),dimension(7,max_elem)    :: refcn
   integer, dimension(7,max_elem)    :: refsys
   real(wp),dimension(23,7,max_elem) :: alphaiw
   real(wp),dimension(17)       :: secq
   real(wp),dimension(17)       :: sscale
   real(wp),dimension(17)       :: seccn
   real(wp),dimension(17)       :: seccnd3
   real(wp),dimension(23,17)    :: secaiw

   include 'reference.inc'

contains


!> Get number of references for a given element symbol
elemental function get_nref_sym(sym) result(n)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Number of references
   integer :: n

   n = get_nref(to_number(sym))

end function get_nref_sym


!> Get number of references for a given atomic number
elemental function get_nref_num(num) result(n)

   !> Atomic number
   integer, intent(in) :: num

   !> Number of references
   integer :: n

   if (num > 0 .and. num <= size(refn)) then
      n = refn(num)
   else
      n = 0
   end if

end function get_nref_num


!> Set the reference coordination numbers for an element symbol
pure subroutine set_refcn_sym(cn, sym)

   !> Reference coordination number
   real(wp), intent(out) :: cn(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refcn(cn, to_number(sym))

end subroutine set_refcn_sym


!> Set the reference coordination numbers for an atomic number
pure subroutine set_refcn_num(cn, num)

   !> Reference coordination number
   real(wp), intent(out) :: cn(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref

   cn(:) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      cn(:ref) = refcovcn(:ref, num)
   end if

end subroutine set_refcn_num


!> Set the number of gaussian weights for an element symbol
pure subroutine set_refgw_sym(ngw, sym)

   !> Number of gaussian weights
   integer, intent(out) :: ngw(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refgw(ngw, to_number(sym))

end subroutine set_refgw_sym


!> Set the number of gaussian weights for an atomic number
pure subroutine set_refgw_num(ngw, num)

   !> Number of gaussian weights
   integer, intent(out) :: ngw(:)

   !> Atomic number
   integer, intent(in) :: num

   integer, parameter :: max_cn = 19
   integer :: icn, ir, ref
   integer :: cnc(0:max_cn)

   ngw(:) = 1
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      cnc(:) = [1, spread(0, 1, max_cn)]
      do ir = 1, ref
         icn = min(nint(refcn(ir, num)), max_cn)
         cnc(icn) = cnc(icn) + 1
      end do
      do ir = 1, ref
         icn = cnc(min(nint(refcn(ir, num)), max_cn))
         ngw(ir) = icn*(icn+1)/2
      end do
   end if

end subroutine set_refgw_num


!> Set the reference partial charges for an element symbol
pure subroutine set_refq_sym(q, sym)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refq(q, to_number(sym))

end subroutine set_refq_sym


!> Set the reference partial charges for an atomic number
pure subroutine set_refq_num(q, num)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref

   q(:) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      q(:ref) = clsq(:ref, num)
   end if

end subroutine set_refq_num


!> Set the reference polarizibility for an element symbol
pure subroutine set_refalpha_sym(alpha, ga, gc, sym)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refalpha(alpha, ga, gc, to_number(sym))

end subroutine set_refalpha_sym


!> Set the reference polarizibility for an atomic number
pure subroutine set_refalpha_num(alpha, ga, gc, num)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref
   integer :: ir, is
   real(wp) :: iz
   real(wp) :: aiw(23)

   alpha(:, :) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      do ir = 1, ref
         is = refsys(ir, num)
         iz = get_effective_charge(is)
         aiw = sscale(is)*secaiw(:, is) &
            &    * zeta(ga, get_hardness(is)*gc, iz, clsh(ir, num)+iz)
         alpha(:, ir) = max(ascale(ir, num)*(alphaiw(:, ir, num) &
            &            - hcount(ir, num)*aiw), 0.0_wp)
      end do
   end if

end subroutine set_refalpha_num


!> charge scaling function
elemental function zeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: zeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1.0_wp - exp( c * ( 1.0_wp - qref/qmod ) ) ) )
   endif

end function zeta


end module dftd4_reference

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_read
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_filetype, only : filetype, get_filetype
   use mctc_io_read_ctfile, only : read_molfile, read_sdf
   use mctc_io_read_gaussian, only : read_gaussian_external
   use mctc_io_read_genformat, only : read_genformat
   use mctc_io_read_pdb, only : read_pdb
   use mctc_io_read_turbomole, only : read_coord
   use mctc_io_read_vasp, only : read_vasp
   use mctc_io_read_xyz, only : read_xyz
   use mctc_io_structure, only : structure_type, new_structure
   implicit none
   private

   public :: read_structure
   public :: structure_reader, get_structure_reader


   interface read_structure
      module procedure :: read_structure_from_file
      module procedure :: read_structure_from_unit
   end interface read_structure


   abstract interface
      !> Read molecular structure data from formatted unit
      subroutine structure_reader(self, unit, error)
         import :: structure_type, error_type

         !> Instance of the molecular structure data
         type(structure_type), intent(out) :: self

         !> File handle
         integer, intent(in) :: unit

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

      end subroutine structure_reader
   end interface


contains


subroutine read_structure_from_file(self, file, error, format)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> Name of the file to read
   character(len=*), intent(in) :: file

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> File type format hint
   integer, intent(in), optional :: format

   logical :: exist
   integer :: unit, stat, ftype

   inquire(file=file, exist=exist)
   if (.not.exist) then
      call fatal_error(error, "File '"//file//"' cannot be found")
      return
   end if

   open(file=file, newunit=unit, status='old', iostat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot open '"//file//"'")
      return
   end if

   if (present(format)) then
      ftype = format
   else
      ftype = get_filetype(file)
   end if

   call read_structure(self, unit, ftype, error)
   close(unit)

end subroutine read_structure_from_file


subroutine read_structure_from_unit(self, unit, ftype, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> File type to read
   integer, intent(in) :: ftype

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   procedure(structure_reader), pointer :: reader

   call get_structure_reader(reader, ftype)
   if (.not.associated(reader)) then
      call fatal_error(error, "Cannot read structure from unknown file format")
      return
   end if

   call reader(self, unit, error)

end subroutine read_structure_from_unit


!> Retrieve reader for corresponding file type
subroutine get_structure_reader(reader, ftype)

   !> Reader for the specified file type
   procedure(structure_reader), pointer, intent(out) :: reader

   !> File type to read
   integer, intent(in) :: ftype

   nullify(reader)

   select case(ftype)
   case(filetype%xyz)
      reader => read_xyz

   case(filetype%molfile)
      reader => read_molfile

   case(filetype%pdb)
      reader => read_pdb

   case(filetype%gen)
      reader => read_genformat

   case(filetype%sdf)
      reader => read_sdf

   case(filetype%vasp)
      reader => read_vasp

   case(filetype%tmol)
      reader => read_coord

   case(filetype%gaussian)
      reader => read_gaussian_external

   end select

end subroutine get_structure_reader


end module mctc_io_read

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module mctc_io_write
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_filetype, only : filetype, get_filetype
   use mctc_io_write_ctfile, only : write_molfile, write_sdf
   use mctc_io_write_gaussian, only : write_gaussian_external
   use mctc_io_write_genformat, only : write_genformat
   use mctc_io_write_pdb, only : write_pdb
   use mctc_io_write_turbomole, only : write_coord
   use mctc_io_write_vasp, only : write_vasp
   use mctc_io_write_xyz, only : write_xyz
   use mctc_io_structure, only : structure_type, new_structure
   implicit none
   private

   public :: write_structure


   interface write_structure
      module procedure :: write_structure_to_file
      module procedure :: write_structure_to_unit
   end interface write_structure


contains


subroutine write_structure_to_file(self, file, error, format)

   !> Instance of the molecular structure data
   class(structure_type), intent(in) :: self

   !> Name of the file to read
   character(len=*), intent(in) :: file

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> File type format hint
   integer, intent(in), optional :: format

   integer :: unit, ftype, stat

   open(file=file, newunit=unit, iostat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot open '"//file//"'")
      return
   end if

   if (present(format)) then
      ftype = format
   else
      ftype = get_filetype(file)
   end if

   ! Unknown file type is inacceptable in this situation,
   ! try to figure something at least something out
   if (ftype == filetype%unknown) then
      if (any(self%periodic)) then
         ftype = filetype%vasp
      else if (allocated(self%sdf)) then
         ftype = filetype%sdf
      else if (allocated(self%pdb)) then
         ftype = filetype%pdb
      else
         ftype = filetype%xyz
      end if
   end if

   call write_structure(self, unit, ftype, error)
   close(unit)

end subroutine write_structure_to_file


subroutine write_structure_to_unit(self, unit, ftype, error)

   !> Instance of the molecular structure data
   class(structure_type), intent(in) :: self

   !> File handle
   integer, intent(in) :: unit

   !> File type to read
   integer, intent(in) :: ftype

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(ftype)
   case default
      call fatal_error(error, "Cannot write unknown file format")

   case(filetype%xyz)
      call write_xyz(self, unit)

   case(filetype%molfile)
      call write_molfile(self, unit)

   case(filetype%pdb)
      call write_pdb(self, unit)

   case(filetype%gen)
      call write_genformat(self, unit)

   case(filetype%sdf)
      call write_sdf(self, unit)

   case(filetype%vasp)
      call write_vasp(self, unit)

   case(filetype%tmol)
      call write_coord(self, unit)

   case(filetype%gaussian)
      call write_gaussian_external(self, unit)

   end select

end subroutine write_structure_to_unit


end module mctc_io_write

! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Input and output module of the tool chain library.
!>
!> This module exports the basic [[structure_type]] as well as routines
!> to read it from a file or formatted unit ([[read_structure]]) or write
!> it to a formatted unit ([[write_structure]]).
!>
!> Both [[read_structure]] and [[write_structure]] take format hints from
!> the filetype enumerator. File names can be translated to the respective
!> enumerator by using the [[get_filetype]] function. This can be useful in
!> case the caller routine wants to open the formatted unit itself or uses
!> a non-standard file extension.
module mctc_io
   use mctc_io_filetype, only : filetype, get_filetype
   use mctc_io_read, only : read_structure
   use mctc_io_structure, only : structure_type, new_structure, new
   use mctc_io_symbols, only : to_symbol, to_number
   use mctc_io_write, only : write_structure
   implicit none
   private

   public :: filetype, get_filetype
   public :: read_structure, write_structure
   public :: structure_type, new_structure, new
   public :: to_symbol, to_number


contains
end module mctc_io

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Generic interface to define damping functions for the DFT-D4 model
module dftd4_damping
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: damping_param, dispersion_interface


   type, abstract :: damping_param
   contains
      procedure(dispersion_interface), deferred :: get_dispersion2
      procedure(dispersion_interface), deferred :: get_dispersion3
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion2
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion3
   end type damping_param


   abstract interface
      !> Evaluation of the dispersion energy expression
      subroutine dispersion_interface(self, mol, trans, cutoff, r4r2, &
            & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

         !> Molecular structure data
         class(structure_type), intent(in) :: mol

         !> Lattice points
         real(wp), intent(in) :: trans(:, :)

         !> Real space cutoff
         real(wp), intent(in) :: cutoff

         !> Expectation values for r4 over r2 operator
         real(wp), intent(in) :: r4r2(:)

         !> C6 coefficients for all atom pairs.
         real(wp), intent(in) :: c6(:, :)

         !> Derivative of the C6 w.r.t. the coordination number
         real(wp), intent(in), optional :: dc6dcn(:, :)

         !> Derivative of the C6 w.r.t. the partial charges
         real(wp), intent(in), optional :: dc6dq(:, :)

         !> Dispersion energy
         real(wp), intent(inout) :: energy(:)

         !> Derivative of the energy w.r.t. the coordination number
         real(wp), intent(inout), optional :: dEdcn(:)

         !> Derivative of the energy w.r.t. the partial charges
         real(wp), intent(inout), optional :: dEdq(:)

         !> Dispersion gradient
         real(wp), intent(inout), optional :: gradient(:, :)

         !> Dispersion virial
         real(wp), intent(inout), optional :: sigma(:, :)
      end subroutine dispersion_interface

      !> Evaluation of the pairwise representation of the dispersion energy
      subroutine pairwise_dispersion_interface(self, mol, trans, cutoff, r4r2, c6, energy)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

         !> Molecular structure data
         class(structure_type), intent(in) :: mol

         !> Lattice points
         real(wp), intent(in) :: trans(:, :)

         !> Real space cutoff
         real(wp), intent(in) :: cutoff

         !> Expectation values for r4 over r2 operator
         real(wp), intent(in) :: r4r2(:)

         !> C6 coefficients for all atom pairs.
         real(wp), intent(in) :: c6(:, :)

         !> Pairwise representation of the dispersion energy
         real(wp), intent(inout) :: energy(:, :)
      end subroutine pairwise_dispersion_interface
   end interface


end module dftd4_damping

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Definition of the D4 dispersion model for the evaluation of C6 coefficients.
module dftd4_model
   use ieee_arithmetic, only : ieee_is_nan
   use dftd4_data, only : get_covalent_rad, get_r4r2_val, get_effective_charge, &
      get_electronegativity, get_hardness
   use dftd4_reference
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: d4_model, new_d4_model


   !> Base D4 dispersion model to evaluate C6 coefficients
   type :: d4_model

      !> Charge scaling height
      real(wp) :: ga

      !> Charge scaling steepness
      real(wp) :: gc

      !> Weighting factor for CN interpolation
      real(wp) :: wf

      !> Effective nuclear charges
      real(wp), allocatable :: zeff(:)

      !> Chemical hardness
      real(wp), allocatable :: eta(:)

      !> Electronegativity
      real(wp), allocatable :: en(:)

      !> Covalent radii for coordination number
      real(wp), allocatable :: rcov(:)

      !> Expectation values for C8 extrapolation
      real(wp), allocatable :: r4r2(:)

      !> Number of reference systems
      integer, allocatable :: ref(:)

      !> Number of Gaussian weights for each reference
      integer, allocatable :: ngw(:, :)

      !> Reference coordination numbers
      real(wp), allocatable :: cn(:, :)

      !> Reference partial charges
      real(wp), allocatable :: q(:, :)

      !> Reference dynamic polarizibilities
      real(wp), allocatable :: aiw(:, :, :)

      !> Reference C6 coefficients
      real(wp), allocatable :: c6(:, :, :, :)

   contains

      !> Generate weights for all reference systems
      procedure :: weight_references

      !> Evaluate C6 coefficient
      procedure :: get_atomic_c6

      !> Evaluate atomic polarizibilities
      procedure :: get_polarizibilities

   end type d4_model


   !> Default maximum charge scaling height for partial charge extrapolation
   real(wp), parameter :: ga_default = 3.0_wp

   !> Default charge scaling steepness for partial charge extrapolation
   real(wp), parameter :: gc_default = 2.0_wp

   !> Default weighting factor for coordination number interpolation
   real(wp), parameter :: wf_default = 6.0_wp


contains


!> Create new dispersion model from molecular structure input
subroutine new_d4_model(self, mol, ga, gc, wf)

   !> Instance of the dispersion model
   type(d4_model), intent(out) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Charge scaling height
   real(wp), intent(in), optional :: ga

   !> Charge scaling steepness
   real(wp), intent(in), optional :: gc

   !> Weighting factor for coordination number interpolation
   real(wp), intent(in), optional :: wf

   integer :: isp, izp, iref, jsp, jzp, jref
   integer :: mref
   real(wp) :: aiw(23), c6
   real(wp), parameter :: thopi = 3.0_wp/pi

   if (present(ga)) then
      self%ga = ga
   else
      self%ga = ga_default
   end if

   if (present(gc)) then
      self%gc = gc
   else
      self%gc = gc_default
   end if

   if (present(wf)) then
      self%wf = wf
   else
      self%wf = wf_default
   end if

   allocate(self%rcov(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%rcov(isp) = get_covalent_rad(izp)
   end do

   allocate(self%en(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%en(isp) = get_electronegativity(izp)
   end do

   allocate(self%zeff(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%zeff(isp) = get_effective_charge(izp)
   end do

   allocate(self%eta(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%eta(isp) = get_hardness(izp)
   end do

   allocate(self%r4r2(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%r4r2(isp) = get_r4r2_val(izp)
   end do

   allocate(self%ref(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%ref(isp) = get_nref(izp)
   end do

   mref = maxval(self%ref)
   allocate(self%cn(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refcn(self%cn(:, isp), izp)
   end do

   allocate(self%q(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refq(self%q(:, isp), izp)
   end do

   allocate(self%ngw(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refgw(self%ngw(:, isp), izp)
   end do

   allocate(self%aiw(23, mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refalpha(self%aiw(:, :, isp), self%ga, self%gc, izp)
   end do

   allocate(self%c6(mref, mref, mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         do iref = 1, self%ref(isp)
            do jref = 1, self%ref(jsp)
               aiw(:) = self%aiw(:, iref, isp) * self%aiw(:, jref, jsp)
               c6 = thopi * trapzd(aiw)
               self%c6(jref, iref, jsp, isp) = c6
               self%c6(iref, jref, isp, jsp) = c6
            end do
         end do
      end do
   end do

end subroutine new_d4_model


!> Calculate the weights of the reference system and the derivatives w.r.t.
!> coordination number for later use.
subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)

   !> Instance of the dispersion model
   class(d4_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Coordination number of every atom
   real(wp), intent(in) :: cn(:)

   !> Partial charge of every atom
   real(wp), intent(in) :: q(:)

   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)

   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out), optional :: gwdcn(:, :)

   !> derivative of the weighting function w.r.t. the charge scaling
   real(wp), intent(out), optional :: gwdq(:, :)

   integer :: iat, izp, iref, igw
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi

   if (present(gwdcn) .and. present(gwdq)) then
      gwvec(:, :) = 0.0_wp
      gwdcn(:, :) = 0.0_wp
      gwdq(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, gwdcn, gwdq, mol, self, cn, q) private(iat, izp, iref, &
      !$omp& igw, norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc
         norm = 0.0_wp
         dnorm = 0.0_wp
         do iref = 1, self%ref(izp)
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
               norm = norm + gw
               dnorm = dnorm + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
            end do
         end do
         norm = 1.0_wp / norm
         do iref = 1, self%ref(izp)
            expw = 0.0_wp
            expd = 0.0_wp
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
               expw = expw + gw
               expd = expd + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
            end do
            gwk = expw * norm
            if (ieee_is_nan(gwk)) then
               if (maxval(self%cn(:self%ref(izp), izp)) == self%cn(iref, izp)) then
                  gwk = 1.0_wp
               else
                  gwk = 0.0_wp
               end if
            end if
            gwvec(iref, iat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
            gwdq(iref, iat) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)

            dgwk = norm * (expd - expw * dnorm * norm)
            if (ieee_is_nan(dgwk)) then
               dgwk = 0.0_wp
            end if
            gwdcn(iref, iat) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
         end do
      end do

   else

      gwvec(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, mol, self, cn, q) &
      !$omp private(iat, izp, iref, igw, norm, gw, expw, gwk, wf, zi, gi)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc
         norm = 0.0_wp
         do iref = 1, self%ref(izp)
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
               norm = norm + gw
            end do
         end do
         norm = 1.0_wp / norm
         do iref = 1, self%ref(izp)
            expw = 0.0_wp
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               expw = expw + weight_cn(wf, cn(iat), self%cn(iref, izp))
            end do
            gwk = expw * norm
            if (ieee_is_nan(gwk)) then
               if (maxval(self%cn(:self%ref(izp), izp)) == self%cn(iref, izp)) then
                  gwk = 1.0_wp
               else
                  gwk = 0.0_wp
               end if
            end if
            gwvec(iref, iat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
         end do
      end do
   end if

end subroutine weight_references


!> Calculate atomic dispersion coefficients and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   !> Instance of the dispersion model
   class(d4_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> Derivative of the weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charge
   real(wp), intent(out), optional :: dc6dq(:, :)

   integer :: iat, jat, izp, jzp, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj

   if (present(gwdcn).and.present(dc6dcn) &
      & .and.present(gwdq).and.present(dc6dq)) then
      c6(:, :) = 0.0_wp
      dc6dcn(:, :) = 0.0_wp
      dc6dq(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(c6, dc6dcn, dc6dq, mol, self, gwvec, gwdcn, gwdq) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6, dc6dqi, dc6dqj, &
      !$omp& dc6dcni, dc6dcnj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            dc6dcni = 0.0_wp
            dc6dcnj = 0.0_wp
            dc6dqi = 0.0_wp
            dc6dqj = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
                  dc6dqi = dc6dqi + gwdq(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dqj = dc6dqj + gwvec(iref, iat) * gwdq(jref, jat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
            dc6dcn(iat, jat) = dc6dcni
            dc6dcn(jat, iat) = dc6dcnj
            dc6dq(iat, jat) = dc6dqi
            dc6dq(jat, iat) = dc6dqj
         end do
      end do

   else

      c6(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(c6, mol, self, gwvec) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
         end do
      end do
   end if

end subroutine get_atomic_c6


!> Calculate atomic polarizibilities and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_polarizibilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dadq)

   !> Instance of the dispersion model
   class(d4_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> Derivative of the weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :)

   !> Static polarizibilities for all atoms.
   real(wp), intent(out) :: alpha(:)

   !> Derivative of the polarizibility w.r.t. the coordination number
   real(wp), intent(out), optional :: dadcn(:)

   !> Derivative of the polarizibility w.r.t. the partial charge
   real(wp), intent(out), optional :: dadq(:)

   integer :: iat, izp, iref
   real(wp) :: refa, da, dadcni, dadqi

   if (present(gwdcn).and.present(dadcn) &
      & .and.present(gwdq).and.present(dadq)) then
      alpha(:) = 0.0_wp
      dadcn(:) = 0.0_wp
      dadq(:) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(alpha, dadcn, dadq, mol, self, gwvec, gwdcn, gwdq) &
      !$omp private(iat, izp, iref, refa, da, dadqi, dadcni)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         da = 0.0_wp
         dadcni = 0.0_wp
         dadqi = 0.0_wp
         do iref = 1, self%ref(izp)
            refa = self%aiw(1, iref, izp)
            da = da + gwvec(iref, iat) * refa
            dadcni = dadcni + gwdcn(iref, iat) * refa
            dadqi = dadqi + gwdq(iref, iat) * refa
         end do
         alpha(iat) = da
         dadcn(iat) = dadcni
         dadq(iat) = dadqi
      end do

   else

      alpha(:) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(alpha, mol, self, gwvec) private(iat, izp, iref, refa, da)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         da = 0.0_wp
         do iref = 1, self%ref(izp)
            da = da + gwvec(iref, iat) * self%aiw(1, iref, izp)
         end do
         alpha(iat) = da
      end do
   end if

end subroutine get_polarizibilities


elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   real(wp) :: cngw
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

!> charge scaling function
elemental function zeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: zeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1.0_wp - exp( c * ( 1.0_wp - qref/qmod ) ) ) )
   endif

end function zeta

!> derivative of charge scaling function w.r.t. charge
elemental function dzeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: dzeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      dzeta = 0.0_wp
   else
      dzeta = - a * c * exp( c * ( 1.0_wp - qref/qmod ) ) &
         & * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

!> numerical Casimir--Polder integration
pure function trapzd(pol)
   real(wp), intent(in) :: pol(23)
   real(wp) :: trapzd

   real(wp), parameter :: freq(23) = [ &
      & 0.000001_wp, 0.050000_wp, 0.100000_wp, &
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.500000_wp, 0.600000_wp, 0.700000_wp, &
      & 0.800000_wp, 0.900000_wp, 1.000000_wp, &
      & 1.200000_wp, 1.400000_wp, 1.600000_wp, &
      & 1.800000_wp, 2.000000_wp, 2.500000_wp, &
      & 3.000000_wp, 4.000000_wp, 5.000000_wp, &
      & 7.500000_wp, 10.00000_wp]
   real(wp), parameter :: weights(23) = 0.5_wp * [ &
      &  ( freq (2) - freq (1) ),  &
      &  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
      &  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
      &  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
      &  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
      &  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
      &  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
      &  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
      &  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
      &  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
      &  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
      &  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
      &  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
      &  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
      &  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
      &  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
      &  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
      &  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
      &  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
      &  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
      &  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
      &  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
      &  ( freq(23) - freq(22) ) ]

   trapzd = sum(pol*weights)

end function trapzd


end module dftd4_model

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: get_coordination_number


   !> Steepness of counting function
   real(wp), parameter :: kcn = 7.5_wp

   !> Parameter for electronegativity scaling
   real(wp), parameter :: k4 = 4.10451_wp

   !> Parameter for electronegativity scaling
   real(wp), parameter :: k5 = 19.08857_wp

   !> Parameter for electronegativity scaling
   real(wp), parameter :: k6 = 2*11.28174_wp**2


contains


!> Geometric fractional coordination number, supports error function counting.
subroutine get_coordination_number(mol, trans, cutoff, rcov, en, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   if (present(dcndr) .and. present(dcndL)) then
      call ncoord_derf(mol, trans, cutoff, rcov, en, cn, dcndr, dcndL)
   else
      call ncoord_erf(mol, trans, cutoff, rcov, en, cn)
   end if

end subroutine get_coordination_number


subroutine ncoord_erf(mol, trans, cutoff, rcov, en, cn)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, cutoff2, den

   cn(:) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) reduction(+:cn) &
   !$omp shared(mol, trans, cutoff2, rcov, en) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf, den)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         den = k4*exp(-(abs(en(izp)-en(jzp)) + k5)**2/k6)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = den*erf_count(kcn, r1, rc)

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

         end do
      end do
   end do

end subroutine ncoord_erf


subroutine ncoord_derf(mol, trans, cutoff, rcov, en, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), sigma(3, 3), cutoff2, den

   cn(:) = 0.0_wp
   dcndr(:, :, :) = 0.0_wp
   dcndL(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) reduction(+:cn, dcndr, dcndL) &
   !$omp shared(mol, trans, cutoff2, rcov, en) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf, countd, sigma, den)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         den = k4*exp(-(abs(en(izp)-en(jzp)) + k5)**2/k6)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = den*erf_count(kcn, r1, rc)
            countd = den*derf_count(kcn, r1, rc) * rij/r1

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

            dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
            dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
            dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
            dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

            sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

            dcndL(:, :, iat) = dcndL(:, :, iat) + sigma
            if (iat /= jat) then
               dcndL(:, :, jat) = dcndL(:, :, jat) + sigma
            end if

         end do
      end do
   end do

end subroutine ncoord_derf


!> Error function counting function for coordination number contributions.
pure function erf_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))

end function erf_count


!> Derivative of the counting function w.r.t. the distance.
pure function derf_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp), parameter :: sqrtpi = sqrt(pi)

   real(wp) :: count

   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

end function derf_count


end module dftd4_ncoord

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Implementation of the Axilrod-Teller-Muto triple dipole dispersion
!> contribution with a modified zero (Chai--Head-Gordon) damping together
!> with the critical radii from the rational (Becke--Johnson) damping.
module dftd4_damping_atm
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: get_atm_dispersion


contains


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion(mol, trans, cutoff, s9, a1, a2, alp, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling parameter for critical radius
   real(wp), intent(in) :: a1

   !> Offset parameter for critical radius
   real(wp), intent(in) :: a2

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charges
   real(wp), intent(in), optional :: dc6dq(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad
   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
   real(wp) :: cutoff2, c9, dE, dGij(3), dGjk(3), dGik(3), dS(3, 3)

   if (abs(s9) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(dc6dq) &
      & .and. present(dEdq) .and. present(gradient) .and. present(sigma)
   cutoff2 = cutoff*cutoff

   if (grad) then
      !$omp parallel do schedule(runtime) default(none) &
      !$omp reduction(+:energy, gradient, sigma, dEdcn, dEdq) &
      !$omp shared(mol, trans, c6, s9, a1, a2, alp, r4r2, cutoff2, dc6dcn, dc6dq) &
      !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
      !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
      !$omp& r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang, c9, dE, dGij, dGjk, &
      !$omp& dGik, dS)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            c6ij = c6(jat, iat)
            r0ij = a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + a2
            do jtr = 1, size(trans, 2)
               vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
               r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
               if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
               do kat = 1, jat
                  kzp = mol%id(kat)
                  c6ik = c6(kat, iat)
                  c6jk = c6(kat, jat)
                  c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
                  r0ik = a1 * sqrt(3*r4r2(kzp)*r4r2(izp)) + a2
                  r0jk = a1 * sqrt(3*r4r2(kzp)*r4r2(jzp)) + a2
                  r0 = r0ij * r0ik * r0jk
                  triple = triple_scale(iat, jat, kat)
                  do ktr = 1, size(trans, 2)
                     vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                     r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                     if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                     vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                        & - trans(:, jtr)
                     r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                     if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                     r2 = r2ij*r2ik*r2jk
                     r1 = sqrt(r2)
                     r3 = r2 * r1
                     r5 = r3 * r2

                     fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                     ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                        & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                     rr = ang*fdmp

                     dfdmp = -2.0_wp * alp * (r0 / r1)**(alp / 3.0_wp) * fdmp**2

                     ! d/drij
                     dang = -0.375_wp * (r2ij**3 + r2ij**2 * (r2jk + r2ik)&
                        & + r2ij * (3.0_wp * r2jk**2 + 2.0_wp * r2jk*r2ik&
                        & + 3.0_wp * r2ik**2)&
                        & - 5.0_wp * (r2jk - r2ik)**2 * (r2jk + r2ik)) / r5
                     dGij(:) = c9 * (-dang*fdmp + ang*dfdmp) / r2ij * vij

                     ! d/drik
                     dang = -0.375_wp * (r2ik**3 + r2ik**2 * (r2jk + r2ij)&
                        & + r2ik * (3.0_wp * r2jk**2 + 2.0_wp * r2jk * r2ij&
                        & + 3.0_wp * r2ij**2)&
                        & - 5.0_wp * (r2jk - r2ij)**2 * (r2jk + r2ij)) / r5
                     dGik(:) = c9 * (-dang * fdmp + ang * dfdmp) / r2ik * vik

                     ! d/drjk
                     dang = -0.375_wp * (r2jk**3 + r2jk**2*(r2ik + r2ij)&
                        & + r2jk * (3.0_wp * r2ik**2 + 2.0_wp * r2ik * r2ij&
                        & + 3.0_wp * r2ij**2)&
                        & - 5.0_wp * (r2ik - r2ij)**2 * (r2ik + r2ij)) / r5
                     dGjk(:) = c9 * (-dang * fdmp + ang * dfdmp) / r2jk * vjk

                     dE = rr * c9 * triple
                     energy(iat) = energy(iat) - dE/3
                     energy(jat) = energy(jat) - dE/3
                     energy(kat) = energy(kat) - dE/3

                     gradient(:, iat) = gradient(:, iat) - dGij - dGik
                     gradient(:, jat) = gradient(:, jat) + dGij - dGjk
                     gradient(:, kat) = gradient(:, kat) + dGik + dGjk

                     dS(:, :) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                        & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                        & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)

                     sigma(:, :) = sigma + dS * triple

                     dEdcn(iat) = dEdcn(iat) - dE * 0.5_wp &
                        & * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik)
                     dEdcn(jat) = dEdcn(jat) - dE * 0.5_wp &
                        & * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk)
                     dEdcn(kat) = dEdcn(kat) - dE * 0.5_wp &
                        & * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk)

                     dEdq(iat) = dEdq(iat) - dE * 0.5_wp &
                        & * (dc6dq(iat, jat) / c6ij + dc6dq(iat, kat) / c6ik)
                     dEdq(jat) = dEdq(jat) - dE * 0.5_wp &
                        & * (dc6dq(jat, iat) / c6ij + dc6dq(jat, kat) / c6jk)
                     dEdq(kat) = dEdq(kat) - dE * 0.5_wp &
                        & * (dc6dq(kat, iat) / c6ik + dc6dq(kat, jat) / c6jk)
                  end do
               end do
            end do
         end do
      end do
   else
      !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
      !$omp shared(mol, trans, c6, s9, a1, a2, alp, r4r2, cutoff2) &
      !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
      !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
      !$omp& r1, r2, r3, r5, rr, fdmp, ang, c9, dE)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            c6ij = c6(jat, iat)
            r0ij = a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + a2
            do jtr = 1, size(trans, 2)
               vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
               r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
               if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
               do kat = 1, jat
                  kzp = mol%id(kat)
                  c6ik = c6(kat, iat)
                  c6jk = c6(kat, jat)
                  c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
                  r0ik = a1 * sqrt(3*r4r2(kzp)*r4r2(izp)) + a2
                  r0jk = a1 * sqrt(3*r4r2(kzp)*r4r2(jzp)) + a2
                  r0 = r0ij * r0ik * r0jk
                  triple = triple_scale(iat, jat, kat)
                  do ktr = 1, size(trans, 2)
                     vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                     r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                     if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                     vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                        & - trans(:, jtr)
                     r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                     if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                     r2 = r2ij*r2ik*r2jk
                     r1 = sqrt(r2)
                     r3 = r2 * r1
                     r5 = r3 * r2

                     fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                     ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                        & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                     rr = ang*fdmp

                     dE = rr * c9 * triple
                     energy(iat) = energy(iat) - dE/3
                     energy(jat) = energy(jat) - dE/3
                     energy(kat) = energy(kat) - dE/3
                  end do
               end do
            end do
         end do
      end do
   end if

end subroutine get_atm_dispersion


!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(triple)

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: triple

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         triple = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         triple = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         triple = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         triple = 0.5_wp
      end if
   end if

end function triple_scale


end module dftd4_damping_atm

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: get_coordination_number, cut_coordination_number


   !> Steepness of counting function
   real(wp), parameter :: kcn = 7.5_wp


contains


!> Geometric fractional coordination number, supports error function counting.
subroutine get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   !> Cut coordination number
   real(wp), intent(in), optional :: cut

   if (present(dcndr) .and. present(dcndL)) then
      call ncoord_derf(mol, trans, cutoff, rcov, cn, dcndr, dcndL)
   else
      call ncoord_erf(mol, trans, cutoff, rcov, cn)
   end if

   if (present(cut)) then
      call cut_coordination_number(cut, cn, dcndr, dcndL)
   end if

end subroutine get_coordination_number


subroutine ncoord_erf(mol, trans, cutoff, rcov, cn)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, cutoff2

   cn(:) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) schedule(runtime) reduction(+:cn) &
   !$omp shared(mol, trans, cutoff2, rcov) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = erf_count(kcn, r1, rc)

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

         end do
      end do
   end do

end subroutine ncoord_erf


subroutine ncoord_derf(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), sigma(3, 3), cutoff2

   cn(:) = 0.0_wp
   dcndr(:, :, :) = 0.0_wp
   dcndL(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:cn, dcndr, dcndL) shared(mol, trans, cutoff2, rcov) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf, countd, sigma)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = erf_count(kcn, r1, rc)
            countd = derf_count(kcn, r1, rc) * rij/r1

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

            dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
            dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
            dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
            dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

            sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

            dcndL(:, :, iat) = dcndL(:, :, iat) + sigma
            if (iat /= jat) then
               dcndL(:, :, jat) = dcndL(:, :, jat) + sigma
            end if

         end do
      end do
   end do

end subroutine ncoord_derf


!> Error function counting function for coordination number contributions.
elemental function erf_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))

end function erf_count


!> Derivative of the counting function w.r.t. the distance.
elemental function derf_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp), parameter :: sqrtpi = sqrt(pi)

   real(wp) :: count

   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

end function derf_count


!> Cutoff function for large coordination numbers
pure subroutine cut_coordination_number(cn_max, cn, dcndr, dcndL)

   !> Maximum CN (not strictly obeyed)
   real(wp), intent(in) :: cn_max

   !> On input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(:)

   !> On input derivative of CN w.r.t. cartesian coordinates,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(:, :, :)

   !> On input derivative of CN w.r.t. strain deformation,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndL(:, :, :)

   real(wp) :: dcnpdcn
   integer  :: iat

   if (present(dcndL)) then
      do iat = 1, size(cn)
         dcnpdcn = dlog_cn_cut(cn(iat), cn_max)
         dcndL(:, :, iat) = dcnpdcn*dcndL(:, :, iat)
      enddo
   endif

   if (present(dcndr)) then
      do iat = 1, size(cn)
         dcnpdcn = dlog_cn_cut(cn(iat), cn_max)
         dcndr(:, :, iat) = dcnpdcn*dcndr(:, :, iat)
      enddo
   endif

   do iat = 1, size(cn)
      cn(iat) = log_cn_cut(cn(iat), cn_max)
   enddo

end subroutine cut_coordination_number

elemental function log_cn_cut(cn, cnmax) result(cnp)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: cnp
   cnp = log(1.0_wp + exp(cnmax)) - log(1.0_wp + exp(cnmax - cn))
end function log_cn_cut

elemental function dlog_cn_cut(cn, cnmax) result(dcnpdcn)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: dcnpdcn
   dcnpdcn = exp(cnmax)/(exp(cnmax) + exp(cn))
end function dlog_cn_cut


end module multicharge_ncoord

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_wignerseitz
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge_cutoff, only : get_lattice_points
   implicit none
   private

   public :: wignerseitz_cell_type, new_wignerseitz_cell

   type :: wignerseitz_cell_type
      integer, allocatable :: nimg(:, :)
      integer, allocatable :: tridx(:, :, :)
      real(wp), allocatable :: trans(:, :)
   end type wignerseitz_cell_type


   !> Small cutoff threshold to create only closest cells
   real(wp), parameter :: thr = sqrt(epsilon(0.0_wp))

   !> Tolerance to consider equivalent images
   real(wp), parameter :: tol = 0.01_wp


contains


subroutine new_wignerseitz_cell(self, mol)

   !> Wigner-Seitz cell instance
   type(wignerseitz_cell_type), intent(out) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: iat, jat, ntr, nimg
   integer, allocatable :: tridx(:)
   real(wp) :: vec(3)
   real(wp), allocatable :: trans(:, :)

   call get_lattice_points(mol%periodic, mol%lattice, thr, trans)
   ntr = size(trans, 2)
   allocate(self%nimg(mol%nat, mol%nat), self%tridx(ntr, mol%nat, mol%nat), &
      & tridx(ntr))

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(mol, trans, self) private(iat, jat, vec, nimg, tridx)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_pairs(nimg, trans, vec, tridx)
         self%nimg(jat, iat) = nimg
         self%tridx(:, jat, iat) = tridx
      end do
   end do

   call move_alloc(trans, self%trans)
   
end subroutine new_wignerseitz_cell


subroutine get_pairs(iws, trans, rij, list)
   integer, intent(out) :: iws
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: trans(:, :)
   integer, intent(out) :: list(:)

   logical :: mask(size(list))
   real(wp) :: dist(size(list)), vec(3), r2
   integer :: itr, img, pos

   iws = 0
   img = 0
   list(:) = 0
   mask(:) = .true.

   do itr = 1, size(trans, 2)
      vec(:) = rij - trans(:, itr)
      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      if (r2 < thr) cycle
      img = img + 1
      dist(img) = r2
   end do

   if (img == 0) return

   pos = minloc(dist(:img), dim=1)

   r2 = dist(pos)
   mask(pos) = .false.

   iws = 1
   list(iws) = pos
   if (img <= iws) return

   do
      pos = minloc(dist(:img), dim=1, mask=mask(:img))
      if (abs(dist(pos) - r2) > tol) exit
      mask(pos) = .false.
      iws = iws + 1
      list(iws) = pos
   end do

end subroutine get_pairs


end module multicharge_wignerseitz

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Implementation of the rational (Becke--Johnson) damping function.
module dftd4_damping_rational
   use dftd4_damping, only : damping_param
   use dftd4_damping_atm, only : get_atm_dispersion
   use dftd4_data, only : get_r4r2_val
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: rational_damping_param


   !> Rational (Becke-Johnson) damping model
   type, extends(damping_param) :: rational_damping_param
      real(wp) :: s6 = 1.0_wp
      real(wp) :: s8
      real(wp) :: s9 = 1.0_wp
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: alp = 16.0_wp
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

      !> Evaluate pairwise representation of additive dispersion energy
      procedure :: get_pairwise_dispersion2

      !> Evaluate pairwise representation of non-additive dispersion energy
      procedure :: get_pairwise_dispersion3

   end type rational_damping_param


contains


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, r4r2, c6, dc6dcn, dc6dq, &
      & energy, dEdcn, dEdq, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charges
   real(wp), intent(in), optional :: dc6dq(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad
   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, d6, d8, edisp, gdisp
   real(wp) :: dE, dG(3), dS(3, 3)

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(dc6dq) &
      & .and. present(dEdq) .and. present(gradient) .and. present(sigma)
   cutoff2 = cutoff*cutoff

   if (grad) then
      !$omp parallel do schedule(runtime) default(none) &
      !$omp reduction(+:energy, gradient, sigma, dEdcn, dEdq) &
      !$omp shared(mol, self, c6, dc6dcn, dc6dq, trans, cutoff2, r4r2) &
      !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, t8, &
      !$omp& d6, d8, edisp, gdisp, dE, dG, dS)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rrij = 3*r4r2(izp)*r4r2(jzp)
            r0ij = self%a1 * sqrt(rrij) + self%a2
            c6ij = c6(jat, iat)
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

               t6 = 1.0_wp/(r2**3 + r0ij**6)
               t8 = 1.0_wp/(r2**4 + r0ij**8)

               d6 = -6*r2**2*t6**2
               d8 = -8*r2**3*t8**2

               edisp = self%s6*t6 + self%s8*rrij*t8
               gdisp = self%s6*d6 + self%s8*rrij*d8

               dE = -c6ij*edisp * 0.5_wp
               dG(:) = -c6ij*gdisp*vec
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

               energy(iat) = energy(iat) + dE
               dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * edisp
               dEdq(iat) = dEdq(iat) - dc6dq(iat, jat) * edisp
               sigma(:, :) = sigma + dS
               if (iat /= jat) then
                  energy(jat) = energy(jat) + dE
                  dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * edisp
                  dEdq(jat) = dEdq(jat) - dc6dq(jat, iat) * edisp
                  gradient(:, iat) = gradient(:, iat) + dG
                  gradient(:, jat) = gradient(:, jat) - dG
                  sigma(:, :) = sigma + dS
               end if
            end do
         end do
      end do
   else
      !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
      !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
      !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, &
      !$omp& t6, t8, edisp, dE)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rrij = 3*r4r2(izp)*r4r2(jzp)
            r0ij = self%a1 * sqrt(rrij) + self%a2
            c6ij = c6(jat, iat)
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

               t6 = 1.0_wp/(r2**3 + r0ij**6)
               t8 = 1.0_wp/(r2**4 + r0ij**8)

               edisp = self%s6*t6 + self%s8*rrij*t8

               dE = -c6ij*edisp * 0.5_wp

               energy(iat) = energy(iat) + dE
               if (iat /= jat) then
                  energy(jat) = energy(jat) + dE
               end if
            end do
         end do
      end do
   end if

end subroutine get_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_dispersion3(self, mol, trans, cutoff, r4r2, c6, dc6dcn, dc6dq, &
      & energy, dEdcn, dEdq, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charges
   real(wp), intent(in), optional :: dc6dq(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   call get_atm_dispersion(mol, trans, cutoff, self%s9, self%a1, self%a2, &
      & self%alp, r4r2, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, &
      & gradient, sigma)

end subroutine get_dispersion3


!> Evaluation of the dispersion energy expression projected on atomic pairs
subroutine get_pairwise_dispersion2(self, mol, trans, cutoff, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, edisp, dE

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, &
   !$omp& t6, t8, edisp, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = self%s6*t6 + self%s8*rrij*t8

            dE = -c6ij*edisp * 0.5_wp

            energy(jat, iat) = energy(jat, iat) + dE
            if (iat /= jat) then
               energy(iat, jat) = energy(iat, jat) + dE
            end if
         end do
      end do
   end do

end subroutine get_pairwise_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_pairwise_dispersion3(self, mol, trans, cutoff, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, ang
   real(wp) :: cutoff2, c9, dE

   if (abs(self%s9) < epsilon(1.0_wp)) return
   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, trans, c6, r4r2, cutoff2, self) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
   !$omp& r1, r2, r3, r5, rr, fdmp, ang, c9, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = self%a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + self%a2
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               c6ik = c6(kat, iat)
               c6jk = c6(kat, jat)
               c9 = -self%s9 * sqrt(abs(c6ij*c6ik*c6jk))
               r0ik = self%a1 * sqrt(3*r4r2(kzp)*r4r2(izp)) + self%a2
               r0jk = self%a1 * sqrt(3*r4r2(kzp)*r4r2(jzp)) + self%a2
               r0 = r0ij * r0ik * r0jk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(self%alp / 3.0_wp))
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dE = rr * c9 * triple/6
                  energy(jat, iat) = energy(jat, iat) - dE
                  energy(kat, iat) = energy(kat, iat) - dE
                  energy(iat, jat) = energy(iat, jat) - dE
                  energy(kat, jat) = energy(kat, jat) - dE
                  energy(iat, kat) = energy(iat, kat) - dE
                  energy(jat, kat) = energy(jat, kat) - dE
               end do
            end do
         end do
      end do
   end do

end subroutine get_pairwise_dispersion3


!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(triple)

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: triple

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         triple = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         triple = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         triple = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         triple = 0.5_wp
      end if
   end if

end function triple_scale


end module dftd4_damping_rational

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_model
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_ewald, only : get_alpha
   use multicharge_lapack, only : sytrf, sytrs, sytri
   use multicharge_wignerseitz, only : wignerseitz_cell_type, new_wignerseitz_cell
   implicit none
   private

   public :: mchrg_model_type, new_mchrg_model


   type :: mchrg_model_type
      real(wp), allocatable :: rad(:)
      real(wp), allocatable :: chi(:)
      real(wp), allocatable :: eta(:)
      real(wp), allocatable :: kcn(:)
   contains
      procedure :: solve
   end type mchrg_model_type


   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains


subroutine new_mchrg_model(self, chi, rad, eta, kcn)
   type(mchrg_model_type), intent(out) :: self
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: chi(:)
   real(wp), intent(in) :: eta(:)
   real(wp), intent(in) :: kcn(:)

   self%rad = rad
   self%chi = chi
   self%eta = eta
   self%kcn = kcn

end subroutine new_mchrg_model


subroutine get_vrhs(self, mol, cn, xvec, dxdcn)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(out) :: xvec(:)
   real(wp), intent(out), optional :: dxdcn(:)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   if (present(dxdcn)) then
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, xvec, dxdcn) private(iat, izp, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
         dxdcn(iat) = 0.5_wp*tmp
      end do
      dxdcn(mol%nat+1) = 0.0_wp
   else
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, xvec) private(iat, izp, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
      end do
   end if
   xvec(mol%nat+1) = mol%charge

end subroutine get_vrhs


subroutine get_dir_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans


subroutine get_amat_0d(self, mol, amat)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2*gam))/sqrt(r2)
         amat(jat, iat) = amat(jat, iat) + tmp
         amat(iat, jat) = amat(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat(iat, iat) = amat(iat, iat) + tmp
   end do

   amat(mol%nat+1, 1:mol%nat+1) = 1.0_wp
   amat(1:mol%nat+1, mol%nat+1) = 1.0_wp
   amat(mol%nat+1, mol%nat+1) = 0.0_wp

end subroutine get_amat_0d

subroutine get_amat_3d(self, mol, wsc, alpha, amat)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   amat(:, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, self, wsc, dtrans, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, gam, wsw, vec, dtmp, rtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat(jat, iat) = amat(jat, iat) + (dtmp + rtmp) * wsw
            amat(iat, jat) = amat(iat, jat) + (dtmp + rtmp) * wsw
         end do
      end do

      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
         amat(iat, iat) = amat(iat, iat) + (dtmp + rtmp) * wsw
      end do

      dtmp = self%eta(izp) + sqrt2pi / self%rad(izp) - 2 * alpha / sqrtpi
      amat(iat, iat) = amat(iat, iat) + dtmp
   end do

   amat(mol%nat+1, 1:mol%nat+1) = 1.0_wp
   amat(1:mol%nat+1, mol%nat+1) = 1.0_wp
   amat(mol%nat+1, mol%nat+1) = 0.0_wp

end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = erf(gam*r1)/r1 - erf(alp*r1)/r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, tmp

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      tmp = cos(dot_product(rij, vec)) * fac * exp(-0.25_wp*g2/(alp*alp))/g2
      amat = amat + tmp
   end do

end subroutine get_amat_rec_3d

subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) shared(mol, self, qvec) &
   !$omp private(iat, izp, jat, jzp, gam, r2, vec, dG, dS, dtmp, arg)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp/sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         arg = gam*gam*r2
         dtmp = 2.0_wp*gam*exp(-arg)/(sqrtpi*r2)-erf(sqrt(arg))/(r2*sqrt(r2))
         dG = dtmp*vec
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)
         atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
         atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat)
         dadr(:, jat, iat) = -dG*qvec(jat)
         dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
         dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
      end do
   end do

end subroutine get_damat_0d

subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vol, gam, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) &
   !$omp shared(mol, self, wsc, alpha, vol, dtrans, rtrans, qvec) &
   !$omp private(iat, izp, jat, jzp, img, gam, wsw, vec, dG, dS, &
   !$omp& dGr, dSr, dGd, dSd)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
         atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = -dG*qvec(jat) + dadr(:, jat, iat)
         dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
         dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
      end do

      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
         dS = dS + (dSd + dSr) * wsw
      end do
      dadL(:, :, iat) = +dS*qvec(iat) + dadL(:, :, iat)
   end do

end subroutine get_damat_3d

subroutine get_damat_dir_3d(rij, gam, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, gam2, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   gam2 = gam*gam
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      gtmp = +2*gam*exp(-r2*gam2)/(sqrtpi*r2) - erf(r1*gam)/(r2*r1)
      atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2) + erf(r1*alp)/(r2*r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

subroutine get_damat_rec_3d(rij, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, etmp, dtmp, alp2
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4*pi/vol
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      etmp = fac * exp(-0.25_wp*g2/alp2)/g2
      dtmp = -sin(gv) * etmp
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + etmp * cos(gv) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_rec_3d

subroutine solve(self, mol, cn, dcndr, dcndL, energy, gradient, sigma, qvec, dqdr, dqdL)
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in), contiguous :: cn(:)
   real(wp), intent(in), contiguous, optional :: dcndr(:, :, :)
   real(wp), intent(in), contiguous, optional :: dcndL(:, :, :)
   real(wp), intent(out), contiguous, optional :: qvec(:)
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)
   real(wp), intent(inout), contiguous, optional :: energy(:)
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   integer :: ic, jc, iat, ndim, info
   logical :: grad, cpq, dcn
   real(wp) :: alpha
   integer, allocatable :: ipiv(:)
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :), ainv(:, :)
   real(wp), allocatable :: dxdcn(:), atrace(:, :), dadr(:, :, :), dadL(:, :, :)
   type(wignerseitz_cell_type) :: wsc

   ndim = mol%nat + 1
   if (any(mol%periodic)) then
      call new_wignerseitz_cell(wsc, mol)
      call get_alpha(mol%lattice, alpha)
   end if

   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn

   allocate(amat(ndim, ndim), xvec(ndim))
   allocate(ipiv(ndim))
   if (grad.or.cpq) then
      allocate(dxdcn(ndim))
   end if

   call get_vrhs(self, mol, cn, xvec, dxdcn)
   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, wsc, alpha, amat)
   else
      call get_amat_0d(self, mol, amat)
   end if

   vrhs = xvec
   ainv = amat

   call sytrf(ainv, ipiv, info=info, uplo='l')

   if (info == 0) then
      if (cpq) then
         call sytri(ainv, ipiv, info=info, uplo='l')
         if (info == 0) then
            call symv(ainv, xvec, vrhs, uplo='l')
            do ic = 1, ndim
               do jc = ic+1, ndim
                  ainv(ic, jc) = ainv(jc, ic)
               end do
            end do
         end if
      else
         call sytrs(ainv, vrhs, ipiv, info=info, uplo='l')
      end if
   end if

   if (present(qvec)) then
      qvec(:) = vrhs(:mol%nat)
   end if

   if (present(energy)) then
      call symv(amat, vrhs, xvec, alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
   end if

   if (grad.or.cpq) then
      allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), atrace(3, mol%nat))
      if (any(mol%periodic)) then
         call get_damat_3d(self, mol, wsc, alpha, vrhs, dadr, dadL, atrace)
      else
         call get_damat_0d(self, mol, vrhs, dadr, dadL, atrace)
      end if
      xvec(:) = -dxdcn * vrhs
   end if

   if (grad) then
      call gemv(dadr, vrhs, gradient, beta=1.0_wp)
      call gemv(dcndr, xvec(:mol%nat), gradient, beta=1.0_wp)
      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dcndL, xvec(:mol%nat), sigma, beta=1.0_wp)
   end if

   if (cpq) then
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
         dadr(:, :, iat) = -dcndr(:, :, iat) * dxdcn(iat) + dadr(:, :, iat)
         dadL(:, :, iat) = -dcndL(:, :, iat) * dxdcn(iat) + dadL(:, :, iat)
      end do

      call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
      call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
   end if

end subroutine solve


end module multicharge_model

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_output
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa, autokcal, autoev
   use mctc_io_math, only : matinv_3x3
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_model, only : d4_model
   use dftd4_version, only : get_dftd4_version
   implicit none
   private

   public :: ascii_atomic_radii, ascii_atomic_references, ascii_system_properties
   public :: ascii_results, ascii_damping_param, ascii_pairwise
   public :: turbomole_gradient, turbomole_gradlatt
   public :: json_results, tagged_result


contains


subroutine ascii_atomic_radii(unit, mol, disp)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   integer :: isp

   write(unit, '(a,":")') "Atomic data, radii in Ångström"
   write(unit, '(54("-"))')
   write(unit, '(a4,5x,*(1x,a10))') &
      "Z", "R(cov)", "r4/r2", "hardness", "EN"
   write(unit, '(54("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4, *(1x,f10.4))') &
         & mol%num(isp), mol%sym(isp), &
         & disp%rcov(isp)*autoaa, &
         & disp%r4r2(isp)*autoaa, &
         & disp%eta(isp), &
         & disp%en(isp)
   end do
   write(unit, '(54("-"))')
   write(unit, '(a)')

end subroutine ascii_atomic_radii


subroutine ascii_atomic_references(unit, mol, disp)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   integer :: isp, iref, mref

   mref = maxval(disp%ref)
   write(unit, '(a,":")') "Atomic reference systems (in atomic units)"
   write(unit, '(70("-"))')
   write(unit, '(a4, 5x)', advance='no') "Z"
   do iref = 1, 2
      write(unit, '(a4, 2(1x, a7), 1x, a9)', advance='no') &
         "#", "CN", "q+Z", "C6(AA)"
   end do
   write(unit, '(a)')
   write(unit, '(70("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4)', advance='no') &
         & mol%num(isp), mol%sym(isp)
      do iref = 1, disp%ref(isp)
         write(unit, '(i4, 2(1x, f7.4), 1x, f9.4)', advance='no') &
            iref, disp%cn(iref, isp), disp%q(iref, isp) + disp%zeff(isp), &
            disp%c6(iref, iref, isp, isp)
         if (iref == 2 .and. disp%ref(isp) > 2) then
            write(unit, '(/,9x)', advance='no')
         end if
         if (iref == 4 .and. disp%ref(isp) > 4) then
            write(unit, '(/,9x)', advance='no')
         end if
         if (iref == 6 .and. disp%ref(isp) > 6) then
            write(unit, '(/,9x)', advance='no')
         end if
      end do
      write(unit, '(a)')
   end do
   write(unit, '(70("-"))')
   write(unit, '(a)')

end subroutine ascii_atomic_references


subroutine ascii_system_properties(unit, mol, disp, cn, q, c6)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Coordination numbers
   real(wp), intent(in) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(in) :: q(:)

   !> Atomic dispersion coefficients
   real(wp), intent(in) :: c6(:, :)

   integer :: iat, isp, jat
   real(wp) :: sum_c8

   sum_c8 = 0.0_wp
   write(unit, '(a,":")') "Atomic properties (in atomic units)"
   write(unit, '(61("-"))')
   write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "CN", "q", "C6(AA)", "C8(AA)"
   write(unit, '(61("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6,1x,i4,1x,a4,*(1x,f10.4))') &
         & iat, mol%num(isp), mol%sym(isp), cn(iat), q(iat), c6(iat, iat), &
         & c6(iat, iat)*3*disp%r4r2(isp)**2
      do jat = 1, mol%nat
         sum_c8 = sum_c8 + 3*c6(jat, iat)*disp%r4r2(mol%id(jat))*disp%r4r2(isp)
      end do
   end do
   write(unit, '(61("-"))')
   write(unit, '(a)')

   write(unit, '(a,":")') "Molecular properties (in atomic units)"
   write(unit, '(40("-"))')
   write(unit, '(1x, a, t20, f19.4)') &
      "molecular C6",  sum(c6), &
      "molecular C8",  sum_c8
   write(unit, '(40("-"))')
   write(unit, '(a)')

end subroutine ascii_system_properties


subroutine ascii_results(unit, mol, energy, gradient, sigma)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   real(wp), intent(in) :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)

   integer :: iat, isp
   logical :: grad
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]

   grad = present(gradient) .and. present(sigma)

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Dispersion energy", energy, "Eh"
   write(unit, '(a)')
   if (grad) then
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
         & "Gradient norm", norm2(gradient), "Eh/a0"
      write(unit, '(50("-"))')
      write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "dE/dx", "dE/dy", "dE/dz"
      write(unit, '(50("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         write(unit, '(i6,1x,i4,1x,a4,*(es11.3))') &
            & iat, mol%num(isp), mol%sym(isp), gradient(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')

      write(unit, '(a,":")') &
         & "Virial"
      write(unit, '(50("-"))')
      write(unit, '(a15,1x,*(1x,a10))') "component", "x", "y", "z"
      write(unit, '(50("-"))')
      do iat = 1, 3
         write(unit, '(2x,4x,1x,a4,1x,4x,*(es11.3))') &
            & comp(iat), sigma(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')
   end if

end subroutine ascii_results


subroutine ascii_pairwise(unit, mol, pair_disp2, pair_disp3)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   real(wp), intent(in) :: pair_disp2(:, :)
   real(wp), intent(in) :: pair_disp3(:, :)

   integer :: iat, jat, isp, jsp
   real(wp) :: disp, e2, e3

   e2 = 0.0_wp
   e3 = 0.0_wp

   write(unit, '(a,":")') "Pairwise representation of dispersion (in kcal/mol)"
   write(unit, '(82("-"))')
   write(unit, '(2(a6,1x,a4,5x),*(1x,a10:,1x,a7))') &
      "#", "Z", "#", "Z", "additive", "(rel.)", "non-add.", "(rel.)", "total"
   write(unit, '(82("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         e2 = e2 + pair_disp2(jat, iat)
         e3 = e3 + pair_disp3(jat, iat)
         disp = pair_disp2(jat, iat) + pair_disp3(jat, iat)
         if (abs(disp) < epsilon(disp)) cycle
         write(unit, '(2(i6,1x,i4,1x,a4),*(1x,es10.2:,1x,"(",i4,"%)"))') &
            & iat, mol%num(isp), mol%sym(isp), &
            & jat, mol%num(jsp), mol%sym(jsp), &
            & pair_disp2(jat, iat) * autokcal, nint(pair_disp2(jat, iat)/disp*100), &
            & pair_disp3(jat, iat) * autokcal, nint(pair_disp3(jat, iat)/disp*100), &
            & disp * autokcal
      end do
   end do
   write(unit, '(82("-"))')
   disp = e2 + e3
   write(unit, '(1x, a, t33,*(1x,es10.2:,1x,"(",i4,"%)"))') &
      & "total dispersion energy", &
      & e2 * autokcal, nint(e2/disp*100), &
      & e3 * autokcal, nint(e3/disp*100), &
      & disp * autokcal
   write(unit, '(82("-"))')
   write(unit, '(a)')

end subroutine ascii_pairwise


subroutine ascii_damping_param(unit, param, method)

   !> Unit for output
   integer, intent(in) :: unit

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Method name
   character(len=*), intent(in), optional :: method

   select type(param)
   type is (rational_damping_param)
      write(unit, '(a, ":", 1x)', advance="no") "Rational (Becke-Johnson) damping"
      if (present(method)) then
         write(unit, '(a, "-")', advance="no") method
      end if
      write(unit, '(a)') &
         & trim(merge("D4-ATM", "D4    ", abs(param%s9) > 0))
      write(unit, '(21("-"))')
      write(unit, '(a4, t10, f10.4)') &
         & "s6", param%s6, &
         & "s8", param%s8, &
         & "s9", param%s9, &
         & "a1", param%a1, &
         & "a2", param%a2, &
         & "alp", param%alp
      write(unit, '(20("-"))')
      write(unit, '(a)')
   end select

end subroutine ascii_damping_param


subroutine turbomole_gradlatt(mol, fname, energy, sigma, stat)
   type(structure_type),intent(in) :: mol
   character(len=*),intent(in) :: fname
   real(wp),intent(in) :: energy
   real(wp),intent(in) :: sigma(3,3)
   integer, intent(out) :: stat
   character(len=:),allocatable :: line
   integer  :: i,j,icycle,line_number
   integer  :: err
   integer  :: igrad ! file handle
   logical  :: exist
   real(wp) :: escf
   real(wp) :: glat(3,3), inv_lat(3,3), gradlatt(3, 3)
   real(wp) :: dlat(3,3)
   stat = 0

   inv_lat = matinv_3x3(mol%lattice)

   do i = 1, 3
      do j = 1, 3
         gradlatt(i,j) = sigma(i,1)*inv_lat(j,1) &
            & + sigma(i,2)*inv_lat(j,2) &
            & + sigma(i,3)*inv_lat(j,3)
      enddo
   enddo

   icycle = 1
   i = 0
   escf = 0.0_wp

   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=igrad,file=fname)
      read_file: do
         call getline(igrad,line,iostat=err)
         if (err.ne.0) exit read_file
         i=i+1
         if (index(line,'cycle') > 0) line_number = i
      enddo read_file
      if (line_number < 2) then
         stat = 1
         return
      endif

      rewind(igrad)
      skip_lines: do i = 1, line_number-1
         read(igrad,'(a)')
      enddo skip_lines
      call getline(igrad,line)
      read(line(10:17),*,iostat=err) icycle
      read(line(33:51),*,iostat=err) escf

      do i = 1, 3
         call getline(igrad,line)
         read(line,*,iostat=err) dlat(1,i),dlat(2,i),dlat(3,i)
      enddo
      if (any(abs(dlat-mol%lattice) > 1.0e-8_wp)) then
         stat = 1
         return
      endif
      do i = 1, 3
         call getline(igrad,line)
         read(line,*,iostat=err) glat(1,i),glat(2,i),glat(3,i)
      enddo
      do i = 1, 3
         backspace(igrad)
         backspace(igrad)
      enddo
      backspace(igrad)
   else
      open(newunit=igrad,file=fname)
      write(igrad,'("$gradlatt")')
   endif

   write(igrad,'(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
                   '"|dE/dlatt| =",f10.6)') &
      icycle, energy+escf, norm2(gradlatt+glat)
   do i = 1, 3
      write(igrad,'(3(F20.14,2x))') mol%lattice(1,i),mol%lattice(2,i),mol%lattice(3,i)
   enddo
   do i = 1, 3
      write(igrad,'(3D22.13)') gradlatt(1,i)+glat(1,i),gradlatt(2,i)+glat(2,i),gradlatt(3,i)+glat(3,i)
   enddo
   write(igrad,'("$end")')
   close(igrad)

end subroutine turbomole_gradlatt


subroutine turbomole_gradient(mol, fname, energy, gradient, stat)
   type(structure_type),intent(in) :: mol
   character(len=*),intent(in) :: fname
   real(wp),intent(in) :: energy
   real(wp),intent(in) :: gradient(:, :)
   integer, intent(out) :: stat
   character(len=:),allocatable :: line
   integer  :: i,icycle,line_number
   integer  :: err
   integer  :: igrad ! file handle
   logical  :: exist
   real(wp) :: escf
   real(wp),allocatable :: gscf(:,:)
   real(wp),allocatable :: xyz (:,:)
   allocate( gscf(3,mol%nat), source = 0.0_wp )
   stat = 0
   icycle = 1
   i = 0
   escf = 0.0_wp

   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=igrad,file=fname)
      read_file: do
         call getline(igrad,line,iostat=err)
         if (err.ne.0) exit read_file
         i=i+1
         if (index(line,'cycle') > 0) line_number = i
      enddo read_file
      if (line_number < 2) then
         stat = 1
         return
      endif

      rewind(igrad)
      skip_lines: do i = 1, line_number-1
         read(igrad,'(a)')
      enddo skip_lines
      call getline(igrad,line)
      read(line(10:17),*,iostat=err) icycle
      read(line(33:51),*,iostat=err) escf

      allocate(xyz(3,mol%nat))
      do i = 1, mol%nat
         call getline(igrad,line)
         read(line,*,iostat=err) xyz(1,i),xyz(2,i),xyz(3,i)
      enddo
      if (any(abs(xyz-mol%xyz) > 1.0e-8_wp)) then
         stat = 1
         return
      endif
      do i = 1, mol%nat
         call getline(igrad,line)
         read(line,*,iostat=err) gscf(1,i),gscf(2,i),gscf(3,i)
      enddo
      do i = 1, mol%nat
         backspace(igrad)
         backspace(igrad)
      enddo
      backspace(igrad)
   else
      open(newunit=igrad,file=fname)
      write(igrad,'("$grad")')
   endif

   write(igrad,'(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
                   '"|dE/dxyz| =",f10.6)') &
      icycle, energy+escf, norm2(gradient+gscf)
   do i = 1, mol%nat
      write(igrad,'(3(F20.14,2x),4x,a2)') mol%xyz(1,i),mol%xyz(2,i),mol%xyz(3,i),mol%sym(i)
   enddo
   do i = 1, mol%nat
      write(igrad,'(3D22.13)') gradient(1,i)+gscf(1,i),gradient(2,i)+gscf(2,i),gradient(3,i)+gscf(3,i)
   enddo
   write(igrad,'("$end")')
   close(igrad)

end subroutine turbomole_gradient


!> reads a line from unit into an allocatable character
subroutine getline(unit,line,iostat)
   integer,intent(in) :: unit
   character(len=:),allocatable,intent(out) :: line
   integer,intent(out),optional :: iostat

   integer,parameter  :: buffersize=256
   character(len=buffersize) :: buffer
   integer :: size
   integer :: stat

   line = ''
   do
      read(unit,'(a)',advance='no',iostat=stat,size=size)  &
      &    buffer
      if (stat.gt.0) then
         if (present(iostat)) iostat=stat
         return ! an error occurred
      endif
      line = line // buffer(:size)
      if (stat.lt.0) then
         if (is_iostat_eor(stat)) stat = 0
         if (present(iostat)) iostat=stat
         return
      endif
   enddo

end subroutine getline


subroutine json_results(unit, indentation, energy, gradient, sigma, cn, q, c6, alpha, &
      & pairwise_energy2, pairwise_energy3)
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: q(:)
   real(wp), intent(in), optional :: c6(:, :)
   real(wp), intent(in), optional :: alpha(:)
   real(wp), intent(in), optional :: pairwise_energy2(:, :)
   real(wp), intent(in), optional :: pairwise_energy3(:, :)
   character(len=:), allocatable :: indent, version_string
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)

   call get_dftd4_version(string=version_string)

   if (present(indentation)) then
      indent = indentation
   end if

   write(unit, '("{")', advance='no')
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, jsonkey, advance='no') 'version'
   write(unit, '(1x,a)', advance='no') '"'//version_string//'"'
   if (present(energy)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'energy'
      write(unit, '(1x,es25.16)', advance='no') energy
   end if
   if (present(sigma)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'virial'
      array = reshape(sigma, [size(sigma)])
      call write_json_array(unit, array, indent)
   end if
   if (present(gradient)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'gradient'
      array = reshape(gradient, [size(gradient)])
      call write_json_array(unit, array, indent)
   end if
   if (present(cn)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'coordination numbers'
      call write_json_array(unit, cn, indent)
   end if
   if (present(q)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'partial charges'
      call write_json_array(unit, q, indent)
   end if
   if (present(c6)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'c6 coefficients'
      array = reshape(c6, [size(c6)])
      call write_json_array(unit, array, indent)
   end if
   if (present(alpha)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'polarizibilities'
      call write_json_array(unit, alpha, indent)
   end if
   if (present(pairwise_energy2)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'additive pairwise energy'
      array = reshape(pairwise_energy2, [size(pairwise_energy2)])
      call write_json_array(unit, array, indent)
   end if
   if (present(pairwise_energy3)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'non-additive pairwise energy'
      array = reshape(pairwise_energy3, [size(pairwise_energy3)])
      call write_json_array(unit, array, indent)
   end if
   if (allocated(indent)) write(unit, '(/)', advance='no')
   write(unit, '("}")')

end subroutine json_results


subroutine write_json_array(unit, array, indent)
   integer, intent(in) :: unit
   real(wp), intent(in) :: array(:)
   character(len=:), allocatable, intent(in) :: indent
   integer :: i
   write(unit, '("[")', advance='no')
   do i = 1, size(array)
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
      write(unit, '(es23.16)', advance='no') array(i)
      if (i /= size(array)) write(unit, '(",")', advance='no')
   end do
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("]")', advance='no')
end subroutine write_json_array


subroutine tagged_result(unit, energy, gradient, sigma)
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   character(len=*), parameter :: tag_header = &
      & '(a,t20,":",a,":",i0,":",*(i0:,","))'

   if (present(energy)) then
      write(unit, tag_header) "energy", "real", 0
      write(unit, '(3es24.16)') energy
   end if
   if (present(gradient)) then
      write(unit, tag_header) "gradient", "real", 2, 3, size(gradient, 2)
      write(unit, '(3es24.16)') gradient
   end if
   if (present(sigma)) then
      write(unit, tag_header) "virial", "real", 2, 3, 3
      write(unit, '(3es24.16)') sigma
   end if

end subroutine tagged_result

end module dftd4_output

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_param
   use mctc_env, only : wp
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   implicit none
   private

   public :: get_rational_damping


   enum, bind(C)
      enumerator :: p_invalid, &
         & p_hf, p_blyp, p_bpbe, p_bp, p_bpw, p_lb94, p_mpwlyp, p_mpwpw, &
         & p_olyp, p_opbe, p_pbe, p_rpbe, p_revpbe, p_pw86pbe, &
         & p_rpw86pbe, p_pw91, p_pwp, p_xlyp, p_b97, p_tpss, p_revtpss, &
         & p_scan, p_rscan, p_r2scan, p_b1lyp, p_b3lyp, p_bhlyp, p_b1p, &
         & p_b3p, p_b1pw, p_b3pw, p_o3lyp, p_revpbe0, p_revpbe38, &
         & p_pbe0, p_pwp1, p_pw1pw, p_mpw1pw, p_mpw1lyp, p_pw6b95, &
         & p_tpssh, p_tpss0, p_x3lyp, p_m06l, p_m06, p_m062x, p_b97d, &
         & p_wb97, p_wb97x, p_b97m, p_wb97m, p_camb3lyp, p_lcblyp, &
         & p_lh07tsvwn, p_lh07ssvwn, p_lh12ctssirpw92, p_lh12ctssifpw92, &
         & p_lh14tcalpbe, p_lh20t, &
         & p_b2plyp, p_b2gpplyp, p_mpw2plyp, p_pwpb95, &
         & p_dsdblyp, p_dsdpbe, p_dsdpbeb95, p_dsdpbep86, p_dsdsvwn, &
         & p_dodblyp, p_dodpbe, p_dodpbeb95, p_dodpbep86, p_dodsvwn, &
         & p_pbe0_2, p_pbe0_dh, p_hf3c, p_hf3cv, p_pbeh3c, p_b973c, &
         & p_hsesol, p_pwgga, p_dftb_3ob, p_dftb_mio, p_dftb_ob2, &
         & p_dftb_matsci, p_dftb_pbc, p_hcth120, p_ptpss, p_lcwpbe, &
         & p_bmk, p_b1b95, p_pwb6k, p_otpss, p_ssb, p_revssb, &
         & p_pbesol, p_hse06, p_pbexalpha, p_pbehpbe, p_hcth407, &
         & p_n12, p_pkzb, p_thcth, p_m11l, p_mn15l, p_mpwb1k, &
         & p_mpw1kcis, p_mpwkcis1k, p_pbeh1pbe, p_pbe1kcis, p_b97_1, &
         & p_b97_2, p_b98, p_hiss, p_hse03, p_revtpssh, p_tpss1kcis, &
         & p_m05, p_m052x, p_m08hx, p_lcwhpbe, p_mn12l, p_tauhcthhyb, &
         & p_sogga11x, p_n12sx, p_mn12sx, p_mn15, p_glyp, p_bop, &
         & p_mpw1b95, p_revpbe0dh, p_revtpss0, p_revdsdpbep86, p_revdsdpbe, &
         & p_revdsdblyp, p_revdodpbep86
   end enum
   integer, parameter :: df_enum = kind(p_invalid)

contains

subroutine get_rational_damping(functional, param, s9)
   character(len=*), intent(in) :: functional
   class(damping_param), allocatable, intent(out) :: param
   real(wp), intent(in), optional :: s9

   character(len=:), allocatable :: fname
   integer :: is, id
   logical :: mbd

   mbd = merge(s9 /= 0.0_wp, .true., present(s9))

   is = index(functional, '/')
   if (is == 0) is = len_trim(functional) + 1
   fname = lowercase(functional(:is-1))

   id = get_functional_id(fname)

   if (mbd) then
      call get_d4eeq_bjatm_parameter(id, param, s9)
      if (.not.allocated(param)) then
         call get_d4eeq_bj_parameter(id, param, s9)
      end if
   else
      call get_d4eeq_bj_parameter(id, param, s9)
      if (.not.allocated(param)) then
         call get_d4eeq_bjatm_parameter(id, param, s9)
      end if
   end if

end subroutine get_rational_damping

subroutine get_d4eeq_bj_parameter(dfnum, param, s9)
   integer(df_enum), intent(in) :: dfnum
   class(damping_param), allocatable, intent(out) :: param
   real(wp), intent(in), optional :: s9
   select case(dfnum)
   case(p_dftb_3ob)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=0.4727337_wp, a1=0.5467502_wp, a2=4.4955068_wp)
   case(p_dftb_matsci)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.7711819_wp, a1=0.4681712_wp, a2=5.2918629_wp)
   case(p_dftb_mio)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=1.1948145_wp, a1=0.6074567_wp, a2=4.9336133_wp)
   case(p_dftb_ob2)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.7611320_wp, a1=0.6037249_wp, a2=5.3900004_wp)
   case(p_dftb_pbc)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=1.7303734_wp, a1=0.5546548_wp, a2=4.7973454_wp)
   end select

contains

   pure function dftd_param(s6, s8, a1, a2, alp) result(param)
      real(wp), intent(in) :: s8, a1, a2
      real(wp), intent(in), optional :: s6, alp
      type(rational_damping_param) :: param
      param = rational_damping_param(&
         & s6=merge(s6, 1.0_wp, present(s6)), &
         & s8=s8, a1=a1, a2=a2, &
         & s9=merge(s9, 0.0_wp, present(s9)), &
         & alp=merge(alp, 16.0_wp, present(alp)))
   end function dftd_param

end subroutine get_d4eeq_bj_parameter

subroutine get_d4eeq_bjatm_parameter(dfnum, param, s9)
   integer(df_enum), intent(in) :: dfnum
   class(damping_param), allocatable, intent(out) :: param
   real(wp), intent(in), optional :: s9
   select case(dfnum)
   case(p_b1b95)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.27701162_wp, a1=0.40554715_wp, a2=4.63323074_wp )
      !  Fitset: MD= 0.22852 MAD= 0.35189 RMSD= 0.46982
   case(p_b1lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.98553711_wp, a1=0.39309040_wp, a2=4.55465145_wp )
      !  Fitset: MD= -0.04797 MAD= 0.25597 RMSD= 0.38778
   case(p_b1p)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.36115015_wp, a1=0.48665293_wp, a2=5.05219572_wp )
      !  Fitset: MD= -0.01406 MAD= 0.27441 RMSD= 0.47328
   case(p_b1pw)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=3.02227550_wp, a1=0.47396846_wp, a2=4.49845309_wp )
      !  Fitset: MD= 0.10485 MAD= 0.32175 RMSD= 0.48508
   case(p_b2gpplyp)
      param = dftd_param ( & ! (SAW190107)
         &  s6=0.5600_wp, s8=0.94633372_wp, a1=0.42907301_wp, a2=5.18802602_wp )
      !  Fitset: MD= -0.05248 MAD= 0.18110 RMSD= 0.27365
   case(p_b2plyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.6400_wp, s8=1.16888646_wp, a1=0.44154604_wp, a2=4.73114642_wp )
      !  Fitset: MD= -0.03761 MAD= 0.18247 RMSD= 0.27109
   case(p_b3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.02929367_wp, a1=0.40868035_wp, a2=4.53807137_wp )
      !  Fitset: MD= -0.05892 MAD= 0.26117 RMSD= 0.40531
   case(p_b3p)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.08822155_wp, a1=0.47324238_wp, a2=4.98682134_wp )
      !  Fitset: MD= -0.02970 MAD= 0.26962 RMSD= 0.46761
   case(p_b3pw)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=2.88364295_wp, a1=0.46990860_wp, a2=4.51641422_wp )
      !  Fitset: MD= 0.06643 MAD= 0.29151 RMSD= 0.45541
   case(p_b97)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.87854260_wp, a1=0.29319126_wp, a2=4.51647719_wp )
      !  Fitset: MD= -0.13017 MAD= 0.24778 RMSD= 0.36116
   case(p_bhlyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.65281646_wp, a1=0.27263660_wp, a2=5.48634586_wp )
      !  Fitset: MD= -0.15832 MAD= 0.34132 RMSD= 0.57342
   case(p_blyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.34076671_wp, a1=0.44488865_wp, a2=4.09330090_wp )
      !  Fitset: MD= 0.04801 MAD= 0.28161 RMSD= 0.38321
   case(p_bpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.64405246_wp, a1=0.52905620_wp, a2=4.11311891_wp )
      !  Fitset: MD= 0.19316 MAD= 0.41912 RMSD= 0.60452
   case(p_bp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.35497927_wp, a1=0.43645861_wp, a2=4.92406854_wp )
      !  Fitset: MD= 0.08252 MAD= 0.32681 RMSD= 0.47063
   case(p_bpw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.24571506_wp, a1=0.50050454_wp, a2=4.12346483_wp )
      !  Fitset: MD= 0.20607 MAD= 0.41941 RMSD= 0.59589
   case(p_camb3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.66041301_wp, a1=0.40267156_wp, a2=5.17432195_wp )
      !  Fitset: MD= -0.19675 MAD= 0.34901 RMSD= 0.59087
   case(p_dodblyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4700_wp, s8=1.31146043_wp, a1=0.43407294_wp, a2=4.27914360_wp )
      !  Fitset: MD= 0.03323 MAD= 0.13858 RMSD= 0.20861
   case(p_dodpbeb95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5600_wp, s8=0.01574635_wp, a1=0.43745720_wp, a2=3.69180763_wp )
      !  Fitset: MD= 0.03704 MAD= 0.13343 RMSD= 0.18278
   case(p_dodpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4800_wp, s8=0.92051454_wp, a1=0.43037052_wp, a2=4.38067238_wp )
      !  Fitset: MD= 0.01065 MAD= 0.13414 RMSD= 0.21424
   case(p_dodpbep86)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4600_wp, s8=0.71405681_wp, a1=0.42408665_wp, a2=4.52884439_wp )
      !  Fitset: MD= -0.03740 MAD= 0.12467 RMSD= 0.18127
   case(p_dodsvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4200_wp, s8=0.94500207_wp, a1=0.47449026_wp, a2=5.05316093_wp )
      !  Fitset: MD= -0.07427 MAD= 0.16970 RMSD= 0.25286
   case(p_dsdblyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5400_wp, s8=0.63018237_wp, a1=0.47591835_wp, a2=4.73713781_wp )
      !  Fitset: MD= -0.01981 MAD= 0.14823 RMSD= 0.21530
   case(p_dsdpbeb95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5400_wp, s8=-0.14668670_wp, a1=0.46394587_wp, a2=3.64913860_wp )
      !  Fitset: MD= 0.02996 MAD= 0.12414 RMSD= 0.16860
   case(p_dsdpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4500_wp, s8=0.70584116_wp, a1=0.45787085_wp, a2=4.44566742_wp )
      !  Fitset: MD= 0.00866 MAD= 0.13406 RMSD= 0.21380
   case(p_dsdpbep86)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4700_wp, s8=0.37586675_wp, a1=0.53698768_wp, a2=5.13022435_wp )
      !  Fitset: MD= -0.05273 MAD= 0.14259 RMSD= 0.21271
   case(p_dsdsvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4100_wp, s8=0.72914436_wp, a1=0.51347412_wp, a2=5.11858541_wp )
      !  Fitset: MD= -0.08974 MAD= 0.32285 RMSD= 0.43146
   case(p_glyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=4.23798924_wp, a1=0.38426465_wp, a2=4.38412863_wp )
      !  Fitset: MD= 0.63466 MAD= 0.89568 RMSD= 1.11309
   case(p_hf)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.61679827_wp, a1=0.44959224_wp, a2=3.35743605_wp )
      !  Fitset: MD= -0.02597 MAD= 0.34732 RMSD= 0.49719
   case(p_lb94)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.59538499_wp, a1=0.42088944_wp, a2=3.28193223_wp )
      !  Fitset: MD= 0.31701 MAD= 0.53196 RMSD= 0.74553
   case(p_lcblyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.60344180_wp, a1=0.45769839_wp, a2=7.86924893_wp )
      !  Fitset: MD= -0.39724 MAD= 0.72327 RMSD= 1.18218
   case(p_lh07ssvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.16675531_wp, a1=0.35965552_wp, a2=4.31947614_wp )
      !  Fitset: MD= 0.32224 MAD= 0.59006 RMSD= 0.86272
   case(p_lh07tsvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.09333001_wp, a1=0.35025189_wp, a2=4.34166515_wp )
      !  Fitset: MD= 0.24243 MAD= 0.43497 RMSD= 0.61671
   case(p_lh12ctssifpw92)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.68467610_wp, a1=0.34190416_wp, a2=3.91039666_wp )
      !  Fitset: MD= 0.55106 MAD= 0.80783 RMSD= 1.11048
   case(p_lh12ctssirpw92)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.48973402_wp, a1=0.34026075_wp, a2=3.96948081_wp )
      !  Fitset: MD= 0.47785 MAD= 0.71188 RMSD= 0.98422
   case(p_lh14tcalpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.28130770_wp, a1=0.38822021_wp, a2=4.92501211_wp )
      !  Fitset: MD= -0.02105 MAD= 0.22968 RMSD= 0.36045
   case(p_lh20t)
      param = dftd_param ( & ! (10.1021/acs.jctc.0c00498)
         & s6=1.000_wp, s8=0.113_wp, a1=0.479_wp, a2=4.635_wp )
   case(p_m06)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.16366729_wp, a1=0.53456413_wp, a2=6.06192174_wp )
      !  Fitset: MD= 0.01788 MAD= 0.24914 RMSD= 0.38604
   case(p_m06l)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.59493760_wp, a1=0.71422359_wp, a2=6.35314182_wp )
      !  Fitset: MD= 0.08395 MAD= 0.24888 RMSD= 0.34879
   case(p_mpw1b95)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=0.50093024_wp, a1=0.41585097_wp, a2=4.99154869_wp )
      !  Fitset: MD= 0.00585 MAD= 0.15695 RMSD= 0.21297
   case(p_mpw1lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.15591153_wp, a1=0.25603493_wp, a2=5.32083895_wp )
      !  Fitset: MD= -0.26979 MAD= 0.41542 RMSD= 0.60678
   case(p_mpw1pw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.80841716_wp, a1=0.42961819_wp, a2=4.68892341_wp )
      !  Fitset: MD= -0.08840 MAD= 0.26815 RMSD= 0.45231
   case(p_mpw2plyp)
      param = dftd_param ( & ! (SAW190107)
         &  s6=0.7500_wp, s8=0.45788846_wp, a1=0.42997704_wp, a2=5.07650682_wp )
      !  Fitset: MD= -0.18921 MAD= 0.30115 RMSD= 0.44049
   case(p_mpwb1k)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=0.57338313_wp, a1=0.44687975_wp, a2=5.21266777_wp )
      !  Fitset: MD= -0.00870 MAD= 0.17226 RMSD= 0.23614
   case(p_mpwlyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.25842942_wp, a1=0.25773894_wp, a2=5.02319542_wp )
      !  Fitset: MD= -0.24426 MAD= 0.39145 RMSD= 0.54503
   case(p_mpwpw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.82596836_wp, a1=0.34526745_wp, a2=4.84620734_wp )
      !  Fitset: MD= -0.06278 MAD= 0.27913 RMSD= 0.43988
   case(p_o3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.75762508_wp, a1=0.10348980_wp, a2=6.16233282_wp )
      !  Fitset: MD= -0.19268 MAD= 0.38577 RMSD= 0.62168
   case(p_olyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.74836820_wp, a1=0.60184498_wp, a2=2.53292167_wp )
      !  Fitset: MD= 0.12352 MAD= 0.37113 RMSD= 0.58291
   case(p_opbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.06917417_wp, a1=0.68267534_wp, a2=2.22849018_wp )
      !  Fitset: MD= 0.26699 MAD= 0.55308 RMSD= 0.85023
   case(p_pbe0_2)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5000_wp, s8=0.64299082_wp, a1=0.76542115_wp, a2=5.78578675_wp )
      !  Fitset: MD= -0.04260 MAD= 0.21186 RMSD= 0.34045
   case(p_pbe0)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.20065498_wp, a1=0.40085597_wp, a2=5.02928789_wp )
      !  Fitset: MD= -0.17892 MAD= 0.30557 RMSD= 0.51050
   case(p_pbe0_dh)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.8750_wp, s8=0.96811578_wp, a1=0.47592488_wp, a2=5.08622873_wp )
      !  Fitset: MD= -0.13857 MAD= 0.27919 RMSD= 0.47256
   case(p_pbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
      !  Fitset: MD= -0.20544 MAD= 0.33635 RMSD= 0.51168
   case(p_pw1pw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.96850170_wp, a1=0.42427511_wp, a2=5.02060636_wp )
      !  Fitset: MD= -0.27325 MAD= 0.42206 RMSD= 0.64119
   case(p_pw6b95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.31926054_wp, a1=0.04142919_wp, a2=5.84655608_wp )
      !  Fitset: MD= -0.04767 MAD= 0.14330 RMSD= 0.18958
   case(p_pw86pbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.21362856_wp, a1=0.40510366_wp, a2=4.66737724_wp )
      !  Fitset: MD= -0.11505 MAD= 0.24691 RMSD= 0.38101
   case(p_pw91)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.77283111_wp, a1=0.39581542_wp, a2=4.93405761_wp )
      !  Fitset: MD= -0.33019 MAD= 0.48611 RMSD= 0.68110
   case(p_pwp1)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.60492565_wp, a1=0.46855837_wp, a2=5.76921413_wp )
      !  Fitset: MD= -0.35321 MAD= 0.54026 RMSD= 0.86629
   case(p_pwpb95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.8200_wp, s8=-0.34639127_wp, a1=0.41080636_wp, a2=3.83878274_wp )
      !  Fitset: MD= 0.02143 MAD= 0.13040 RMSD= 0.17599
   case(p_pwp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.32801227_wp, a1=0.35874687_wp, a2=6.05861168_wp )
      !  Fitset: MD= -0.42482 MAD= 0.62607 RMSD= 0.91840
   case(p_revpbe0)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.57185414_wp, a1=0.38705966_wp, a2=4.11028876_wp )
      !  Fitset: MD= 0.02724 MAD= 0.21587 RMSD= 0.36040
   case(p_revpbe0dh)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.8750_wp, s8=1.24456037_wp, a1=0.36730560_wp, a2=4.71126482_wp )
      !  Fitset: MD= -0.01089 MAD= 0.20910 RMSD= 0.33564
   case(p_revpbe38)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.66597472_wp, a1=0.39476833_wp, a2=4.39026628_wp )
      !  Fitset: MD= -0.01326 MAD= 0.22598 RMSD= 0.36210
   case(p_revpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.74676530_wp, a1=0.53634900_wp, a2=3.07261485_wp )
      !  Fitset: MD= 0.05649 MAD= 0.25212 RMSD= 0.40863
   case(p_revtpss0)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.54664499_wp, a1=0.45890964_wp, a2=4.78426405_wp )
      !  Fitset: MD= -0.05298 MAD= 0.19965 RMSD= 0.32081
   case(p_revtpss)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.53089454_wp, a1=0.44880597_wp, a2=4.64042317_wp )
      !  Fitset: MD= -0.01904 MAD= 0.19568 RMSD= 0.29618
   case(p_revtpssh)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.52740307_wp, a1=0.45161957_wp, a2=4.70779483_wp )
      !  Fitset: MD= -0.03731 MAD= 0.19133 RMSD= 0.29091
   case(p_rpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.31183787_wp, a1=0.46169493_wp, a2=3.15711757_wp )
      !  Fitset: MD= -0.07156 MAD= 0.26348 RMSD= 0.38671
   case(p_rpw86pbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.12624034_wp, a1=0.38151218_wp, a2=4.75480472_wp )
      !  Fitset: MD= -0.12740 MAD= 0.26294 RMSD= 0.40614
   case(p_scan)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.46126056_wp, a1=0.62930855_wp, a2=6.31284039_wp )
      !  Fitset: MD= -0.13170 MAD= 0.28640 RMSD= 0.51183
   case(p_rscan)
      param = dftd_param ( & ! (10.1063/5.0041008)
         &  s6=1.0000_wp, s8=0.87728975_wp, a1=0.49116966_wp, a2=5.75859346_wp )
   case(p_r2scan)
      param = dftd_param ( & ! (10.1063/5.0041008)
         &  s6=1.0000_wp, s8=0.60187490_wp, a1=0.51559235_wp, a2=5.77342911_wp )
   case(p_tpss0)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.62438102_wp, a1=0.40329022_wp, a2=4.80537871_wp )
      !  Fitset: MD= -0.09569 MAD= 0.26733 RMSD= 0.44767
   case(p_tpss)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
      !  Fitset: MD= -0.09296 MAD= 0.27505 RMSD= 0.42537
   case(p_tpssh)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.85897750_wp, a1=0.44286966_wp, a2=4.60230534_wp )
      !  Fitset: MD=  0.02238 MAD= 0.16042 RMSD= 0.33519
   case(p_b97d)
      param = dftd_param ( & ! (SAW201029)
         &  s6=1.0000_wp, s8=1.69460052_wp, a1=0.28904684_wp, a2=4.13407323_wp )
      !  Fitset: MD= -0.09858 MAD= 0.26757 RMSD= 0.42380
   case(p_wb97)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=6.55792598_wp, a1=0.76666802_wp, a2=8.36027334_wp )
      !  Fitset: MD= -0.12779 MAD= 0.36152 RMSD= 0.49991
   case(p_wb97x)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.07519516_wp, a1=0.45094893_wp, a2=6.78425255_wp )
      !  S22x5: MD= 0.05 MAD= 0.16 RMSD= 0.22
      !  S66x8: MD= 0.06 MAD= 0.16 RMSD= 0.21
      !  NCI10: MD= 0.08 MAD= 0.15 RMSD= 0.25
   case(p_b97m)
      param = dftd_param ( & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.6633_wp, a1=0.4288_wp, a2=3.9935_wp )
      !  S22x5: MD= 0.03 MAD= 0.12 RMSD= 0.18
      !  S66x8: MD= 0.09 MAD= 0.17 RMSD= 0.22
      !  NCI10: MD= 0.09 MAD= 0.15 RMSD= 0.32
   case(p_wb97m)
      param = dftd_param ( & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.7761_wp, a1=0.7514_wp, a2=2.7099_wp )
      !  Fitset: MD= -0.20216 MAD= 0.34696 RMSD= 0.53641
   case(p_x3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.54701429_wp, a1=0.20318443_wp, a2=5.61852648_wp )
      !  Fitset: MD= -0.15607 MAD= 0.31342 RMSD= 0.49546
   case(p_xlyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.62972054_wp, a1=0.11268673_wp, a2=5.40786417_wp )
      !  Fitset: MD= -0.03900 MAD= 0.27562 RMSD= 0.38491
   case(p_revdsdpbep86)
      param = dftd_param ( & ! (WTMAD2)
         &  s6=0.5132_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_revdsdpbe)
      param = dftd_param ( & ! (WTMAD2)
         &  s6=0.6706_wp, s8=0.00000000_wp, a1=0.40000000_wp, a2=3.60000000_wp )
   case(p_revdsdblyp)
      param = dftd_param ( & !(WTMAD2)
         &  s6=0.6141_wp, s8=0.00000000_wp, a1=0.38000000_wp, a2=3.52000000_wp )
   case(p_revdodpbep86)
      param = dftd_param ( & !(WTMAD2)
         &  s6=0.5552_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_dftb_3ob)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=0.6635015_wp, a1=0.5523240_wp, a2=4.3537076_wp)
   case(p_dftb_matsci)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=3.3157614_wp, a1=0.4826330_wp, a2=5.3811976_wp)
   case(p_dftb_mio)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=1.2916225_wp, a1=0.5965326_wp, a2=4.8778602_wp)
   case(p_dftb_ob2)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.9692689_wp, a1=0.6068916_wp, a2=5.4476789_wp)
   case(p_dftb_pbc)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.1667394_wp, a1=0.5646391_wp, a2=4.9576353_wp)
   end select

contains

   pure function dftd_param(s6, s8, a1, a2, alp) result(param)
      real(wp), intent(in) :: s8, a1, a2
      real(wp), intent(in), optional :: s6, alp
      type(rational_damping_param) :: param
      param = rational_damping_param(&
         & s6=merge(s6, 1.0_wp, present(s6)), &
         & s8=s8, a1=a1, a2=a2, &
         & s9=merge(s9, 1.0_wp, present(s9)), &
         & alp=merge(alp, 16.0_wp, present(alp)))
   end function dftd_param

end subroutine get_d4eeq_bjatm_parameter

!> Get the unique identifier for most functionals, returns none if
!> the functional was not known at the time I implemented this mapping
pure function get_functional_id(df) result(num)
   integer(df_enum) :: num
   character(len=*), intent(in) :: df
   select case(df)
   case default
      num = p_invalid
   case('hf')
      num = p_hf
   case('b-lyp', 'blyp')
      num = p_blyp
   case('bpbe')
      num = p_bpbe
   case('b-p', 'bp86', 'bp', 'b-p86')
      num = p_bp
   case('bpw', 'b-pw')
      num = p_bpw
   case('lb94')
      num = p_lb94
   case('mpwlyp', 'mpw-lyp')
      num = p_mpwlyp
   case('mpwpw', 'mpw-pw', 'mpwpw91')
      num = p_mpwpw
   case('o-lyp', 'olyp')
      num = p_olyp
   case('opbe')
      num = p_opbe
   case('pbe')
      num = p_pbe
   case('rpbe')
      num = p_rpbe
   case('revpbe')
      num = p_revpbe
   case('pw86pbe')
      num = p_pw86pbe
   case('rpw86pbe')
      num = p_rpw86pbe
   case('pw91')
      num = p_pw91
   case('pwp', 'pw-p', 'pw91p86')
      num = p_pwp
   case('x-lyp', 'xlyp')
      num = p_xlyp
   case('b97')
      num = p_b97
   case('tpss')
      num = p_tpss
   case('revtpss')
      num = p_revtpss
   case('scan')
      num = p_scan
   case('rscan')
      num = p_rscan
   case('r2scan', 'r²scan')
      num = p_r2scan
   case('b1lyp', 'b1-lyp')
      num = p_b1lyp
   case('b3-lyp', 'b3lyp')
      num = p_b3lyp
   case('bh-lyp', 'bhlyp')
      num = p_bhlyp
   case('b1p', 'b1-p', 'b1p86')
      num = p_b1p
   case('b3p', 'b3-p', 'b3p86')
      num = p_b3p
   case('b1pw', 'b1-pw', 'b1pw91')
      num = p_b1pw
   case('b3pw', 'b3-pw', 'b3pw91')
      num = p_b3pw
   case('o3-lyp', 'o3lyp')
      num = p_o3lyp
   case('revpbe0')
      num = p_revpbe0
   case('revpbe38')
      num = p_revpbe38
   case('pbe0')
      num = p_pbe0
   case('pwp1')
      num = p_pwp1
   case('pw1pw', 'pw1-pw')
      num = p_pw1pw
   case('mpw1pw', 'mpw1-pw', 'mpw1pw91')
      num = p_mpw1pw
   case('mpw1lyp', 'mpw1-lyp')
      num = p_mpw1lyp
   case('pw6b95')
      num = p_pw6b95
   case('tpssh')
      num = p_tpssh
   case('tpss0')
      num = p_tpss0
   case('x3-lyp', 'x3lyp')
      num = p_x3lyp
   case('m06l')
      num = p_m06l
   case('m06')
      num = p_m06
   case('m06-2x', 'm062x')
      num = p_m062x
   case('wb97', 'ωb97', 'omegab97')
      num = p_wb97
   case('wb97x', 'ωb97x', 'omegab97x')
      num = p_wb97x
   case('cam-b3lyp')
      num = p_camb3lyp
   case('lc-blyp')
      num = p_lcblyp
   case('lh07tsvwn', 'lh07t-svwn')
      num = p_lh07tsvwn
   case('lh07ssvwn', 'lh07s-svwn')
      num = p_lh07ssvwn
   case('lh12ctssirpw92', 'lh12ct-ssirpw92')
      num = p_lh12ctssirpw92
   case('lh12ctssifpw92', 'lh12ct-ssifpw92')
      num = p_lh12ctssifpw92
   case('lh14tcalpbe', 'lh14t-calpbe')
      num = p_lh14tcalpbe
   case('lh20t')
      num = p_lh20t
   case('b2plyp', 'b2-plyp')
      num = p_b2plyp
   case('b2gpplyp', 'b2gp-plyp')
      num = p_b2gpplyp
   case('mpw2plyp')
      num = p_mpw2plyp
   case('pwpb95')
      num = p_pwpb95
   case('dsdblyp', 'dsd-blyp')
      num = p_dsdblyp
   case('dsdpbe', 'dsd-pbe')
      num = p_dsdpbe
   case('dsdpbeb95', 'dsd-pbeb95')
      num = p_dsdpbeb95
   case('dsdpbep86', 'dsd-pbep86')
      num = p_dsdpbep86
   case('dsdsvwn', 'dsd-svwn')
      num = p_dsdsvwn
   case('dodblyp', 'dod-blyp')
      num = p_dodblyp
   case('dodpbe', 'dod-pbe')
      num = p_dodpbe
   case('dodpbeb95', 'dod-pbeb95')
      num = p_dodpbeb95
   case('dodpbep86', 'dod-pbep86')
      num = p_dodpbep86
   case('dodsvwn', 'dod-svwn')
      num = p_dodsvwn
   case('pbe02', 'pbe0-2')
      num = p_pbe0_2
   case('pbe0dh', 'pbe0-dh')
      num = p_pbe0_dh
   case('hf-3c', 'hf3c')
      num = p_hf3c
   case('hf-3cv', 'hf3cv')
      num = p_hf3cv
   case('pbeh3c', 'pbeh-3c')
      num = p_pbeh3c
   case('b973c', 'b97-3c')
      num = p_b973c
   case('hsesol')
      num = p_hsesol
   case('pwgga')
      num = p_pwgga
   case('dftb3', 'dftb(3ob)')
      num = p_dftb_3ob
   case('dftb(mio)')
      num = p_dftb_mio
   case('dftb(pbc)')
      num = p_dftb_pbc
   case('dftb(matsci)')
      num = p_dftb_matsci
   case('lc-dftb', 'dftb(ob2)')
      num = p_dftb_ob2
   case('hcth120')
      num = p_hcth120
   case('ptpss')
      num = p_ptpss
   case('lc-wpbe', 'lcwpbe')
      num = p_lcwpbe
   case('bmk')
      num = p_bmk
   case('b1b95')
      num = p_b1b95
   case('bwb6k')
      num = p_pwb6k
   case('otpss')
      num = p_otpss
   case('ssb')
      num = p_ssb
   case('revssb')
      num = p_revssb
   case('pbesol')
      num = p_pbesol
   case('hse06')
      num = p_hse06
   case('pbexalpha')
      num = p_pbexalpha
   case('pbehpbe')
      num = p_pbehpbe
   case('hcth407')
      num = p_hcth407
   case('n12')
      num = p_n12
   case('pkzb')
      num = p_pkzb
   case('thcth', 'tauhctc')
      num = p_thcth
   case('m11l')
      num = p_m11l
   case('mn15l')
      num = p_mn15l
   case('mpwb1k')
      num = p_mpwb1k
   case('mpw1kcis')
      num = p_mpw1kcis
   case('mpwkcis1k')
      num = p_mpwkcis1k
   case('pbeh1pbe')
      num = p_pbeh1pbe
   case('pbe1kcis')
      num = p_pbe1kcis
   case('b97-1')
      num = p_b97_1
   case('b97-2')
      num = p_b97_2
   case('b98')
      num = p_b98
   case('hiss')
      num = p_hiss
   case('hse03')
      num = p_hse03
   case('revtpssh')
      num = p_revtpssh
   case('tpss1kcis')
      num = p_tpss1kcis
   case('m05')
      num = p_m05
   case('m052x', 'm05-2x')
      num = p_m052x
   case('m08hx', 'm08-hx')
      num = p_m08hx
   case('lcwhpbe', 'lc-whpbe')
      num = p_lcwhpbe
   case('mn12l')
      num = p_mn12l
   case('tauhcthhyb')
      num = p_tauhcthhyb
   case('sogga11x')
      num = p_sogga11x
   case('n12sx')
      num = p_n12sx
   case('mn12sx')
      num = p_mn12sx
   case('mn15')
      num = p_mn15
   case('glyp', 'g-lyp')
      num = p_glyp
   case('revpbe0dh', 'revpbe0-dh')
      num = p_revpbe0dh
   case('revtpss0')
      num = p_revtpss0
   case('revdsd-pbep86', 'revdsdpbep86')
      num = p_revdsdpbep86
   case('revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe')
      num = p_revdsdpbe
   case('revdsd-blyp', 'revdsdblyp')
      num = p_revdsdblyp
   case('revdod-pbep86', 'revdodpbep86')
      num = p_revdodpbep86
   case('b97m')
      num = p_b97m
   case('wb97m', 'ωb97m', 'omegab97m')
      num = p_wb97m
   end select
end function get_functional_id

!> Convert string to lower case
pure function lowercase(str) result(lcstr)
   character(len=*), intent(in)  :: str
   character(len=len_trim(str)) :: lcstr
   integer :: ilen, ioffset, iquote, i, iav, iqc

   ilen=len_trim(str)
   ioffset=iachar('A')-iachar('a')
   iquote=0
   lcstr=str
   do i=1, ilen
      iav=iachar(str(i:i))
      if(iquote==0 .and. (iav==34 .or.iav==39)) then
         iquote=1
         iqc=iav
        cycle
      endif
      if(iquote==1 .and. iav==iqc) then
         iquote=0
         cycle
      endif
      if (iquote==1) cycle
      if(iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr(i:i)=achar(iav-ioffset)
      else
         lcstr(i:i)=str(i:i)
      endif
   enddo

end function lowercase


end module dftd4_param

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_output
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   use mctc_io_constants, only : pi
   use multicharge_model, only : mchrg_model_type
   implicit none
   private

   public :: write_ascii_model, write_ascii_properties, write_ascii_results

contains

subroutine write_ascii_model(unit, mol, model)

   !> Formatted unit
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   integer :: isp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)

   write(unit, '(a, ":")') "Charge model parameter"
   write(unit, '(54("-"))')
   write(unit, '(a4,5x,*(1x,a10))') "Z", "chi/Eh", "kcn/Eh", "eta/Eh", "rad/AA"
   write(unit, '(54("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4, *(1x,f10.4))') &
         & mol%num(isp), mol%sym(isp), model%chi(isp), model%kcn(isp), &
         & model%eta(isp) + sqrt2pi/model%rad(isp), model%rad(isp) * autoaa
   end do
   write(unit, '(54("-"),/)')

end subroutine write_ascii_model

subroutine write_ascii_properties(unit, mol, model, cn, qvec)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Coordination numbers
   real(wp), intent(in) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(in) :: qvec(:)

   integer :: iat, isp

   write(unit, '(a,":")') "Electrostatic properties (in atomic units)"
   write(unit, '(50("-"))')
   write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "CN", "q", "chi"
   write(unit, '(50("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6,1x,i4,1x,a4,*(1x,f10.4))') &
         & iat, mol%num(isp), mol%sym(isp), cn(iat), qvec(iat), &
         & model%chi(isp) - model%kcn(isp) * sqrt(cn(iat))
   end do
   write(unit, '(50("-"),/)')

end subroutine write_ascii_properties

subroutine write_ascii_results(unit, mol, energy, gradient, sigma)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   real(wp), intent(in) :: energy(:)
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)

   integer :: iat, isp
   logical :: grad
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]

   grad = present(gradient) .and. present(sigma)

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Electrostatic energy", sum(energy), "Eh"
   write(unit, '(a)')
   if (grad) then
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
         & "Gradient norm", norm2(gradient), "Eh/a0"
      write(unit, '(50("-"))')
      write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "dE/dx", "dE/dy", "dE/dz"
      write(unit, '(50("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         write(unit, '(i6,1x,i4,1x,a4,*(es11.3))') &
            & iat, mol%num(isp), mol%sym(isp), gradient(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')

      write(unit, '(a,":")') &
         & "Virial"
      write(unit, '(50("-"))')
      write(unit, '(a15,1x,*(1x,a10))') "component", "x", "y", "z"
      write(unit, '(50("-"))')
      do iat = 1, 3
         write(unit, '(2x,4x,1x,a4,1x,4x,*(es11.3))') &
            & comp(iat), sigma(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')
   end if

end subroutine write_ascii_results

end module multicharge_output

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge_model, only : mchrg_model_type, new_mchrg_model
   use multicharge_param_eeq2019, only : get_eeq_chi, get_eeq_eta, &
      & get_eeq_rad, get_eeq_kcn
   implicit none
   private

   public :: new_eeq2019_model

contains

subroutine new_eeq2019_model(mol, model)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Electronegativity equilibration model
   type(mchrg_model_type), intent(out) :: model

   real(wp), allocatable :: chi(:), eta(:), kcn(:), rad(:)

   chi = get_eeq_chi(mol%num)
   eta = get_eeq_eta(mol%num)
   kcn = get_eeq_kcn(mol%num)
   rad = get_eeq_rad(mol%num)

   call new_mchrg_model(model, chi=chi, rad=rad, eta=eta, kcn=kcn)

end subroutine new_eeq2019_model

end module multicharge_param

! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_data, only : get_covalent_rad
   use multicharge_model, only : mchrg_model_type
   use multicharge_ncoord, only : get_coordination_number
   use multicharge_output, only : write_ascii_model, write_ascii_properties, &
      & write_ascii_results
   use multicharge_param, only : new_eeq2019_model
   use multicharge_version, only : get_multicharge_version
   implicit none
   public


end module multicharge

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Interface to the charge model
module dftd4_charge
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge, only : mchrg_model_type, new_eeq2019_model, &
      & write_ascii_model, write_ascii_properties, write_ascii_results, &
      & get_coordination_number, get_covalent_rad, get_lattice_points
   implicit none
   private

   public :: get_charges


contains


!> Obtain charges from electronegativity equilibration model
subroutine get_charges(mol, qvec, dqdr, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   logical :: grad
   type(mchrg_model_type) :: model
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: rcov(:), trans(:, :)

   grad = present(dqdr) .and. present(dqdL)

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   end if

   rcov = get_covalent_rad(mol%num)
   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   call model%solve(mol, cn, dcndr, dcndL, qvec=qvec, dqdr=dqdr, dqdL=dqdL)

end subroutine get_charges


end module dftd4_charge

! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> High-level wrapper to obtain the dispersion energy for a DFT-D4 calculation
module dftd4_disp
   use dftd4_blas, only : d4_gemv
   use dftd4_charge, only : get_charges
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd4_damping, only : damping_param
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : d4_model
   use dftd4_ncoord, only : get_coordination_number
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   implicit none
   private

   public :: get_dispersion, get_properties, get_pairwise_dispersion


contains


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   allocate(cn(mol%nat))
   if (grad) allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, &
      & cn, dcndr, dcndL)

   allocate(q(mol%nat))
   if (grad) allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
   call get_charges(mol, q, dqdr, dqdL)

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat), gwdq(mref, mol%nat))
   call disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   allocate(energies(mol%nat))
   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
      call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call d4_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      call d4_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

   energy = sum(energies)

end subroutine get_dispersion


!> Wrapper to handle the evaluation of properties related to this dispersion model
subroutine get_properties(mol, disp, cutoff, cn, q, c6, alpha)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Coordination number
   real(wp), intent(out) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(out) :: q(:)

   !> C6 coefficients
   real(wp), intent(out) :: c6(:, :)

   !> Static polarizibilities
   real(wp), intent(out) :: alpha(:)

   integer :: mref
   real(wp), allocatable :: gwvec(:, :), lattr(:, :)

   mref = maxval(disp%ref)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   call get_charges(mol, q)

   allocate(gwvec(mref, mol%nat))
   call disp%weight_references(mol, cn, q, gwvec)

   call disp%get_atomic_c6(mol, gwvec, c6=c6)
   call disp%get_polarizibilities(mol, gwvec, alpha=alpha)

end subroutine get_properties


!> Wrapper to handle the evaluation of pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion(mol, disp, param, cutoff, energy2, energy3)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   integer :: mref
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :), c6(:, :), lattr(:, :)

   mref = maxval(disp%ref)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   call get_charges(mol, q)

   allocate(gwvec(mref, mol%nat))
   call disp%weight_references(mol, cn, q, gwvec)

   allocate(c6(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   energy2(:, :) = 0.0_wp
   energy3(:, :) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_pairwise_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, &
      & c6, energy2)

   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_pairwise_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
      & c6, energy3)

end subroutine get_pairwise_dispersion


end module dftd4_disp

module dftd4
   use mctc_env, only : wp
   use mctc_io, only : structure_type, new
   use dftd4_blas, only : d4_gemv
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd4_data, only : get_covalent_rad
   use dftd4_ncoord, only : get_coordination_number
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_model, only : d4_model, new_d4_model
   use dftd4_param, only : get_rational_damping
   use dftd4_version, only : get_dftd4_version
   implicit none
   private

   public :: dftd4_dispersion

contains


!> Evaluate DFT-D4 dispersion energy
subroutine dftd4_dispersion(num, xyz, charges, s6, s8, s9, a1, a2, ga, gc, energy)

   !> Atomic number for each atom
   integer, intent(in) :: num(:)

   !> Cartesian coordinates for each atom
   real(wp), intent(in) :: xyz(:, :)

   !> Atomic partial charges for each atom
   real(wp), contiguous, intent(in) :: charges(:)

   !> Scaling parameter for dipole-dipole terms
   real(wp), intent(in) :: s6

   !> Scaling parameter for quadrupole-dipole terms
   real(wp), intent(in) :: s8

   !> Scaling parameter for tripole dipole terms
   real(wp), intent(in) :: s9

   !> Scaling parameter for critical radii
   real(wp), intent(in) :: a1

   !> Offset parameter for critical radii
   real(wp), intent(in) :: a2

   !> Charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Dispersion energy
   real(wp), intent(out) :: energy

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d4_model) :: disp

   ! Create new molecular structure
   call new(mol, num, xyz, charge=sum(charges))

   ! Create dispersion model
   call new_d4_model(disp, mol, ga=ga, gc=gc)

   ! Create damping parameters
   param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2)

   ! Perform actual calculation of dispersion energy
   call get_dispersion(mol, disp, param, realspace_cutoff(), charges, energy)
end subroutine dftd4_dispersion


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion(mol, disp, param, cutoff, charges, energy, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: charges(:)

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   allocate(cn(mol%nat))
   if (grad) allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, &
      & cn, dcndr, dcndL)

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat), gwdq(mref, mol%nat))
   call disp%weight_references(mol, cn, charges, gwvec, gwdcn, gwdq)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   allocate(energies(mol%nat))
   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)

   allocate(q(mol%nat))
   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call d4_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      call d4_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

   energy = sum(energies)

end subroutine get_dispersion

end module dftd4
