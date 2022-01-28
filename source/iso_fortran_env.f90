MODULE ISO_FORTRAN_ENV
    INTEGER (KIND=4), PARAMETER                :: ATOMIC_INT_KIND              = 4
    INTEGER (KIND=4), PARAMETER                :: ATOMIC_LOGICAL_KIND          = 4
    INTEGER (KIND=4), PARAMETER, DIMENSION (1) :: CHARACTER_KINDS              = (/ 1 /)
    INTEGER (KIND=4), PARAMETER                :: CHARACTER_STORAGE_SIZE       = 8
    INTEGER (KIND=4), PARAMETER                :: ERROR_UNIT                   = 0
    INTEGER (KIND=4), PARAMETER                :: FILE_STORAGE_SIZE            = 8 ! -assume byterecl
    INTEGER (KIND=4), PARAMETER                :: INPUT_UNIT                   = 5
    INTEGER (KIND=4), PARAMETER, DIMENSION (4) :: INTEGER_KINDS                = (/1, 2, 4, 8/)
    INTEGER (KIND=4), PARAMETER                :: INT8                         = 1
    INTEGER (KIND=4), PARAMETER                :: INT16                        = 2
    INTEGER (KIND=4), PARAMETER                :: INT32                        = 4
    INTEGER (KIND=4), PARAMETER                :: INT64                        = 8
    INTEGER (KIND=4), PARAMETER                :: IOSTAT_END                   = -1
    INTEGER (KIND=4), PARAMETER                :: IOSTAT_EOR                   = -2
    INTEGER (KIND=4), PARAMETER, DIMENSION (4) :: LOGICAL_KINDS                = (/1, 2, 4, 8/)
    INTEGER (KIND=4), PARAMETER                :: NUMERIC_STORAGE_SIZE         = 32
    INTEGER (KIND=4), PARAMETER                :: OUTPUT_UNIT                  = 6
    INTEGER (KIND=4), PARAMETER, DIMENSION (3) :: REAL_KINDS                   = (/4, 8, 16/)
    INTEGER (KIND=4), PARAMETER                :: REAL32                       = 4
    INTEGER (KIND=4), PARAMETER                :: REAL64                       = 8
    INTEGER (KIND=4), PARAMETER                :: REAL128                      = 16
!
! These values are the error-message indices for the related error messages.
!
    INTEGER (KIND=4), PARAMETER                :: STAT_LOCKED                  = 775
    INTEGER (KIND=4), PARAMETER                :: STAT_UNLOCKED                = 776
    INTEGER (KIND=4), PARAMETER                :: STAT_LOCKED_OTHER_IMAGE      = 777
    INTEGER (KIND=4), PARAMETER                :: STAT_STOPPED_IMAGE           = 778

! -------------------------------------------------------------------------
! LOCK_TYPE definition
! It's a special derived type for LOCK / UNLOCK statements
! -------------------------------------------------------------------------

END MODULE ISO_FORTRAN_ENV
