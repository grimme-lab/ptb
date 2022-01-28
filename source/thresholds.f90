module thresholds
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp),parameter :: gamma_thr = 1.d-7

end module thresholds   
