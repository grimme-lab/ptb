module molcom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      integer  n

      integer, allocatable :: at(:)  
      real(wp),allocatable :: z(:)  
      real(wp),allocatable :: xyz(:,:)  
      real(wp),allocatable :: rab(:)  

end module molcom       
