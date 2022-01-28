module mocom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp),allocatable :: cmo_ref(:,:)  
      real(wp),allocatable :: epsref (:)  
      real(wp),allocatable :: totmatch    
      real(wp)             :: fitdat(10000)
      integer              :: fitcount         

end module mocom       
